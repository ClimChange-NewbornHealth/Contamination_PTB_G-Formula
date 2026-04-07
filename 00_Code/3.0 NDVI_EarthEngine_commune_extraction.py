"""
Extract NDVI (MODIS MOD13Q1, ~16-day composites) per district geometry.
Builds a daily calendar from 2010-01-01 through 2020-12-31; days without
a composite are NA.

Requires: shapefile at 01_Data/Input/district_geo/district_geo.shp (WGS84).

    pip uninstall ee  # remove wrong package if needed
    pip install earthengine-api geopandas pandas

First run: ee.Authenticate() then ee.Initialize(project=...).
"""

import os

import ee
import geopandas as gpd
import pandas as pd

PROJECT_ID = "quadrant-rm"
SHAPE_DIR = "01_Data/Input/district_geo/"
OUT_DIR = "01_Data/Input/Clime_series/"
OUT_CSV = os.path.join(OUT_DIR, "ndvi_daily_district.csv")

# EE end is exclusive
START_DATE = "2010-01-01"
END_DATE_EE = "2021-01-01"  # includes all of 2020-12-31
DAILY_START = "2010-01-01"
DAILY_END = "2020-12-31"

try:
    ee.Initialize(project=PROJECT_ID)
except Exception:
    ee.Authenticate()
    ee.Initialize(project=PROJECT_ID)

os.makedirs(OUT_DIR, exist_ok=True)

shapefile_path = os.path.join(SHAPE_DIR, "district_geo.shp")
if not os.path.exists(shapefile_path):
    raise FileNotFoundError(
        f"Shapefile not found: {shapefile_path}\n"
        "Place district_geo.shp (and sidecar files) in that folder."
    )

commune_gdf = gpd.read_file(shapefile_path)
print(f"Loaded {len(commune_gdf)} district geometries")

if commune_gdf.crs is None or str(commune_gdf.crs) != "EPSG:4326":
    print(f"Reprojecting to EPSG:4326 (was {commune_gdf.crs})")
    commune_gdf = commune_gdf.to_crs("EPSG:4326")

if "geometry_id" not in commune_gdf.columns:
    commune_gdf["geometry_id"] = range(len(commune_gdf))

if "codigo_comuna" not in commune_gdf.columns:
    commune_gdf["codigo_comuna"] = commune_gdf["geometry_id"]


def gdf_to_ee_featurecollection(gdf):
    features = []
    for idx, row in gdf.iterrows():
        geom = row.geometry
        try:
            if geom.geom_type == "Polygon":
                coords = [[[p[0], p[1]] for p in geom.exterior.coords]]
                ee_geom = ee.Geometry.Polygon(coords)
            elif geom.geom_type == "MultiPolygon":
                coords = [
                    [[[p[0], p[1]] for p in poly.exterior.coords] for poly in geom.geoms]
                ]
                ee_geom = ee.Geometry.MultiPolygon(coords)
            else:
                print(f"Skipping geometry {idx}: type {geom.geom_type}")
                continue

            geometry_id = int(row.get("geometry_id", idx))
            codigo_comuna = int(row.get("codigo_comuna", geometry_id))

            features.append(
                ee.Feature(
                    ee_geom,
                    {"geometry_id": geometry_id, "codigo_comuna": codigo_comuna},
                )
            )
        except Exception as e:
            print(f"Error processing geometry {idx}: {e}")
            continue

    return ee.FeatureCollection(features)


def extract_ndvi_sparse(feature_collection, start_date, end_date_exclusive):
    """One row per geometry per MOD13Q1 composite date (16-day product).

    NDVI is kept as stored in the product (no scale factor); MOD13Q1 NDVI is
    typically in raw integer units (e.g. −2000 to 10000).
    """
    modis = (
        ee.ImageCollection("MODIS/061/MOD13Q1")
        .filterDate(start_date, end_date_exclusive)
        .select(["NDVI"])
    )

    n_images = modis.size().getInfo()
    print(f"Found {n_images} MODIS MOD13Q1 composites (16-day)")

    def process_image(image):
        date = ee.Date(image.get("system:time_start"))
        date_str = date.format("YYYY-MM-dd")
        ndvi = image.select("NDVI")

        def extract_stats(feature):
            geom = feature.geometry()
            stats = ndvi.reduceRegion(
                reducer=ee.Reducer.mean(),
                geometry=geom,
                scale=250,
                maxPixels=1e9,
            )
            ndvi_value = stats.get("NDVI")
            return ee.Feature(
                None,
                {
                    "geometry_id": feature.get("geometry_id"),
                    "codigo_comuna": feature.get("codigo_comuna"),
                    "date": date_str,
                    "ndvi": ndvi_value,
                },
            )

        return feature_collection.map(extract_stats)

    return modis.map(process_image).flatten()


def feature_collection_to_dataframe(fc, page_size=5000):
    """Download FeatureCollection to pandas (paginated)."""
    n = fc.size().getInfo()
    print(f"Downloading {n} features from Earth Engine...")
    rows = []
    for start in range(0, n, page_size):
        sub = ee.FeatureCollection(fc.toList(page_size, start))
        payload = sub.getInfo()
        for f in payload["features"]:
            p = f["properties"]
            rows.append(
                {
                    "geometry_id": p.get("geometry_id"),
                    "codigo_comuna": p.get("codigo_comuna"),
                    "date": p.get("date"),
                    "ndvi": p.get("ndvi"),
                }
            )
    return pd.DataFrame(rows)


def expand_to_daily_calendar(df_sparse, daily_start, daily_end):
    """Full daily index per geometry_id; NDVI NA between composite dates."""
    if df_sparse.empty:
        return df_sparse

    df_sparse = df_sparse.copy()
    df_sparse["date"] = pd.to_datetime(df_sparse["date"], errors="coerce")
    df_sparse = df_sparse.dropna(subset=["date"])

    meta = (
        df_sparse[["geometry_id", "codigo_comuna"]]
        .drop_duplicates(subset=["geometry_id"])
        .set_index("geometry_id")["codigo_comuna"]
    )

    all_dates = pd.date_range(daily_start, daily_end, freq="D")
    geometry_ids = sorted(df_sparse["geometry_id"].unique())

    idx = pd.MultiIndex.from_product(
        [geometry_ids, all_dates], names=["geometry_id", "date"]
    )
    wide = (
        df_sparse.set_index(["geometry_id", "date"])["ndvi"]
        .reindex(idx)
        .rename("ndvi")
        .reset_index()
    )
    wide["codigo_comuna"] = wide["geometry_id"].map(meta)
    cols = ["geometry_id", "codigo_comuna", "date", "ndvi"]
    return wide[cols]


# --- run ---
ee_fc = gdf_to_ee_featurecollection(commune_gdf)
print(f"Earth Engine FeatureCollection size: {ee_fc.size().getInfo()}")

print("\n=== Extracting NDVI (MOD13Q1 composites) ===")
ndvi_fc = extract_ndvi_sparse(ee_fc, START_DATE, END_DATE_EE)

df_sparse = feature_collection_to_dataframe(ndvi_fc)
print(f"Sparse rows: {len(df_sparse)}")

df_daily = expand_to_daily_calendar(df_sparse, DAILY_START, DAILY_END)
df_daily["date"] = df_daily["date"].dt.strftime("%Y-%m-%d")

df_daily.to_csv(OUT_CSV, index=False)
print(f"Wrote {len(df_daily)} rows to {OUT_CSV}")
