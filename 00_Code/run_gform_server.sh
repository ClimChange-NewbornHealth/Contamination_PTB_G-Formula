#!/usr/bin/env bash
# Lanzar G-Formula en servidor Linux (paralelo, cohorte completa).
#
# Uso (desde la raíz del proyecto):
#   bash "00_Code/run_gform_server.sh"           # primer plano (default; ves la salida en SSH)
#   bash "00_Code/run_gform_server.sh" --tmux     # segundo plano con tmux (sobrevive a SSH)
#   bash "00_Code/run_gform_server.sh" --attach  # reconectar a sesión tmux
#
# Sesión tmux por defecto: gform-server (cambiar con GFORM_TMUX_SESSION=nombre)

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
SESSION="${GFORM_TMUX_SESSION:-gform-server}"

run_gform_job() {
  cd "$ROOT"

  export GFORM_EXEC_MODE=server
  export OMP_NUM_THREADS=1
  export OPENBLAS_NUM_THREADS=1
  export MKL_NUM_THREADS=1
  export VECLIB_MAXIMUM_THREADS=1
  unset R_PARALLELLY_FORK_ENABLE 2>/dev/null || true
  export R_FUTURE_FORK_ENABLE=true
  export GFORM_RESERVE_RAM_GB="${GFORM_RESERVE_RAM_GB:-56}"
  export GFORM_BOOTSTRAP_WORKERS="${GFORM_BOOTSTRAP_WORKERS:-1}"
  export GFORM_BOOTSTRAP_MAX_WORKERS="${GFORM_BOOTSTRAP_MAX_WORKERS:-8}"
  export GFORM_BOOTSTRAP_RAM_PER_WORKER_GB="${GFORM_BOOTSTRAP_RAM_PER_WORKER_GB:-8}"
  export GFORM_BOOTSTRAP_PARALLEL="${GFORM_BOOTSTRAP_PARALLEL:-false}"
  export GFORM_HEATMAP_PARALLEL="${GFORM_HEATMAP_PARALLEL:-false}"
  export GFORM_HEATMAP_WORKERS="${GFORM_HEATMAP_WORKERS:-4}"
  export GFORM_HEATMAP_BATCH_SIZE="${GFORM_HEATMAP_BATCH_SIZE:-2}"
  export GFORM_HEATMAP_MAX_WORKERS="${GFORM_HEATMAP_MAX_WORKERS:-4}"
  export GFORM_HEATMAP_RAM_PER_WORKER_GB="${GFORM_HEATMAP_RAM_PER_WORKER_GB:-14}"
  export GFORM_BOOT_ITER="${GFORM_BOOT_ITER:-250}"
  export GFORM_GLOBALS_MAX_GB="${GFORM_GLOBALS_MAX_GB:-96}"
  export GFORM_SKIP_COMPLETED="${GFORM_SKIP_COMPLETED:-true}"
  export GFORM_RUN_BOOTSTRAP="${GFORM_RUN_BOOTSTRAP:-true}"
  export GFORM_MAX_BATCH_HOURS="${GFORM_MAX_BATCH_HOURS:-168}"

  mkdir -p "02_Output/G-Form/Timing"
  LOG="02_Output/G-Form/Timing/server_run_$(date +%Y%m%d_%H%M%S).log"
  echo "$LOG" > "02_Output/G-Form/Timing/batch_run.logpath"

  echo "=== G-Formula servidor ==="
  echo "Directorio: $ROOT"
  echo "Modo: primer plano (salida en esta terminal)"
  echo "Log: $LOG"
  echo "CPUs: $(nproc 2>/dev/null || echo '?')"
  echo "Inicio: $(date)"
  echo "Nota: al cerrar SSH se detiene la corrida. Usa --tmux para dejarla en segundo plano."
  echo ""

  Rscript "00_Code/10.2 G-Form_models.R" 2>&1 | tee "$LOG"
  EXIT=${PIPESTATUS[0]}
  echo "EXIT=$EXIT" >> "$LOG"
  echo "Fin: $(date) | EXIT=$EXIT"

  if [ "$EXIT" -ne 0 ]; then
    echo "La corrida terminó con error (EXIT=$EXIT)."
  else
    echo "La corrida terminó correctamente."
    echo "Descargar outputs al Mac: bash \"00_Code/sync_from_server.sh\""
  fi

  return "$EXIT"
}

launch_tmux() {
  if ! command -v tmux >/dev/null 2>&1; then
    echo "tmux no está instalado. Ejecutando en primer plano." >&2
    run_gform_job
    return $?
  fi

  if tmux has-session -t "$SESSION" 2>/dev/null; then
    echo "La sesión tmux '$SESSION' ya está activa." >&2
    echo "  Reconectar: bash $(printf '%q' "$SCRIPT") --attach" >&2
    echo "  Log:        tail -f \"\$(cat $ROOT/02_Output/G-Form/Timing/batch_run.logpath)\"" >&2
    return 1
  fi

  tmux new-session -d -s "$SESSION" -c "$ROOT" \
    "bash $(printf '%q' "$SCRIPT"); echo; echo 'Pulsa Enter para cerrar la ventana tmux...'; read -r"

  echo "Corrida iniciada en tmux (sesión: $SESSION)."
  echo "  Reconectar:  bash $(printf '%q' "$SCRIPT") --attach"
  echo "  Ver log:     tail -f \"\$(cat $ROOT/02_Output/G-Form/Timing/batch_run.logpath)\""
  echo "  Listar:      tmux ls"
}

case "${1:-}" in
  --tmux)
    launch_tmux
    ;;
  --attach)
    if ! command -v tmux >/dev/null 2>&1; then
      echo "tmux no está instalado." >&2
      exit 1
    fi
    if ! tmux has-session -t "$SESSION" 2>/dev/null; then
      echo "No existe la sesión tmux '$SESSION'." >&2
      echo "Iniciar con: bash $(printf '%q' "$SCRIPT") --tmux" >&2
      exit 1
    fi
    exec tmux attach -t "$SESSION"
    ;;
  --help|-h)
    cat <<EOF
Uso: bash $(basename "$SCRIPT") [opción]

  (sin opción)   Ejecuta en primer plano (ves la salida en SSH; se detiene al cerrar SSH)
  --tmux         Ejecuta en segundo plano con tmux (sobrevive al cierre de SSH)
  --attach       Reconecta a la sesión tmux '$SESSION'
  --help         Muestra esta ayuda

Variables:
  GFORM_TMUX_SESSION              Nombre de la sesión tmux (default: gform-server)
  GFORM_OUTPUT_SUBDIR             Subcarpeta bajo 02_Output/G-Form/ (vacío = full)
  GFORM_SAMPLE_FRAC               Submuestra 0-1 (ej. 0.01 = 1%)
  GFORM_INTERVENTIONS             IDs a correr (ej. 1 o 1,2,3)
  GFORM_BOOT_ITER                 Réplicas bootstrap (default: 250)
  GFORM_BOOTSTRAP_PARALLEL        true/false (default: false = secuencial)
  GFORM_BOOTSTRAP_WORKERS         Workers si paralelo (default: 1)
  GFORM_HEATMAP_PARALLEL          true/false (default: false = secuencial)
  GFORM_HEATMAP_WORKERS           Workers heatmap si paralelo (default: 4)
  GFORM_HEATMAP_BATCH_SIZE        Columnas por lote heatmap (default: 2)
  GFORM_HEATMAP_MAX_WORKERS       Tope workers heatmap (default: 4)
  GFORM_HEATMAP_RAM_PER_WORKER_GB RAM estimada por worker heatmap (default: 14)
  GFORM_RESERVE_RAM_GB            RAM reservada proceso padre (default: 56)
  GFORM_SKIP_COMPLETED            Omitir intervenciones ya completas (default: true)
  GFORM_GLOBALS_MAX_GB            Límite future.globals.maxSize (default: 96)

Descargar outputs al Mac:
  bash "00_Code/sync_from_server.sh"
  bash "00_Code/sync_from_server.sh" --stub pm25_krg_pct20

Ver log en otra terminal:
  tail -f "\$(cat 02_Output/G-Form/Timing/batch_run.logpath)"
EOF
    exit 0
    ;;
  ""|--foreground|-f)
    run_gform_job
    exit $?
    ;;
  *)
    echo "Opción desconocida: $1" >&2
    echo "Use --help para ver las opciones." >&2
    exit 1
    ;;
esac
