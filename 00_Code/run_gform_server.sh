#!/usr/bin/env bash
# Lanzar G-Formula en servidor Linux (paralelo, cohorte completa).
#
# Uso (desde la raíz del proyecto):
#   bash "00_Code/run_gform_server.sh"          # inicia en tmux (sobrevive a SSH)
#   bash "00_Code/run_gform_server.sh" --attach # reconectar a la sesión
#   bash "00_Code/run_gform_server.sh" --inside # solo uso interno (dentro de tmux)
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
  export R_FUTURE_FORK_ENABLE=1

  mkdir -p "02_Output/G-Form/Timing"
  LOG="02_Output/G-Form/Timing/server_run_$(date +%Y%m%d_%H%M%S).log"
  echo "$LOG" > "02_Output/G-Form/Timing/batch_run.logpath"

  echo "=== G-Formula servidor ==="
  echo "Directorio: $ROOT"
  echo "Sesión tmux: $SESSION"
  echo "Log: $LOG"
  echo "CPUs: $(nproc 2>/dev/null || echo '?')"
  echo "Inicio: $(date)"

  Rscript "00_Code/10.2 G-Form_models.R" 2>&1 | tee "$LOG"
  EXIT=${PIPESTATUS[0]}
  echo "EXIT=$EXIT" >> "$LOG"
  echo "Fin: $(date) | EXIT=$EXIT"

  if [ "$EXIT" -ne 0 ]; then
    echo "La corrida terminó con error (EXIT=$EXIT)."
  else
    echo "La corrida terminó correctamente."
  fi

  return "$EXIT"
}

case "${1:-}" in
  --inside)
    run_gform_job
    exit $?
    ;;
  --attach)
    if ! command -v tmux >/dev/null 2>&1; then
      echo "tmux no está instalado." >&2
      exit 1
    fi
    if ! tmux has-session -t "$SESSION" 2>/dev/null; then
      echo "No existe la sesión tmux '$SESSION'." >&2
      echo "Iniciar con: bash $(printf '%q' "$SCRIPT")" >&2
      exit 1
    fi
    exec tmux attach -t "$SESSION"
    ;;
  --help|-h)
    cat <<EOF
Uso: bash $(basename "$SCRIPT") [opción]

  (sin opción)   Inicia la corrida en tmux (desacoplada; sobrevive al cierre de SSH)
  --attach       Reconecta a la sesión tmux '$SESSION'
  --inside       Ejecuta la corrida en el terminal actual (uso interno)
  --help         Muestra esta ayuda

Variables:
  GFORM_TMUX_SESSION   Nombre de la sesión tmux (default: gform-server)

Monitoreo del log (sin tmux):
  tail -f "\$(cat 02_Output/G-Form/Timing/batch_run.logpath)"
EOF
    exit 0
    ;;
  "")
    if ! command -v tmux >/dev/null 2>&1; then
      echo "tmux no encontrado; ejecutando en foreground (se detiene al cerrar SSH)." >&2
      run_gform_job
      exit $?
    fi

    if tmux has-session -t "$SESSION" 2>/dev/null; then
      echo "La sesión tmux '$SESSION' ya está activa." >&2
      echo "  Reconectar: bash $(printf '%q' "$SCRIPT") --attach" >&2
      echo "  Log:        tail -f \"\$(cat $ROOT/02_Output/G-Form/Timing/batch_run.logpath)\"" >&2
      exit 1
    fi

    tmux new-session -d -s "$SESSION" -c "$ROOT" \
      "bash $(printf '%q' "$SCRIPT") --inside; echo; echo 'Pulsa Enter para cerrar la ventana tmux...'; read -r"

    echo "Corrida iniciada en tmux (sesión: $SESSION)."
    echo "  Reconectar:  bash $(printf '%q' "$SCRIPT") --attach"
    echo "  Ver log:     tail -f \"\$(cat $ROOT/02_Output/G-Form/Timing/batch_run.logpath)\""
    echo "  Listar:      tmux ls"
    exit 0
    ;;
  *)
    echo "Opción desconocida: $1" >&2
    echo "Use --help para ver las opciones." >&2
    exit 1
    ;;
esac
