#!/usr/bin/env bash
# Ejecuta todas las simulaciones de la revision y agrega los resultados.
# Uso: bash run_all.sh            (presets por defecto)
#      NBATCH=16 bash run_all.sh  (duplica el numero de realizaciones)
set -e
cd "$(dirname "$0")"
EXTRA=""
if [ -n "$NBATCH" ]; then EXTRA="--batches $NBATCH"; fi
python run_sde.py --preset smoke
python run_sde.py --preset tau20 $EXTRA
python run_sde.py --preset long200 $EXTRA
python run_sde.py --preset long1000 $EXTRA
python aggregate.py --presets tau20 long200 long1000
