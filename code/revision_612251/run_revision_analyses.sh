#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python rev_firstorder.py
python rev_windows.py
python rev_tolerance.py
python rev_mechanism.py
python rev_longcavity_strong.py
# Run ECM scans separately to constrain memory and isolate failures; safe to rerun.
for a in 2 4 6; do
 for t in 10 20 40; do
  python rev_ecm_audit.py --single "$a" "$t"
 done
done
python consolidate_ecm.py
python rev_sde.py --results results
python revision_qa.py
if [[ -f ../Postrun_tools/check_results.py ]]; then
 python ../Postrun_tools/check_results.py
else
 echo 'External raw result structural checker not in this deployment; rev_sde.py manifest completed.'
fi
printf '\nSUCCESS: all scientific revision checks passed.\n'
