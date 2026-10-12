#!/usr/bin/env python3
"""Combine separately executed alpha-tau audits without changing their output values."""
from pathlib import Path
import csv
O=Path(__file__).resolve().parent/'revision_outputs'
for key in ['summary','states']:
    rows=[]
    for a in (2,4,6):
        for t in (10,20,40):
            p=O/f'ecm_audit_{key}_a{a}_t{t}.csv'
            if not p.is_file(): raise SystemExit(f'Missing {p}; run ECM audit for this case')
            with p.open() as f:rows.extend(csv.DictReader(f))
    out=O/f'ecm_audit_{key}.csv'
    with out.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print(out.name,len(rows))
