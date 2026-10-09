#!/usr/bin/env python3
"""Check the frozen 18-country panel without requiring R."""
import csv
import hashlib
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
manifest = {}
with (ROOT / 'docs/data_manifest.csv').open(newline='') as f:
    for row in csv.DictReader(f):
        manifest[row['series']] = row
assert len(manifest) == 18
data = {}
with (ROOT / 'data/processed/core_18.csv').open(newline='') as f:
    for row in csv.DictReader(f):
        data.setdefault(row['series'], []).append(row)
assert set(data) == set(manifest)
for name, rows in data.items():
    months = [x['month'] for x in rows]
    assert months == sorted(months) and len(months) == len(set(months))
    assert len(rows) == int(manifest[name]['rows'])
    assert months[0] == manifest[name]['first'] and months[-1] == manifest[name]['last']
for filename in ('panel.csv', 'cpi_raw.csv'):
    path = ROOT / 'data/raw' / filename
    assert path.exists() and path.stat().st_size > 100_000
print('Verified 18 countries,',sum(map(len,data.values())),'processed rows')
expected = {
    'data/processed/core_18.csv':'98239bc9c8100840e8e41c3f4e93461347b15f9ebbd530d1a7cc4b002f4d62c0',
    'data/raw/panel.csv':'d50770b5a7f7c202523dd58aaae622f306ab71d7dd47b60f7fdc8f7a85c05753',
    'data/raw/cpi_raw.csv':'892734cca67865e4d989a9b283511a06d137d3b6e3277e8ccea0a676da526e23',
}
for filename, digest in expected.items():
    path = ROOT / filename
    observed = hashlib.sha256(path.read_bytes()).hexdigest()
    assert observed == digest, f'Checksum mismatch: {filename}'
    print(filename,observed)
