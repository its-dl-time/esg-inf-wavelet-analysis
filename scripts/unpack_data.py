#!/usr/bin/env python3
"""Expand the frozen CSV snapshots after cloning (Python standard library only)."""
import gzip
import pathlib
import shutil

root = pathlib.Path(__file__).resolve().parents[1]
for path in sorted((root/'data').rglob('*.csv.gz')):
    target = path.with_suffix('')
    if target.exists():
        print('Already present:',target.relative_to(root))
        continue
    with gzip.open(path,'rb') as source, target.open('wb') as destination:
        shutil.copyfileobj(source,destination)
    print('Expanded:',target.relative_to(root))
