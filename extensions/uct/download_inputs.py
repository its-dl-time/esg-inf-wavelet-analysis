"""Fetch versioned-at-run-time ancillary series; requires internet on the user's machine.

UCT: Rogers, Sun & Sun via policyuncertainty.com (CC BY 4.0).
Brent: EIA via FRED MCOILBRENTEU. See upstream use and attribution terms.
"""
import csv
import io
import pathlib
import urllib.request

HERE = pathlib.Path(__file__).resolve().parent
DATA = HERE / "data"
DATA.mkdir(exist_ok=True)
SOURCES = {
    "UCT.csv": "https://www.policyuncertainty.com/media/UCT.csv",
    "MCOILBRENTEU.csv": "https://fred.stlouisfed.org/graph/fredgraph.csv?id=MCOILBRENTEU",
}
for filename, url in SOURCES.items():
    with urllib.request.urlopen(url, timeout=30) as response:
        payload = response.read()
    text = payload.decode("utf-8-sig")
    rows = list(csv.reader(io.StringIO(text)))
    if len(rows) < 5:
        raise ValueError(f"Unexpected short source file: {url}")
    (DATA / filename).write_bytes(payload)
    print(filename, len(rows), "rows", "header:", rows[0])
print("Source snapshots downloaded. Record the retrieval date and raw file hashes in any publication.")
