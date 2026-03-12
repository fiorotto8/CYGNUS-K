#!/usr/bin/env python3
"""
Download tabulated solar-neutrino source spectra from the Bahcall SNdata archive.

Downloads:
- pp
- 8B
- hep
- 13N
- 15O
- 17F
- 7Be line shape

Notes:
- pep is not included here because it is a line source, not a standard continuous dPhi/dE table.
- 7Be is downloaded as a line-shape file.
"""

from pathlib import Path
import csv
import requests

OUTDIR = Path("solar_neutrino_tables")
OUTDIR.mkdir(exist_ok=True)

FILES = {
    "pp": "https://www.sns.ias.edu/~jnb/SNdata/Export/PPenergyspectrum/ppenergytab",
    "b8": "https://www.sns.ias.edu/~jnb/SNdata/Export/B8spectrum/b8spectrum.txt",
    "hep": "https://www.sns.ias.edu/~jnb/SNdata/Export/Hepspectrum/hepspectrum.dat",
    "n13": "https://www.sns.ias.edu/~jnb/SNdata/Export/CNOspectra/n13.dat",
    "o15": "https://www.sns.ias.edu/~jnb/SNdata/Export/CNOspectra/o15.dat",
    "f17": "https://www.sns.ias.edu/~jnb/SNdata/Export/CNOspectra/f17.dat",
    #"be7_lineshape": "https://www.sns.ias.edu/~jnb/SNdata/Export/7Belineshape/7belineshape.dat",
}


def download_file(name: str, url: str, timeout: int = 60) -> Path:
    response = requests.get(url, timeout=timeout)
    response.raise_for_status()

    # Keep a simple extension policy
    if url.endswith(".txt"):
        suffix = ".txt"
    elif url.endswith(".dat"):
        suffix = ".dat"
    else:
        suffix = ".dat"

    outpath = OUTDIR / f"{name}{suffix}"
    outpath.write_bytes(response.content)
    return outpath


def main():
    manifest_rows = []

    for name, url in FILES.items():
        try:
            path = download_file(name, url)
            size = path.stat().st_size
            manifest_rows.append(
                {
                    "source": name,
                    "filename": path.name,
                    "bytes": size,
                    "status": "ok",
                    "url": url,
                }
            )
            print(f"[OK] {name:12s} -> {path}")
        except Exception as exc:
            manifest_rows.append(
                {
                    "source": name,
                    "filename": "",
                    "bytes": 0,
                    "status": f"failed: {exc}",
                    "url": url,
                }
            )
            print(f"[FAIL] {name:12s} -> {exc}")

    manifest = OUTDIR / "manifest.csv"
    with manifest.open("w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["source", "filename", "bytes", "status", "url"]
        )
        writer.writeheader()
        writer.writerows(manifest_rows)

    print(f"\nWrote manifest: {manifest}")


if __name__ == "__main__":
    main()