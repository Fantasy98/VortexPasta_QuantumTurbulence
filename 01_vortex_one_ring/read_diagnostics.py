"""Read and summarize the CSV written by vortex_ring.jl.

Run: python3 read_diagnostics.py
Or:  python3 read_diagnostics.py path/to/diagnostics.csv

You can also import read_diagnostics() to get a dictionary of numerical lists.
This script uses only Python's standard library and does not modify the CSV.
"""

import argparse
import csv
import math
from pathlib import Path


DEFAULT_CSV = Path(__file__).resolve().parent / "output" / "diagnostics.csv"
COLUMNS = ("time", "energy", "impulse_norm", "vortex_length")


def read_diagnostics(path=DEFAULT_CSV):
    """Return data such as data['time'] and data['energy'] as lists of floats."""
    data = {name: [] for name in COLUMNS}

    with Path(path).open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        missing = set(COLUMNS) - set(reader.fieldnames or [])
        if missing:
            raise ValueError("Missing CSV columns: " + ", ".join(sorted(missing)))

        for row in reader:
            for name in COLUMNS:
                try:
                    value = float(row[name])
                except (TypeError, ValueError):
                    raise ValueError(
                        f"Invalid number in column {name!r}, CSV line {reader.line_num}"
                    ) from None
                if not math.isfinite(value):
                    raise ValueError(
                        f"Non-finite value in column {name!r}, CSV line {reader.line_num}"
                    )
                data[name].append(value)

    if not data["time"]:
        raise ValueError("The CSV contains no data rows")
    return data


def print_summary(data):
    times = data["time"]
    print(f"Samples: {len(times)}")
    print(f"Time range: {times[0]:.12g} to {times[-1]:.12g}")
    print()
    print(f"{'Quantity':<16} {'Initial':>18} {'Final':>18} {'Final rel. change':>19} {'Max |rel. change|':>19}")
    for name in COLUMNS[1:]:
        values = data[name]
        initial, final = values[0], values[-1]
        if initial == 0:
            final_change = maximum_change = "n/a (initial = 0)"
        else:
            changes = [(value - initial) / initial for value in values]
            final_change = f"{changes[-1]:.6e}"
            maximum_change = f"{max(abs(change) for change in changes):.6e}"
        print(f"{name:<16} {initial:18.12g} {final:18.12g} {final_change:>19} {maximum_change:>19}")
    print("Relative change = (value - initial) / initial; these are fractions, not percentages.")

    print("\nFirst five rows (or all rows if fewer):")
    print(" ".join(f"{name:>18}" for name in COLUMNS))
    for i in range(min(5, len(times))):
        print(" ".join(f"{data[name][i]:18.12g}" for name in COLUMNS))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path, nargs="?", default=DEFAULT_CSV)
    args = parser.parse_args()
    try:
        data = read_diagnostics(args.csv)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    print(f"File: {args.csv.resolve()}")
    print_summary(data)


if __name__ == "__main__":
    main()
