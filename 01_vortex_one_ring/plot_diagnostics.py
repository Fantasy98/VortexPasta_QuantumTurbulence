"""Save vortex_ring.jl diagnostics as a PNG without a graphical display.

Run after the simulation: python3 plot_diagnostics.py
Or choose a CSV file: python3 plot_diagnostics.py path/to/diagnostics.csv
Uses the existing NumPy and Matplotlib installations; installs nothing.
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "csv", type=Path, nargs="?",
        default=Path(__file__).resolve().parent / "output" / "diagnostics.csv",
    )
    args = parser.parse_args()
    if not args.csv.is_file():
        parser.error(f"Missing {args.csv}; run vortex_ring.jl --no-plot first")
    data = np.atleast_1d(np.genfromtxt(args.csv, delimiter=",", names=True))
    if data.size == 0:
        parser.error("The diagnostics CSV contains no samples")

    fig, ax = plt.subplots(figsize=(8, 5), layout="constrained")
    for name, label, color in (
        # ("energy", "Energy", "#087f8c"),
        ("impulse_norm", "Impulse", "#d17622"),
        ("vortex_length", "Vortex length", "#7651a8"),
    ):
        values = data[name]
        if not np.all(np.isfinite(values)) or values[0] == 0:
            parser.error(f"Cannot normalize {name}: non-finite values or zero initial value")
        ax.plot(data["time"], values / values[0], label=label, color=color, linewidth=1.5)
    ax.set(xlabel="Time", ylabel="Value / initial value", title="Vortex ring diagnostics")
    # ax.ticklabel_format(axis="y", style="plain", useOffset=False,)
    ax.grid(alpha=0.25)
    ax.legend()
    output = args.csv.with_name("energy.png")
    fig.savefig(output, dpi=200)
    plt.close(fig)
    print(f"Saved figure: {output.resolve()}")


if __name__ == "__main__":
    main()
