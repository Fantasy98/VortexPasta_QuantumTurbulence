"""Render saved VortexPasta snapshots to PNGs and a GIF without a display.

Run: python3 render_vortex_ring.py
Requires existing Julia/VortexPasta, NumPy, Matplotlib and Pillow; installs nothing.
The simulation is not rerun. Times and snapshot order come from the .series file.
"""

import argparse
import json
import os
from pathlib import Path
import subprocess

import matplotlib

matplotlib.use("Agg")  # Render to files without X11, OpenGL or a DISPLAY.
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


PROJECT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--series", type=Path,
        default=PROJECT / "output" / "vortex_ring.vtkhdf.series",
    )
    parser.add_argument("--output", type=Path, help="Default: beside the .series file")
    parser.add_argument("--fps", type=int, default=15, help="GIF playback frames/second")
    parser.add_argument("--julia", default="julia", help="Existing Julia executable")
    args = parser.parse_args()
    if not 1 <= args.fps <= 50:
        parser.error("--fps must be between 1 and 50")

    series_path = args.series.resolve()
    with series_path.open() as stream:
        entries = sorted(json.load(stream)["files"], key=lambda entry: entry["time"])
    if not entries:
        parser.error("The .series file contains no snapshots")
    filenames = [series_path.parent / entry["name"] for entry in entries]
    for filename in filenames:
        if not filename.is_file():
            parser.error(f"Snapshot is missing: {filename}")

    print(f"Reading {len(filenames)} saved snapshots with Julia...", flush=True)
    command = [
        args.julia, "--startup-file=no", f"--project={PROJECT}",
        str(PROJECT / "read_vortex_snapshots.jl"),
        *map(str, filenames),
    ]
    result = subprocess.run(
        command, check=True, text=True, stdout=subprocess.PIPE,
        env={**os.environ, "JULIA_PKG_OFFLINE": "true"},
    )
    frames = [
        [np.asarray(points, dtype=float) for points in json.loads(line)]
        for line in result.stdout.splitlines()
    ]
    if len(frames) != len(entries):
        raise RuntimeError("The number of decoded snapshots does not match the series")

    outdir = (args.output or series_path.parent).resolve()
    frames_dir = outdir / "frames"
    frames_dir.mkdir(parents=True, exist_ok=True)

    # Keep the camera, coordinate limits, and spatial scale identical in all frames.
    all_points = np.concatenate([points for frame in frames for points in frame])
    low, high = all_points.min(axis=0), all_points.max(axis=0)
    centre = (low + high) / 2
    half_width = max(float((high - low).max()) * 0.62, 1e-6)
    fig = plt.figure(figsize=(8, 6), dpi=120, facecolor="white")
    ax = fig.add_subplot(111, projection="3d")
    for setter, value in zip((ax.set_xlim, ax.set_ylim, ax.set_zlim), centre):
        setter(value - half_width, value + half_width)
        # setter(value - half_width, value + half_width)
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=23, azim=-55)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    for i, points in enumerate(frames[0]):
        ax.plot(*points.T, color="0.65", linestyle="--", linewidth=1.2,
                label="Initial position" if i == 0 else None)
    lines = [
        ax.plot([], [], [], color="#087f8c", linewidth=2.8,
                label="Current position" if i == 0 else None)[0]
        for i in range(max(map(len, frames)))
    ]
    ax.legend(loc="upper left", frameon=False, fontsize=9)
    fig.subplots_adjust(left=0.02, right=0.94, bottom=0.09, top=0.89)
    title = fig.suptitle("", fontsize=15)
    fig.text(0.5, 0.035, "Fixed camera and equal coordinate scales",
             ha="center", color="0.4", fontsize=9)

    gif_frames = []
    palette = None
    for index, (entry, frame) in enumerate(zip(entries, frames)):
        for line, points in zip(lines, frame):
            line.set_data_3d(*points.T)
            line.set_visible(True)
        for line in lines[len(frame):]:
            line.set_visible(False)
        title.set_text(f"Vortex ring  |  t = {entry['time']:.4f}")
        png_path = frames_dir / f"frame_{index:04d}.png"
        fig.savefig(png_path)
        with Image.open(png_path) as source:
            rgb = source.convert("RGB")
            if palette is None:
                palette = rgb.quantize(colors=256)
            gif_frames.append(rgb.quantize(palette=palette, dither=Image.Dither.NONE))
        if index % 10 == 0 or index == len(frames) - 1:
            print(f"Rendered {index + 1}/{len(frames)}", flush=True)
    plt.close(fig)

    # GIF times have 10 ms resolution. Playback speed is independent of simulation time.
    duration_ms = max(10, round(100 / args.fps) * 10)
    gif_path = outdir / "vortex_ring_vfm.gif"
    gif_frames[0].save(
        gif_path, save_all=True, append_images=gif_frames[1:],
        duration=duration_ms, loop=0, disposal=2, optimize=False,
    )
    print(f"PNG snapshots: {frames_dir}")
    print(f"GIF: {gif_path}")


if __name__ == "__main__":
    main()
