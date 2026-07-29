#!/usr/bin/env python3
"""Convert raw RAMSES AMR outputs into reusable regular HDF5 cube caches."""

from __future__ import annotations

import argparse
import hashlib
import os
import re
import sys
import tempfile
from pathlib import Path

try:
    import h5py
    import numpy as np
    import yt
except ImportError as error:
    missing = getattr(error, "name", "yt/h5py/numpy")
    raise SystemExit(
        f"Missing Python package {missing!r}. Install the converter dependencies "
        "with: python3 -m pip install --user -r requirements-ramses.txt"
    ) from error


DENSITY_STORAGE_UNIT = 1.0e-12
OUTPUT_PATTERN = re.compile(r"^output_(\d+)$")


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Use yt to resample RAMSES outputs on a regular grid and cache the "
            "fields required by the Julia diagnostics."
        )
    )
    parser.add_argument("--simulation", required=True, type=Path)
    parser.add_argument("--output-directory", required=True, type=Path)
    parser.add_argument("--resolution", required=True, type=int)
    parser.add_argument(
        "--outputs",
        default="all",
        help="Comma-separated output_XXXXX directory names, or 'all'.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Rebuild cache files even when their source fingerprint is unchanged.",
    )
    return parser.parse_args()


def discover_outputs(simulation: Path, requested: str) -> list[Path]:
    outputs = []
    requested_names = None if requested.lower() == "all" else {
        value.strip() for value in requested.split(",") if value.strip()
    }
    for path in sorted(simulation.iterdir()):
        matched = OUTPUT_PATTERN.match(path.name)
        if not path.is_dir() or matched is None:
            continue
        if requested_names is not None and path.name not in requested_names:
            continue
        info = path / f"info_{matched.group(1)}.txt"
        if info.is_file():
            outputs.append(path)
    if requested_names is not None:
        found = {path.name for path in outputs}
        missing = sorted(requested_names - found)
        if missing:
            raise SystemExit(
                "Requested RAMSES outputs were not found or lack info files: "
                + ", ".join(missing)
            )
    if not outputs:
        raise SystemExit(f"No readable output_XXXXX directories found in {simulation}")
    return outputs


def available_fields(dataset) -> set[tuple[str, str]]:
    return set(dataset.field_list) | set(dataset.derived_field_list)


def select_field(dataset, aliases: tuple[str, ...]) -> tuple[str, str]:
    available = available_fields(dataset)
    for field_type in ("gas", "ramses"):
        for alias in aliases:
            candidate = (field_type, alias)
            if candidate in available:
                return candidate
    normalized_aliases = {
        re.sub(r"[^a-z0-9]", "", alias.lower()) for alias in aliases
    }
    matches = [
        field
        for field in available
        if re.sub(r"[^a-z0-9]", "", field[1].lower()) in normalized_aliases
    ]
    if len(matches) == 1:
        return matches[0]
    raise RuntimeError(
        f"Could not uniquely resolve aliases {aliases}. "
        f"Available fluid fields: {sorted(available)}"
    )


def grid_for_resolution(dataset, resolution: int):
    dimensions = np.full(3, resolution, dtype=np.int64)
    base_dimensions = np.asarray(dataset.domain_dimensions, dtype=np.int64)
    ratios = dimensions / base_dimensions
    if (
        np.all(dimensions % base_dimensions == 0)
        and np.allclose(ratios, ratios[0])
        and ratios[0] >= 1
        and float(np.log2(ratios[0])).is_integer()
    ):
        level = int(round(np.log2(ratios[0])))
        return dataset.covering_grid(
            level=level,
            left_edge=dataset.domain_left_edge,
            dims=dimensions,
        )
    return dataset.arbitrary_grid(
        dataset.domain_left_edge,
        dataset.domain_right_edge,
        dimensions,
    )


def field_array(grid, dataset, aliases: tuple[str, ...], unit: str) -> np.ndarray:
    field = select_field(dataset, aliases)
    return np.asarray(grid[field].to_value(unit), dtype=np.float32)


def magnetic_array(grid, dataset, axis: str) -> np.ndarray:
    try:
        return field_array(
            grid,
            dataset,
            (f"magnetic_field_{axis}", f"B_{axis}"),
            "gauss",
        )
    except RuntimeError:
        left = field_array(
            grid,
            dataset,
            (f"magnetic_field_{axis}_left", f"B_{axis}_left"),
            "gauss",
        )
        right = field_array(
            grid,
            dataset,
            (f"magnetic_field_{axis}_right", f"B_{axis}_right"),
            "gauss",
        )
        return np.float32(0.5) * (left + right)


def source_fingerprint(output: Path) -> str:
    digest = hashlib.blake2b(digest_size=20)
    relevant_prefixes = (
        "amr_",
        "hydro_",
        "info_",
        "hydro_file_descriptor",
    )
    for path in sorted(output.iterdir()):
        if not path.is_file() or not path.name.startswith(relevant_prefixes):
            continue
        stat = path.stat()
        digest.update(path.name.encode())
        digest.update(str(stat.st_size).encode())
        digest.update(str(stat.st_mtime_ns).encode())
    return digest.hexdigest()


def cache_is_current(
    target: Path, output: Path, info_path: Path, resolution: int, force: bool
) -> bool:
    if force or not target.is_file():
        return False
    fingerprint = source_fingerprint(output)
    try:
        with h5py.File(target, "r") as handle:
            return (
                int(handle.attrs.get("cache_format_version", 0)) == 1
                and int(handle.attrs.get("resolution", -1)) == resolution
                and str(handle.attrs.get("source_fingerprint", "")) == fingerprint
                and str(handle.attrs.get("source_info", "")) == str(info_path.resolve())
            )
    except OSError:
        return False


def write_dataset(handle, name: str, values: np.ndarray) -> None:
    handle.create_dataset(
        name,
        data=values,
        compression="lzf",
        shuffle=True,
        chunks=True,
    )


def convert_output(output: Path, target: Path, resolution: int, force: bool) -> str:
    number = OUTPUT_PATTERN.match(output.name).group(1)
    info_path = output / f"info_{number}.txt"
    if cache_is_current(target, output, info_path, resolution, force):
        return "cached"

    dataset = yt.load(str(info_path))
    grid = grid_for_resolution(dataset, resolution)
    time_myr = float(dataset.current_time.to_value("Myr"))
    fingerprint = source_fingerprint(output)

    target.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{target.name}.", suffix=".tmp", dir=target.parent
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        with h5py.File(temporary, "w") as handle:
            field_builders = (
                (
                    "rho",
                    lambda: field_array(
                        grid, dataset, ("density", "Density"), "g/cm**3"
                    )
                    / np.float32(DENSITY_STORAGE_UNIT),
                ),
                (
                    "P",
                    lambda: field_array(
                        grid,
                        dataset,
                        ("pressure", "Pressure", "thermal_pressure"),
                        "erg/cm**3",
                    ),
                ),
                (
                    "vx",
                    lambda: field_array(
                        grid, dataset, ("velocity_x", "x-velocity"), "km/s"
                    ),
                ),
                (
                    "vy",
                    lambda: field_array(
                        grid, dataset, ("velocity_y", "y-velocity"), "km/s"
                    ),
                ),
                (
                    "vz",
                    lambda: field_array(
                        grid, dataset, ("velocity_z", "z-velocity"), "km/s"
                    ),
                ),
                ("bx", lambda: magnetic_array(grid, dataset, "x")),
                ("by", lambda: magnetic_array(grid, dataset, "y")),
                ("bz", lambda: magnetic_array(grid, dataset, "z")),
            )
            for name, build in field_builders:
                values = build()
                write_dataset(handle, name, values)
                del values
                clear_data = getattr(grid, "clear_data", None)
                if clear_data is not None:
                    clear_data()
            handle.create_dataset("L", data=np.full(3, 100.0, dtype=np.float64))
            handle.create_dataset("t", data=np.float64(time_myr))
            handle.attrs["cache_format_version"] = 1
            handle.attrs["source_format"] = "RAMSES/yt"
            handle.attrs["source_info"] = str(info_path.resolve())
            handle.attrs["source_fingerprint"] = fingerprint
            handle.attrs["resolution"] = resolution
            handle.attrs["yt_version"] = yt.__version__
            handle.attrs["rho_stored_unit_g_cm3"] = DENSITY_STORAGE_UNIT
            handle.flush()
        os.replace(temporary, target)
    finally:
        temporary.unlink(missing_ok=True)
    return "converted"


def main() -> int:
    arguments = parse_arguments()
    if arguments.resolution <= 0:
        raise SystemExit("--resolution must be a positive integer")
    simulation = arguments.simulation.expanduser().resolve()
    output_directory = arguments.output_directory.expanduser().resolve()
    if not simulation.is_dir():
        raise SystemExit(f"Simulation directory not found: {simulation}")
    outputs = discover_outputs(simulation, arguments.outputs)
    output_directory.mkdir(parents=True, exist_ok=True)

    print(f"RAMSES source     : {simulation}")
    print(f"Cache directory  : {output_directory}")
    print(f"Uniform grid     : {arguments.resolution}^3")
    print(f"Selected outputs : {len(outputs)}")
    for index, output in enumerate(outputs, start=1):
        number = OUTPUT_PATTERN.match(output.name).group(1)
        target = output_directory / f"info_{number}.h5"
        status = convert_output(
            output, target, arguments.resolution, arguments.force
        )
        print(f"[{index:>4}/{len(outputs)}] {output.name} -> {target.name} ({status})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
