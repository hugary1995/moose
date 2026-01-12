#!/usr/bin/env python3

import argparse
import json
from pathlib import Path
from tqdm import tqdm


def extract_header_times(trace_files: list[Path], output_csv: Path, MOOSE_ROOT: Path):
    """Extract and accumulate compile times for individual header files from compile time trace JSON files Dump the results to a CSV file. Header files outside the MOOSE root directory are ignored. :param trace_files: List of compile time trace JSON files :param output_csv: Path to the output CSV file :param MOOSE_ROOT: Path to the MOOSE root directory"""

    header_times = {}

    def fixed_width_desc(s: str, width: int = 30) -> str:
        if len(s) <= width:
            return s.ljust(width)  # pad so bar doesn't jitter
        return s[: width - 3] + "..."

    for trace_file in (pbar := tqdm(trace_files)):
        pbar.set_description(f"Processing {fixed_width_desc(trace_file.name)}")
        with open(trace_file, "r") as f:
            trace_data = json.load(f)
            for event in trace_data.get("traceEvents", []):
                dur = event.get("dur")
                if dur is None:
                    continue

                args = event.get("args", {})
                detail = args.get("detail", "")
                if not detail:
                    continue

                # The JSON trace file does not have a fixed schema across compiler versions
                # The rule of thumb we are using here is to look at each trace event and find those
                #   - with a "dur" field (duration)
                #   - with at least one "args.detail" field that contains a header file path
                header_file = Path(detail)
                if header_file.is_file() and header_file.suffix in {".h"}:
                    header_file_str = str(header_file)
                    if not header_file_str.startswith(str(MOOSE_ROOT)):
                        continue
                    if header_file_str not in header_times:
                        header_times[header_file_str] = 0
                    header_times[header_file_str] += int(dur)

    # Dump the results to a CSV file
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(output_csv, "w") as f:
        f.write("header_file,compile_time\n")
        for header_file, compile_time in header_times.items():
            f.write(f"{header_file},{compile_time}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.description = "Analyze compile time trace JSON files"
    parser.epilog = (
        "To use this script, first generate compile time trace JSON files using the appropriate compiler flags. Note that the clang toolchain natively supports generating such traces. GCC users may need to use additional tools to capture compile time traces. Example command to generate compile time trace:"
        "\n  MOOSE_UNITY=0 MOOSE_HEADER_SYMLINKS=0 make -j N CXXFLAGS='-ftime-trace'\n"
        "Note the use of MOOSE_UNITY=0 to disable unity builds and MOOSE_HEADER_SYMLINKS=0 to disable header symlinks, which can obscure individual file compile times. Trace files will be located in .libs directories in each source folder. This script then parses the specified JSON trace files and provides useful insights into compile times."
    )
    parser.add_argument(
        "--directories",
        "-d",
        nargs="+",
        type=Path,
        help="Directories inside which to search for compile time trace JSON files. Use either absolute or relative paths. If a relative path is provided, it is considered relative to the MOOSE root directory.",
        default=[Path("framework"), Path("modules")],
    )
    parser.add_argument(
        "--exclude-dirs",
        "-e",
        nargs="+",
        type=Path,
        help="Directories to exclude from the search for compile time trace JSON files. Use either absolute or relative paths. If a relative path is provided, it is considered relative to the MOOSE root directory.",
        default=[Path("framework/contrib")],
    )
    parser.add_argument(
        "--list-files",
        "-l",
        action="store_true",
        help="List all found compile time trace JSON files without further analysis.",
    )
    parser.add_argument(
        "--extract-header-times",
        nargs=1,
        type=Path,
        help="Extract compile times accumulated for individual header files and dump to the specified CSV file. If a relative path is provided, it is considered relative to the current working directory.",
        default=None,
    )
    parser.add_argument(
        "--n-threads",
        type=int,
        help="Number of threads to use for various tasks such as processing trace files. Default is min(32, os.cpu_count() + 4)",
        default=None,
    )

    args = parser.parse_args()

    # MOOSE root directory
    MOOSE_ROOT = Path(__file__).resolve().parent.parent

    # Find all compile time trace JSON files
    trace_files = []
    exclude_dirs = [
        MOOSE_ROOT / d if not d.is_absolute() else d for d in args.exclude_dirs
    ]
    for dir in args.directories:
        if not dir.is_absolute():
            dir = MOOSE_ROOT / dir
        if not dir.exists():
            raise FileNotFoundError(f"Directory {dir} does not exist.")
        files = list(dir.rglob(".libs/*.json"))
        for file in files:
            if not any(excl in file.parents for excl in exclude_dirs):
                trace_files.append(file)

    # If --list-files is specified, print the found files and exit
    if args.list_files:
        for file in trace_files:
            print(file)
        exit(0)

    # Extract a database of compile times
    if args.extract_header_times is not None:
        ht_csv = Path(args.extract_header_times[0])
        if not ht_csv.is_absolute():
            ht_csv = Path.cwd() / ht_csv
        extract_header_times(trace_files, ht_csv, MOOSE_ROOT)
