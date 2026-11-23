# ==========================================================================
# SAGEhydrology | Python source
# File: run_sage.py
#
# Portable headless SAGE entry point.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Run a saved SAGE session without the graphical application."""

from pathlib import Path
import argparse
import json
import multiprocessing as mp
import sys

ROOT = Path(__file__).resolve().parent
for folder in (ROOT / "python", ROOT / "backend"):
    sys.path.insert(0, str(folder))


class Events:
    """Consume worker notifications without retaining diagnostic arrays."""

    def __init__(self, verbose=False):
        """Select whether to print completed iteration summaries."""
        self.verbose = verbose
        self.error = None

    def put(self, event):
        """Discard transient arrays and optionally report progress."""
        if event.get("kind") == "error" and event.get("fatal"):
            self.error = event["message"]
        if self.verbose and event.get("kind") == "sage_iteration":
            print("iteration", event["iteration"], "objective",
                  format(event["objective"], ".6g"), flush=True)


def build_native():
    """Compile public numerical libraries with an installed C++ compiler."""
    from build_crr import build
    from build_native import build as hymod
    from build_hbv_routing import build as routing
    from build_observables import build as observables

    for action in (build, hymod, routing, observables):
        action(ROOT.parent / "Native")


def execute(config, output, progress=False):
    """Run a public model and return its persisted result report.

    Parameters
    ----------
    config : dict
        SAGE session settings, including dataset paths and basin IDs.
    output : path-like
        New directory for this run. Existing runs are never overwritten.

    Returns
    -------
    dict
        Completed run report, including training and evaluation results.
    """
    from sage_engine import run

    if int(config.get("model", 1)) not in (*range(1, 8), 99):
        raise ValueError(
            "The Python session runner supports public models 1-7 and user_model (99)."
        )
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(
            "Choose a new output directory: " + str(output)
        )
    events = Events(progress)
    run(config, str(output), events, mp.get_context("spawn").Event())
    if events.error is not None:
        raise RuntimeError(events.error)
    report = json.loads((output / "run.json").read_text(encoding="utf-8"))
    if report["status"] != "complete":
        raise RuntimeError(json.dumps(report.get("errors", [])))
    return report


def main():
    """Build native cores or execute a saved JSON session headlessly."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-native", action="store_true")
    parser.add_argument("--session", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.build_native:
        build_native()
    if args.session:
        if args.output is None:
            parser.error("--session requires --output")
        from sage_schema import load_session

        report = execute(load_session(args.session), args.output, progress=True)
        print("Completed", len(report["history"]), "iterations.")
    elif not args.build_native:
        parser.error("Choose --build-native or --session.")


if __name__ == "__main__":
    mp.freeze_support()
    main()
