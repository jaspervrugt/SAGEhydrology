# ==========================================================================
# SAGEhydrology | Python source
# File: gui_python/sage_log.py
#
# Stable central log formatting shared by the worker and GUI.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Stable central log formatting shared by the worker and GUI."""


def iteration(e):
    """Format one SAGE iteration log entry."""
    return (
        f" {e['iteration']:9d}  {e['objective']:14.6g} "
        f" {e['gradient_norm']:14.6g}  {e['learning_rate']:14.6g} "
        f" {e['seconds']:11.2f}"
    )


ITERATION_HEADER = (
    " Iteration       Objective      Gradient norm    Learning rate    "
    " Time (s)"
)


def progress_message(message, done, total):
    if done < total:
        message += f" [{done}/{total},{100*done/total:5.1f}% done]"
    return message


def postprocessing_line(stage, done, total):
    return progress_message(f"({stage:02d}) postprocess basins ...", done, total)


def event(e):
    """Format a structured run event for the log."""
    kind = e["kind"]
    if kind == "log":
        return e["message"]
    if kind == "error":
        return "ERROR: " + e["message"]
    if kind == "sage_iteration":
        return ""  # Dedicated Iter view and demo_SAGE_iterations.txt.
    if kind == 'training_complete':
        from sage_model_selection import progress_text, selection_text
        return (progress_text(e['done'],e['total'],e.get('stage',14),
                              True,e.get('status','completed'))+'\n'
                +selection_text(e.get('selection'),e['done']))
    if kind == "data_progress":
        done, total = e["done"], e["total"]
        stage = e.get("stage", 7)
        label = "meteorological" if stage == 7 else "discharge"
        if done == total:
            return f"({stage:02d}) read {label} data ..."
    if kind == "postprocessing_progress" and e["done"] == e["total"]:
        return postprocessing_line(e["stage"], e["done"], e["total"])
    if kind == "pool":
        return (
            f'      Basin workers: {e["workers"]}; maximum outstanding'
            f' jobs: {e["max_pending"]}; available logical CPUs:'
            f' {e["logical"]}.'
        )
    if kind == "finished":
        return (
            "... SAGE run "
            + e["status"]
            + "\n      Output directory: "
            + e["directory"]
        )
    return ""


BANNER_WIDTH = 60
BANNER_INDENT = 11
LOG_WIDTH = BANNER_INDENT + BANNER_WIDTH + 2


def header(c, *, compact=False, width=BANNER_WIDTH):
    """Port of MATLAB src/write_header.m's fixed-width run banner."""
    import textwrap
    from site_schema import MODEL_NAMES
    from datetime import datetime, timedelta
    from site_schema import region_catalogue, LOSSES

    meta = region_catalogue(c["sage_root"])[c["region"]]
    count = len(c["gauges"])
    training = (
        max(1, min(count, int(count * c["basin_train_fraction"] + 1e-10)))
        if count
        else 0
    )
    resolution = c["resolution"]

    def date(key, end=False):
        """Format a run date for display."""
        raw = c[key]
        value = datetime.fromisoformat(raw)
        if resolution > 1 and end and len(raw) == 10:
            value += timedelta(days=1, minutes=-1440 / resolution)
        return value

    def span(first, last):
        """Format a run date range for display."""
        fmt = "%d %b %Y" if resolution == 1 else "%d %b %Y %H:%M"
        return first.strftime(fmt) + " → " + last.strftime(fmt)

    method = c["split"]
    manual = method == "manual"
    has_eval = c.get("evaluation_enabled", True)
    if manual:
        first, last = date("train_start"), date("train_end", True)
        training_period = span(first, last)
        evaluation_period = "Not configured"
        if has_eval:
            ef, el = date("eval_start"), date("eval_end", True)
            evaluation_period = span(ef, el)
            first, last = min(first, ef), max(last, el)
    else:
        first, last = date("common_start"), date("common_end", True)
    methods = dict(
        manual="Manual dates",
        deterministic_block="Deterministic block",
        random_block="Random block",
        random="Random",
        deterministic_kfold="Deterministic k-fold",
        random_kfold="Random k-fold",
        rainfall_block="Rainfall block",
    )
    optimizer = (
        ("AdamW" if c.get("weight_decay", 0) > 0 else "Adam")
        if c["optimizer"] == 2
        else "Gradient descent"
    )
    rows = [
        ("Model", MODEL_NAMES.get(c["model"], str(c["model"])).lower()),
        ("Language", "Python / NumPy + native C++"),
        ("Region", meta["name"].split(" (")[0]),
        ("Dataset", meta["dataset"]),
        (
            "Resolution",
            {1: "Daily", 24: "Hourly", 96: "15-minute"}.get(
                resolution, f"{resolution} samples/day"
            ),
        ),
        (
            "Basins",
            f"{count} (training: {training}; evaluation: {count-training})",
        ),
        ("Period", span(first, last)),
        ("Time split", methods.get(method, method)),
    ]
    if manual:
        rows.extend(
            [
                ("Training", training_period),
                ("Evaluation", evaluation_period),
            ]
        )
    project = c.get("publication_project", {})
    if project:
        details = [
            ("Project", project.get("name", "")),
            ("Title", project.get("paper_title", "")),
            (
                "Summary",
                project.get("abstract", project.get("summary", "")),
            ),
        ]
        rows[1:1] = [(key, value) for key, value in details if value]
    names = {
        "Q": "Discharge",
        "SWE": "Snow-water equivalent",
        "SM": "Soil-water storage",
    }
    rows.extend(
        [
            ("Loss function", LOSSES[c["loss"]]),
            (
                "Data types",
                ", ".join(
                    names.get(name, name)
                    for name in c.get("observations", ["Q"])
                ),
            ),
            ("Optimizer", optimizer),
        ]
    )
    width = max(36, int(width))
    indent = "" if compact else " " * BANNER_INDENT
    label_width = max(len(label) for label, _ in rows)

    def line(text="", center=False):
        # Consolas falls back to a nearly two-cell glyph for the envelope.
        """Format a labeled line in the run banner."""
        display = len(text) + text.count(chr(0x2709))
        left = (width - display) // 2 if center else 2
        return (
            indent
            + "|"
            + " " * left
            + text
            + " " * max(0, width - left - display)
            + "|"
        )

    result = [
        indent + "┌" + "─" * width + "┐",
        line("S A G E", True),
        line(),
        line("Sensitivity-Aware Learning for Rapid Training of", True),
        line(
            "Process-Based Hydrologic Models across Large Basin Samples",
            True,
        ),
        line(),
    ]
    if compact:
        result = [
            indent + "┌" + "─" * width + "┐",
            line("SAGEhydrology run summary", True),
            line(),
        ]
    for label, value in rows:
        wrapped = textwrap.wrap(
            str(value), width=width - 4 - label_width - 3
        ) or [""]
        for j, part in enumerate(wrapped):
            result.append(
                line(
                    f"{label:<{label_width}} : {part}"
                    if j == 0
                    else " " * (label_width + 3) + part
                )
            )
    result.extend(
        [
            line(),
            line("© 2026 Jasper A. Vrugt"),
            line("  All rights reserved"),
            line(),
            line("jasper@uci.edu" if compact else "✉ jasper@uci.edu"),
            indent + "└" + "─" * width + "┘",
        ]
    )
    return "\n".join(result)


def bounded_lines(message):
    """Keep transcript rows within the standard run banner."""
    from textwrap import wrap

    for line in str(message).split("\n"):
        if len(line) <= LOG_WIDTH:
            yield line
        else:
            indent = " " * min(len(line) - len(line.lstrip()), 11)
            yield from wrap(line.lstrip(), width=LOG_WIDTH,
                            initial_indent=indent, subsequent_indent=indent or "      ",
                            break_on_hyphens=False)


def timestamped(message, now=None):
    """
    Timestamp every displayed line, preserving already stamped
    transcripts.
    """
    import re
    from datetime import datetime

    stamp = (now or datetime.now()).strftime("[%Y-%m-%d %H:%M:%S] ")
    rows = []
    for line in str(message).split("\n"):
        existing = re.match(r"^\[\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}\] ?", line)
        prefix = existing.group().rstrip() + " " if existing else stamp
        body = line[existing.end():] if existing else line
        rows.extend(prefix + part for part in bounded_lines(body))
    return "\n".join(rows)


def automatic_actions(note):
    """Align every action with its heading after the log timestamp."""
    lines = ["Automatic data actions", *str(note).splitlines()]
    return "\n".join("      " + line.lstrip() for line in lines)
