# SAGEhydrology | Python source
# SITEhydrology | Python source
# File: python/installer_success_cleanup.py
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
"""Remove installer-owned working files only after installation succeeds."""

from pathlib import Path
import json
import shutil


def cleanup_success(directory, archives=(), cache=None):
    directory = Path(directory).resolve()
    markers = list(directory.glob("python_installation*.json"))
    if not any(
        json.loads(p.read_text()).get("status") == "installed"
        for p in markers
    ):
        raise ValueError(
            "Cleanup requires a successful installation manifest"
        )
    count = 0
    payload = directory / "source_archives"
    if payload.exists():
        if payload.is_symlink() or payload.resolve().parent != directory:
            raise ValueError("Unsafe extraction cleanup target")
        shutil.rmtree(payload)
    # Transfer sidecars are never runtime dataset inputs.
    for p in directory.rglob("*"):
        if p.is_file() and p.name.endswith(
            (".lock", ".download.json", ".resume.json")
        ):
            p.unlink()
            count += 1
    if cache is None:
        return count
    cache = Path(cache).resolve()
    if cache.parent.name != "_downloads":
        raise ValueError("Unexpected installer download cache")
    for source in map(Path, archives):
        if source.is_symlink():
            continue
        source = source.resolve()
        if not source.is_relative_to(cache):
            # User-selected archives outside our cache are
            # retained.
            continue
        files = list(source.rglob("*")) if source.is_dir() else [source]
        for file in files:
            if not file.is_file() or file.is_symlink():
                continue
            if file.name.endswith((".part", ".resume.json")):
                continue  # Unfinished work remains resumable.
            file.unlink()
            count += 1
            for suffix in (".lock", ".download.json"):
                companion = file.with_name(file.name + suffix)
                if companion.is_file():
                    companion.unlink()
    if cache.exists():
        for folder in sorted(
            (p for p in cache.rglob("*") if p.is_dir()),
            key=lambda p: len(p.parts),
            reverse=True,
        ):
            if not folder.is_symlink() and not any(folder.iterdir()):
                folder.rmdir()
        if not any(cache.iterdir()):
            cache.rmdir()
        if cache.parent.exists() and not any(cache.parent.iterdir()):
            cache.parent.rmdir()
    return count
