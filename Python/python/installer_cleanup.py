# SAGEhydrology | Python source
# SITEhydrology | Python source
# File: python/installer_cleanup.py
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
(
    'Success-only cleanup of installer-owned working files and '
    'verified caches.'
)

from pathlib import Path
import json, shutil


def cleanup_extraction(directory, relative="source_archives"):
    directory = Path(directory).resolve()
    markers = list(directory.glob("python_installation*.json"))
    if not any(
        json.loads(p.read_text()).get("status") == "installed"
        for p in markers
    ):
        return False
    owned = directory / "source_archives"
    target = (directory / relative).resolve()
    if not target.exists():
        return True
    if not target.is_relative_to(owned) or owned.resolve() != owned:
        raise ValueError(
            "Refusing cleanup outside installer extraction directory"
        )
    shutil.rmtree(target)
    if owned.exists() and not any(owned.iterdir()):
        owned.rmdir()
    return True


def cleanup_downloads(archives, cache):
    (
        "Remove only completed verified downloads under this worker's "
        'cache.'
    )
    from dataset_install import _verified_archive

    cache = Path(cache).resolve()
    if cache.parent.name != "_downloads":
        raise ValueError("Unexpected download cache")
    candidates = set()
    for source in map(Path, archives):
        source = source.resolve()
        if not source.is_relative_to(cache):
            continue
        if source.is_file():
            candidates.add(source)
        elif source.is_dir():
            for receipt in source.rglob("*.download.json"):
                candidates.add(
                    receipt.with_name(
                        receipt.name[: -len(".download.json")]
                    )
                )
    removed = 0
    for file in candidates:
        if not file.is_file() or _verified_archive(file) is None:
            continue
        try:
            file.unlink()
            removed += 1
            for suffix in (".download.json", ".lock"):
                companion = file.with_name(file.name + suffix)
                if companion.exists():
                    companion.unlink()
        except OSError:
            # A file still in use is retained, never forced
            # closed.
            continue
    return removed
