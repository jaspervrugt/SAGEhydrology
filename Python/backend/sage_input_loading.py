# ==========================================================================
# SAGEhydrology | Python source
# File: backend/sage_input_loading.py
#
# Headless SAGE numerical backend.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Headless numerical routines extracted from sage_input_loading."""
from collections import deque
from concurrent.futures import ThreadPoolExecutor
from engine import Cancelled

def loading_workers(config):
    """Keep daily or single-worker runs serial; bound hourly I/O memory."""
    if config.get('resolution', 1) == 1:
        return 1
    if not config.get('parallel', True):
        return 1
    return max(1, min(4, int(config.get('workers', 4))))

def prepared_inputs(config, stop, prepare):
    """Yield basin, arrays and error in configured order, reading once.

    Start serially to initialize shared metadata and code identities. Stay
    serial while entries are cached. After a source read, overlap at most
    four remaining reads. NetCDF calls retain their existing native lock;
    source manifests are captured inside each loader's thread context.
    """
    gauges = iter(config['gauges'])
    workers = loading_workers(config)

    def load(gauge):
        """Load."""
        if stop.is_set():
            raise Cancelled()
        try:
            return (gauge, prepare(config, gauge), None)
        except Cancelled:
            raise
        except Exception as error:
            return (gauge, None, error)
    for gauge in gauges:
        if stop.is_set():
            raise Cancelled()
        item = load(gauge)
        yield item
        if workers > 1 and item[1] is not None:
            if not item[1].get('_prepared_cache_hit', False):
                break
    else:
        return
    executor = ThreadPoolExecutor(max_workers=workers, thread_name_prefix='SAGE-input')
    pending = deque()
    try:
        for _ in range(workers):
            gauge = next(gauges, None)
            if gauge is None:
                break
            pending.append(executor.submit(load, gauge))
        while pending:
            if stop.is_set():
                raise Cancelled()
            yield pending.popleft().result()
            if stop.is_set():
                raise Cancelled()
            gauge = next(gauges, None)
            if gauge is not None:
                pending.append(executor.submit(load, gauge))
    finally:
        for future in pending:
            future.cancel()
        executor.shutdown(wait=True, cancel_futures=True)
