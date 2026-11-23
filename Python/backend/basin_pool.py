# ==========================================================================
# SAGEhydrology | Python source
# File: backend/basin_pool.py
#
# Headless SAGE numerical backend.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Headless numerical routines extracted from basin_pool."""
import os

def cpu_inventory():
    """Cpu inventory."""
    logical = os.process_cpu_count() if hasattr(os, 'process_cpu_count') else os.cpu_count()
    physical = None
    if os.name == 'nt':
        try:
            import ctypes as C
            n = C.c_ulong()
            fn = C.windll.kernel32.GetLogicalProcessorInformationEx
            fn(0, None, C.byref(n))
            buffer = C.create_string_buffer(n.value)
            if fn(0, buffer, C.byref(n)):
                import struct
                offset = 0
                physical = 0
                while offset < n.value:
                    relation, size = struct.unpack_from('II', buffer.raw, offset)
                    if size < 8:
                        break
                    if relation == 0:
                        physical += 1
                    offset += size
        except Exception:
            physical = None
    return dict(logical=max(1, logical or 1), physical=physical)
