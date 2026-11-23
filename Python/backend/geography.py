# ==========================================================================
# SAGEhydrology | Python source
# File: backend/geography.py
#
# Headless SAGE numerical backend.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Headless numerical routines extracted from geography."""
import csv
from functools import lru_cache
from pathlib import Path
import struct
import numpy as np
from matplotlib.path import Path as MPath
from matplotlib.patches import PathPatch
LABELS = {1: 'Cold Snow', 2: 'Cold Mountain', 3: 'Mountain', 4: 'Dry Mountain', 6: 'Humid Upland', 7: 'Semi-arid Upland', 8: 'Dry Upland', 9: 'Wet Upland', 10: 'Hot Dry', 11: 'Dry', 12: 'Semi-arid', 13: 'Tropical Wet', 14: 'Wet', 15: 'Tropical Humid', 16: 'Subtropical Humid', 17: 'Temperate Humid', 18: 'Boreal Humid'}

def read_attributes(root):
    """Read attributes."""
    result = {}
    for name in ('camels_topo.txt', 'camels_clim.txt'):
        file = Path(root) / name
        if file.exists():
            with file.open(encoding='utf-8-sig') as handle:
                for row in csv.DictReader(handle, delimiter=';'):
                    result.setdefault(row['gauge_id'], {}).update(row)
    return result

def hydroclimatic_zones(attributes, min_per_zone=5, max_zones=8):
    """Hydroclimatic zones."""
    ids = sorted(attributes)
    labels = []
    features = []

    def num(row, key, default=np.nan):
        """Num."""
        try:
            return float(row.get(key, default))
        except (TypeError, ValueError):
            return default
    for gauge in ids:
        r = attributes[gauge]
        a = num(r, 'aridity')
        e = num(r, 'elev_mean')
        snow = num(r, 'frac_snow')
        lat = abs(num(r, 'gauge_lat'))
        moisture = 'Dry' if a >= 2 else 'Semi-arid' if a >= 1.2 else 'Humid' if a >= 0.65 or not np.isfinite(a) else 'Wet'
        climate = 'Tropical' if lat < 23 else 'Subtropical' if lat < 35 else 'Temperate' if lat < 55 or not np.isfinite(lat) else 'Boreal'
        cold = snow > 0.25 or num(r, 't_mean', np.inf) < 1
        if e >= 1500:
            zone = 4 if moisture == 'Dry' else 2 if cold else 3
        elif cold:
            zone = 1
        elif e >= 500:
            zone = {'Dry': 8, 'Semi-arid': 7, 'Wet': 9}.get(moisture, 6)
        else:
            zone = 10 if moisture == 'Dry' and climate == 'Tropical' else 11 if moisture == 'Dry' else 12 if moisture == 'Semi-arid' else (13 if climate == 'Tropical' else 14) if moisture == 'Wet' else {'Tropical': 15, 'Subtropical': 16, 'Temperate': 17, 'Boreal': 18}[climate]
        labels.append(LABELS[zone])
        features.append([a, e, snow, lat])
    if not ids:
        return {}
    labels = np.array(labels, object)
    X = np.array(features, float)
    valid = np.isfinite(X)
    count = valid.sum(axis=0)
    mean = np.divide(np.where(valid, X, 0).sum(axis=0), count, out=np.zeros(X.shape[1]), where=count > 0)
    residual = np.where(valid, X - mean, 0)
    variance = np.divide((residual ** 2).sum(axis=0), count - 1, out=np.zeros(X.shape[1]), where=count > 1)
    std = np.sqrt(variance)
    std[std == 0] = 1.0
    X = np.where(valid, (X - mean) / std, np.nan)

    def center_of(values):
        """Center of."""
        finite = np.isfinite(values)
        n = finite.sum(axis=0)
        return np.divide(np.where(finite, values, 0).sum(axis=0), n, out=np.full(values.shape[1], np.nan), where=n > 0)

    def unique():
        """Unique."""
        return list(dict.fromkeys(labels))
    while True:
        names = unique()
        counts = [int(np.sum(labels == n)) for n in names]
        if len(names) <= 1 or (min(counts) >= min_per_zone and len(names) <= max(1, len(ids) // min_per_zone)):
            break
        small = int(np.argmin(counts))
        candidates = [i for i in range(len(names)) if i != small]
        large = [i for i in candidates if counts[i] >= min_per_zone]
        candidates = large or candidates
        center = center_of(X[labels == names[small]])
        distances = []
        for i in candidates:
            other = center_of(X[labels == names[i]])
            valid = np.isfinite(center) & np.isfinite(other)
            distances.append(float(np.linalg.norm(center[valid] - other[valid])) if valid.any() else np.inf)
        labels[labels == names[small]] = names[candidates[int(np.argmin(distances))]]
    while len(unique()) > max_zones:
        names = unique()
        counts = [int(np.sum(labels == n)) for n in names]
        small = int(np.argmin(counts))
        candidates = [i for i in range(len(names)) if i != small]
        big = max(candidates, key=lambda i: counts[i])
        labels[labels == names[small]] = names[big]
    return dict(zip(ids, labels.tolist()))

def _dbf_names(path):
    """Dbf names."""
    with path.open('rb') as f:
        header = f.read(32)
        count = struct.unpack_from('<I', header, 4)[0]
        length, row_length = struct.unpack_from('<HH', header, 8)
        fields = []
        while f.tell() < length - 1:
            descriptor = f.read(32)
            if len(descriptor) != 32 or descriptor[0] == 13:
                break
            fields.append((descriptor[:11].split(b'\x00')[0].decode('ascii'), descriptor[16]))
        f.seek(length)
        names = []
        for _ in range(count):
            record = f.read(row_length)
            offset = 1
            values = {}
            for key, size in fields:
                values[key] = record[offset:offset + size].decode('utf-8', errors='replace').strip(' \t\r\n\x00')
                offset += size
            names.append(values.get('ADMIN', values.get('NAME', '')))
        return names

@lru_cache(maxsize=3)
def polygons(shapefile):
    """Polygons."""
    file = Path(shapefile)
    names = _dbf_names(file.with_suffix('.dbf'))
    result = []
    with file.open('rb') as f:
        f.seek(100)
        index = 0
        while (header := f.read(8)):
            if len(header) != 8:
                break
            _, words = struct.unpack('>II', header)
            data = f.read(words * 2)
            kind = struct.unpack_from('<I', data)[0]
            name = names[index] if index < len(names) else ''
            index += 1
            if kind not in (5, 15, 25):
                continue
            parts, points = struct.unpack_from('<II', data, 36)
            starts = list(struct.unpack_from('<' + 'I' * parts, data, 44)) + [points]
            coords = np.frombuffer(data, dtype='<f8', count=points * 2, offset=44 + 4 * parts).reshape(-1, 2).copy()
            rings = [coords[a:b] for a, b in zip(starts[:-1], starts[1:]) if b - a >= 3]
            result.append((name, rings))
    return result

def draw_us(ax, sage_root):
    """Draw us."""
    ax.set_facecolor('#f5f5f5')
    file = Path(sage_root) / 'maps' / '50' / 'ne_50m_admin_0_countries.shp'
    if not file.exists():
        ax.text(0.5, 0.95, 'Natural Earth map file unavailable', transform=ax.transAxes, ha='center', color='#a05b25')
    else:
        for country, rings in polygons(str(file.resolve())):
            if country != 'United States of America':
                continue
            vertices = []
            codes = []
            for ring in rings:
                if np.ptp(ring[:, 0]) > 200 or ring[:, 0].max() < -126 or ring[:, 0].min() > -65 or (ring[:, 1].max() < 23) or (ring[:, 1].min() > 51):
                    continue
                vertices.extend(ring)
                codes.extend([MPath.MOVETO] + [MPath.LINETO] * (len(ring) - 2) + [MPath.CLOSEPOLY])
            if vertices:
                ax.add_patch(PathPatch(MPath(np.asarray(vertices), codes), facecolor='#e3e3e3', edgecolor='#222222', linewidth=1.0, zorder=0))
    ax.set_xlim(-127, -65)
    ax.set_ylim(20, 55)
    ax.set_aspect(1 / np.cos(np.deg2rad(38)), adjustable='box')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.grid(color='#acb9c4', alpha=0.3, linewidth=0.5)
    from matplotlib.ticker import FuncFormatter
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{abs(x):g}°W'))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f'{y:g}°N'))
