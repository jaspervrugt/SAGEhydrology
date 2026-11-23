# ==========================================================================
# SAGEhydrology | Python source
# File: backend/regional_maps.py
#
# Headless SAGE numerical backend.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Headless numerical routines extracted from regional_maps."""
from pathlib import Path
import re
from io import BytesIO
import numpy as np
import pandas as pd
from geography import read_attributes, hydroclimatic_zones
from regional_data import _table

def identifier(value, region):
    """Identifier."""
    value = str(value).strip().removesuffix('.0')
    for prefix in ('GRDC_', 'hysets_', 'IL_', 'BULL_'):
        if value.lower().startswith(prefix.lower()):
            value = value[len(prefix):]
    if region in ('CAMELS_US', 'CAMELSH_US', 'MACH_US'):
        return value.zfill(8)
    if region == 'CAMELS_JP':
        return value.zfill(3)
    if region in ('CAMELS_CH', 'HYD_RESPONSES') and value.isdigit():
        return value.zfill(4)
    if region == 'CAMELS_LUX':
        return value.zfill(2)
    if region == 'CAMELS_CA':
        return value.removeprefix('CAN_')
    if region == 'CAMELS_CZ':
        return re.sub('^camelscz_', '', value, flags=re.I)
    if region == 'CAMELS_DE':
        return value.removeprefix('DE')
    return value

def attributes(root, region, sage):
    """Attributes."""
    import xia_camels_us
    if region == 'CAMELS_GB' and xia_camels_us.project_region(root) == region:
        return attributes(xia_camels_us.standard_root(root, region), region, sage)
    if region == 'CAMELS_US':
        if xia_camels_us.matches(root):
            rows = xia_camels_us.attributes(root)
            standard = xia_camels_us.standard_root(root, region)
            for gauge, geo in read_attributes(standard).items():
                if gauge in rows:
                    for key in ('gauge_lat', 'gauge_lon'):
                        if key in geo:
                            rows[gauge][key] = geo[key]
            return rows
        return read_attributes(root)
    root = Path(root)
    result = {}
    france_sites = {}
    france_station_sites = {}
    files = []
    for folder in (root, root / 'daily', root / 'hourly', root / 'attributes'):
        if not folder.exists():
            continue
        for file in folder.iterdir():
            name = file.name.lower()
            if file.suffix.lower() not in ('.csv', '.txt', '.xlsx'):
                continue
            if any((token in name for token in ('attribute', 'gauge_information', 'topograph', 'location', 'metadata', 'catchment', 'gauging', 'static_data', 'climate'))):
                files.append(file)
    files.sort(key=lambda f: 'gauge_information' in f.name.lower())
    aliases = {'gauge_lat': ['gauge_lat', 'gauge_lat_dd', 'latitude_wgs84', 'lat_gage', 'station_lat', 'latitude (wgs 84)', 'cwc_lat', 'lat', 'latitude'], 'gauge_lon': ['gauge_lon', 'gauge_lon_dd', 'longitude_wgs84', 'lon_gage', 'station_lon', 'longitude(wgs 84)', 'cwc_lon', 'lon', 'longitude'], 'elev_mean': ['elev_mean', 'gauge_elev', 'elev_gauge', 'elevation_mabsl', 'elevation', 'ele_mt_sav', 'elev'], 'aridity': ['aridity', 'aridity_era5l', 'aridity_index'], 'frac_snow': ['frac_snow', 'frac_snow_era5l', 'snow_fraction'], 'area': ['area_km2', 'area', 'area_calc', 'area_metadata', 'basin_area_km2', 'drain_sqkm', 'area_km2', 'uparea'], 'p_mean': ['p_mean', 'p_stn_sum', 'p_mean_mm', 'pre_mean', 'precip_mean', 'precipitation_mean', 'mean_precip', 'p_mean_annual', 'precip_mm_syr', 'annual_precipitation', 'pr0_mean', 'prec_mean', 'pmean_mm_year'], 'pet_mean': ['pet_mean', 'pet_pm_sum', 'et0_mean', 'pet_mean_mm', 'potential_evaporation_mean', 'pet_mean_annual', 'pet_mm_syr', 'eto_mean', 'mean_pet', 'annual_pet', 'pet1_mean', 'pet2_mean'], 't_mean': ['tmean_c', 't_mean', 'temp_mean', 'temperature_mean', 'mean_temperature', 'tmean', 'air_temperature_mean', 'temperature_2m_mean', 'tt_mean', 't_stn_mean']}
    aliases['gauge_lat'] += ['lat_outlet', 'lat_wgs84', 'ghi_lat', 'centroid_lat', 'gauge_latitude', 'latitude_deg', 'sta_y_w84', 'sta_y_w84_snap']
    aliases['gauge_lon'] += ['long_outlet', 'lon_wgs84', 'ghi_lon', 'centroid_lon', 'gauge_longitude', 'longitude_deg', 'sta_x_w84', 'sta_x_w84_snap']
    aliases['elev_mean'] += ['mean_ele', 'elev_med', 'z_mean', 'elv_mean', 'top_altitude_mean', 'sit_altitude', 'elevation_mean', 'mean_elevation', 'elev_m']
    aliases['aridity'] += ['aridity_catch', 'arid_1', 'arid_2', 'aridity_p_pet', 'ai_mean', 'aridity_ratio', 'pet_p_ratio', 'p_pet_pm', 'p_pet_ratio', 'aridity1_mean', 'aridity2_mean']
    aliases['frac_snow'] += ['snow_frac', 'fraction_snow', 'snowfall_fraction', 'frac_precip_snow', 'fracsnow1_mean', 'fracsnow2_mean']
    for file in files:
        keys = ['gauge_id', 'staid', 'id', 'basin_id', 'station_id', 'sta_code_h3', 'catchment_id', 'station_number', 'site_code', 'code_station', 'sit_code_h3', 'gauge_number', 'river_id', 'rid']
        wanted = set(keys + ['sit_longitude', 'sit_latitude'])
        wanted.update((c for candidates in aliases.values() for c in candidates))
        if region == 'CAMELS_CL' and file.name == '1_CAMELScl_attributes.txt':
            frame = _table(file).set_index('gauge_id').T.rename_axis('gauge_id').reset_index()
        elif file.suffix.lower() == '.xlsx':
            try:
                frame = pd.read_excel(file)
            except ImportError:
                from sage_xlsx import read_sheet
                frame = read_sheet(file)
        else:
            from local_data import attribute_bytes
            raw = attribute_bytes(file)
            first = raw.splitlines()[0].decode('utf-8-sig', errors='replace') if raw else ''
            delimiter = ';' if first.count(';') > first.count(',') else ',' if ',' in first else '\t' if '\t' in first else '\\s+'
            options = dict(sep=delimiter, usecols=lambda c: str(c).strip().lower() in wanted)
            try:
                frame = pd.read_csv(BytesIO(raw), encoding='utf-8-sig', **options)
            except UnicodeError:
                frame = pd.read_csv(BytesIO(raw), encoding='cp1252', **options)
        lookup = {str(c).strip().lower(): c for c in frame.columns}
        if region == 'CAMELS_NZ':
            keys = ['rid'] + keys
        key = next((lookup[k] for k in keys if k in lookup), None)
        if key is None:
            continue
        projection = None
        if region in ('LAMAH_CE', 'LAMAH_CE', 'LAMAH_AT', 'CAMELS_IS') and file.name == 'Gauge_attributes.csv':
            from pyproj import Transformer
            projection = Transformer.from_crs(3035 if region != 'CAMELS_IS' else 3057, 4326, always_xy=True)
        if region == 'CAMELS_FR' and all((k in lookup for k in ('sit_longitude', 'sit_latitude'))):
            from pyproj import Transformer
            projection = Transformer.from_crs(2154, 4326, always_xy=True)
        for _, row in frame.iterrows():
            if pd.isna(row[key]):
                continue
            gauge = identifier(row[key], region)
            if region == 'CAMELS_FR' and key == 'sit_code_h3':
                record = france_sites.setdefault(gauge, {})
            else:
                record = result.setdefault(gauge, {})
                if region == 'CAMELS_FR' and 'sit_code_h3' in lookup:
                    site = row[lookup['sit_code_h3']]
                    if pd.notna(site):
                        france_station_sites[gauge] = identifier(site, region)
            for target, candidates in aliases.items():
                col = next((lookup[c] for c in candidates if c in lookup), None)
                if col is None:
                    continue
                value = pd.to_numeric(row[col], errors='coerce')
                if np.isfinite(value):
                    record[target] = float(value)
            if region == 'CAMELS_PL' and 'topographic' in file.name.lower() and all((k in lookup for k in ('gauge_lon', 'gauge_lat'))):
                record['gauge_lat'] = float(row[lookup['gauge_lon']])
                record['gauge_lon'] = float(row[lookup['gauge_lat']])
            if projection:
                x, y = ('sit_longitude', 'sit_latitude') if region == 'CAMELS_FR' else ('lon', 'lat')
                lon, lat = projection.transform(float(row[lookup[x]]), float(row[lookup[y]]))
                record.update(gauge_lon=lon, gauge_lat=lat)
    for gauge, site in france_station_sites.items():
        for key, value in france_sites.get(site, {}).items():
            result[gauge].setdefault(key, value)
    return {g: r for g, r in result.items() if abs(r.get('gauge_lat', np.inf)) <= 90 and abs(r.get('gauge_lon', np.inf)) <= 180}

def zones(records, region, minimum, maximum):
    """
    Apply the explicit regional fallbacks in
    classify_camels_all_zones.m.

    Derived aridity values below are classification codes, not PET/P
    measurements.
    Records and the downloaded metadata are never modified.
    """
    prepared = {g: dict(r) for g, r in records.items()}
    for r in prepared.values():
        if 'aridity' not in r and r.get('p_mean', 0) > 0 and ('pet_mean' in r):
            r['aridity'] = r['pet_mean'] / r['p_mean']
        if 'aridity' not in r and region == 'CAMELS_PL' and ('p_mean' in r):
            r['aridity'] = 2 if r['p_mean'] < 1.6 else 1 if r['p_mean'] < 2.3 else 0.5
        if 'aridity' not in r and region in ('CAMELS_FR', 'CAMELS_SE', 'CAMELS_ET'):
            r['aridity'] = 1
        if 'frac_snow' not in r and 't_mean' in r:
            r['frac_snow'] = np.clip((2 - r['t_mean']) / 8, 0, 1)
        if 'frac_snow' not in r and region == 'CAMELS_ET':
            r['frac_snow'] = 0
    if region in ('CAMELS_DE', 'CAMELS_DEH'):
        values = [r['p_mean'] for r in prepared.values() if 'aridity' not in r and 'p_mean' in r]
        if values:
            q = np.quantile(values, [1 / 3, 2 / 3], method='hazen')
            for r in prepared.values():
                if 'aridity' not in r and 'p_mean' in r:
                    r['aridity'] = 2 if r['p_mean'] <= q[0] else 1.3 if r['p_mean'] <= q[1] else 0.8
    return hydroclimatic_zones(prepared, minimum, maximum)
