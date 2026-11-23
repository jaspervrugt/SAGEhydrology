# ==========================================================================
# SAGEhydrology | Python source
# File: backend/cartplus_basin_inputs.py
#
# Headless SAGE numerical backend.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Headless numerical routines extracted from cartplus_basin_inputs."""
from functools import lru_cache
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

@lru_cache(maxsize=3)
def _load_module(source):
    """Load module."""
    spec = spec_from_file_location('_cartplus_inputs', source)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def infiltration_floor(region):
    """Return Xia's fixed regional infiltration capacity floor in mm/d."""
    return 0.0001 if region == 'CAMELS_AU' else 0.0

def uses_elevation_bands(config):
    """Uses elevation bands."""
    return config.get('model') in (18, 19, 20)

def prepare_cartplus_bands(config, gauge, dat):
    """Attach elevation bands for the selected CaRT+ basin."""
    if not uses_elevation_bands(config):
        return
    region = config.get('region', 'CAMELS_US')
    inputs = _load_module(_source_root() / 'cartplus_inputs.py')
    validate_band_region(config)
    from xia_camels_us import gb_bands, project_region
    if region == 'CAMELS_GB' and project_region(config['data_root']) == region:
        import numpy as np
        dat['cartplus_elevation_bands'] = dict(elevation_m=gb_bands(config['data_root'], gauge), area_fraction=np.full(5, 0.2))
    elif region in inputs.ATTRIBUTE_REGIONS:
        import numpy as np
        from regional_maps import identifier
        rows = _attribute_rows(config['data_root'], region, config['sage_root'])
        row = rows.get(str(gauge)) or rows.get(identifier(gauge, region), {})
        elevation = inputs.bands_from_attributes(region, row)
        if elevation is None:
            raise ValueError(f'Missing usable basin elevation data for {region} basin {gauge}.')
        dat['cartplus_elevation_bands'] = dict(elevation_m=elevation, area_fraction=np.full(5, 0.2))
    else:
        inputs.attach_elevation_bands(dat, config['data_root'], gauge, region)
    dat['cartplus_infiltration_floor'] = infiltration_floor(region)

def validate_band_region(config, *, metadata_only=False):
    """Validate band region."""
    inputs = _load_module(_source_root() / 'cartplus_inputs.py')
    region = config.get('region', 'CAMELS_US')
    if metadata_only and region in inputs.SUBDAILY_REGIONS:
        return
    if config.get('resolution', 1) not in (1, 24, 96):
        raise ValueError('CaRT+ supports daily, hourly and 15-minute forcing.')
    if region in inputs.TERRAIN_REGIONS:
        root = config.get('data_root')
        if root and (Path(root) / 'cartplus_terrain_bands.csv').is_file():
            return
        raise ValueError(f'Terrain elevation bands are missing for {region}.')
    if region not in inputs.REGIONS and region not in inputs.ATTRIBUTE_REGIONS and (region not in inputs.SOURCES):
        raise ValueError(f'No CaRT+ elevation-band source is configured for {region}.')

@lru_cache(maxsize=8)
def _attribute_rows(data_root, region, sage_root):
    """Attribute rows."""
    from sage_attribute_reader import read
    return read(dict(data_root=data_root, region=region, sage_root=sage_root, resolution=1))
PARAMETER_LABELS = {'lC': ('\\ln C', 'Root-zone capacity', 'mm'), 'X4': ('X_{4}', 'Unit hydrograph time base', 'd'), 'lX3': ('\\ln X_{3}', 'Quick-store scale', 'mm'), 'ltau': ('\\ln \\tau', 'Groundwater residence time', 'd'), 'X2': ('X_{2}', 'Exchange loss coefficient', 'mm/d'), 'rKp': ('\\ln r_{\\mathrm{Kp}}', 'Relative percolation scale', '–'), 'TTs': ('T_{\\mathrm{snow}}', 'Rain–snow threshold', '°C'), 'lCF': ('\\ln c_{\\mathrm{m}}', 'Degree-day melt factor', 'mm/d/°C'), 'SFCF': ('f_{\\mathrm{sc}}', 'Snowfall correction factor', '–'), 'lGam': ('\\ln \\Gamma', 'Temperature lapse rate', '°C/km'), 'wz': ('w_{\\mathrm{z}}', 'Near-stream area fraction', '–'), 'lCz': ('\\ln C_{\\mathrm{z}}', 'Near-stream capacity', 'mm'), 'ltz': ('\\ln \\tau_{\\mathrm{z}}', 'Near-stream drainage time', 'd'), 'kv': ('\\kappa_{\\mathrm{v}}', 'Vegetation demand suppression', '–'), 'lamE': ('\\lambda_{\\mathrm{E}}', 'Groundwater uptake fraction', '–'), 'lbP': ('\\ln b', 'Root-zone Pareto shape', '–'), 'lCu': ('\\ln C_{\\mathrm{u}}', 'Surface layer capacity', 'mm'), 'ltu': ('\\ln \\tau_{\\mathrm{u}}', 'Surface drainage time', 'd'), 'lLc': ('\\ln L_{\\mathrm{c}}', 'Channel loss capacity', 'mm/d'), 'lCd': ('\\ln C_{\\mathrm{d}}', 'Mean depression capacity', 'mm'), 'lbd': ('\\ln b_{\\mathrm{d}}', 'Depression Pareto shape', '–'), 'lkd': ('\\ln k_{\\mathrm{d}}', 'Depression seepage coefficient', '1/d'), 'gR': ('\\gamma', 'Connected-area bypass fraction', '–'), 'lJin': ('\\ln J', 'External aquifer inflow', 'mm/d'), 'lks': ('\\ln \\kappa_{\\mathrm{s}}', 'Snow sublimation fraction', '–'), 'lrfa': ('\\ln r_{\\mathrm{p}}', 'Snowmelt per unit PET', '–'), 'lPs': ('\\ln p^{*}', 'Reference infiltration capacity', 'mm/d'), 'gie': ('g_{\\mathrm{F}}', 'Infiltration sensitivity to wetness', '–'), 'aF': ('a_{\\mathrm{F}}', 'Frost effect on quick release', '–'), 'k': ('k', 'Quick release shape exponent', '–'), 'Amin': ('A_{\\mathrm{min}}', 'Minimum connected fraction', '–'), 'lCt': ('\\ln C_{\\theta}', 'Wetness index capacity', 'mm'), 'ltt': ('\\ln \\tau_{\\theta}', 'Wetness index loss time', 'd'), 'lams': ('\\lambda', 'Root-zone connectivity effect', '–'), 'lamW': ('\\lambda_{\\mathrm{W}}', 'Groundwater connectivity effect', '–')}

def catalogue(resolution=1):
    """Catalogue."""
    if resolution not in (1, 24, 96) or not (_source_root() / 'parameters.json').is_file():
        return {}
    spec = _parameter_spec()
    parameters = [dict(symbol=symbol, description=description, unit=f'log of {unit}' if symbol.startswith('\\ln ') and unit != '–' else unit, lower=1.0 / 24 if name == 'X4' else bounds[0], upper=bounds[1]) for group in ('static', 'law') for name, bounds in spec[group].items() for symbol, description, unit in (PARAMETER_LABELS[name],)]
    result = {'18': dict(name='CaRT+', states=0, initial_states=[], parameters=parameters)}
    if (cart_root() / 'cart_daily.dll').is_file():
        from math import log
        bounds = dict(lC=(log(25), log(1600)), X4=(1.0 / 24, 6), lX3=(log(5), log(1000)), ltau=(log(3), log(1460)), X2=(-20, 0), rKp=(log(1 / 30), log(30)), TTs=(-2.5, 2.5), lCF=(log(0.5), log(10)), SFCF=(0.5, 1.5), lGam=(log(3), log(9)), lks=(log(0.0001), 0))
        result['19'] = dict(name='CaRT', states=0, initial_states=[], parameters=[dict(symbol=symbol, description=description, unit=f'log of {unit}' if symbol.startswith('\\ln ') and unit != '–' else unit, lower=lo, upper=hi) for name, (lo, hi) in bounds.items() for symbol, description, unit in (PARAMETER_LABELS[name],)])
    if False:
        from copy import deepcopy
        result['20'] = deepcopy(result['19'])
        result['20']['name'] = 'CaRTsplit'
        result['20']['parameters'].append(dict(symbol='f', description='UH1 / quick-store routing fraction', unit='-', lower=0.0, upper=1.0))
    return result

def cart_root():
    """Return the public baseline model source directory."""
    return _source_root().with_name('cart')

def _source_root():
    """Source root."""
    return Path(__file__).resolve().parents[3] / 'SAGEhydrology/Native/models/CaRTPlus'

@lru_cache(maxsize=1)
def _parameter_spec():
    """Parameter spec."""
    import json
    return json.loads((_source_root() / 'parameters.json').read_text())
