# ==========================================================================
# SAGEhydrology | Python source
# File: backend/sage_iteration_format.py
#
# Headless SAGE numerical backend.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Headless numerical routines extracted from sage_iteration_format."""
import math
FDC = {1: ('S_fdc', 'd_fdc'), 2: ('S_p', 'd_p'), 3: ('S_logp', 'd_logp')}

def number(value):
    """Number."""
    return float('nan') if value is None else float(value)

def fixed(value, decimals):
    """Fixed."""
    value = number(value)
    if math.isnan(value):
        return 'NaN'
    if math.isinf(value):
        return 'Inf' if value > 0 else '-Inf'
    return f'{value:.{decimals}f}'

def round_away(value, decimals=0):
    """Round away."""
    scale = 10 ** decimals
    return math.copysign(math.floor(abs(value) * scale + 0.5) / scale, value)

def metric(value):
    """Metric."""
    value = number(value)
    if not math.isfinite(value):
        return fixed(value, 3)
    decimals = 3 if abs(value) < 10 else max(0, 3 - math.floor(math.log10(abs(value))))
    rounded = abs(round_away(value, decimals))
    if rounded >= 10:
        decimals = max(0, 3 - math.floor(math.log10(rounded)))
    return fixed(value, decimals)

def sib(value):
    """Sib."""
    value = number(value)
    if math.isfinite(value) and abs(value) >= 100000.0:
        power = math.floor(math.log10(abs(value)))
        mantissa = round_away(value / 10 ** power, 1)
        if abs(mantissa) >= 10:
            mantissa /= 10
            power += 1
        return f'{mantissa:.1f}e{power}'
    return metric(value)

def scale(values):
    """Scale."""
    values = [abs(number(v)) for v in values if math.isfinite(number(v))]
    maximum = max(values, default=0)
    power = math.floor(math.log10(maximum)) if maximum >= 100 else 0
    return (power, 10 ** power)

def center(value, width, footprint=5):
    """Center."""
    return (' ' * max(0, (min(width, footprint) - len(value)) // 2) + value).ljust(width)[:width]

def line(values, widths, group_before):
    """Line."""
    parts = [values[0].ljust(widths[0])]
    for k, (value, width) in enumerate(zip(values[1:], widths[1:]), 2):
        if k <= 2:
            part = value.rjust(width)
        elif value in ('nse', 'kge', 'jkge', 'S_fdc', 'S_p', 'S_logp'):
            part = value.ljust(width) if len(value) > 5 else ' ' + center(value, width - 1)
        elif value == 'ΣGLS':
            part = '  ' + center(value, width - 2)
        else:
            try:
                positive = not math.isnan(float(value)) and (not value.startswith(('-', '+')))
            except ValueError:
                positive = False
            part = ((' ' if positive else '') + value).ljust(width)
        separator = '  | ' if k == 3 else '  ' if k == group_before else ' '
        parts.append(separator + part)
    return ''.join(parts) + '\n'

def groups(widths, prefix, nmetric, loss_power, rss_power, use_rss):
    """Groups."""
    starts = []
    stops = []
    pos = 0
    for k, width in enumerate(widths, 1):
        starts.append(pos)
        stops.append(pos + width - 1)
        gap = 4 if k == 2 else 2 if k == prefix + nmetric else 1
        pos += width + gap
    top = list(' ' * (stops[-1] + 1))
    lower = top.copy()

    def put(chars, start, text):
        """Put."""
        chars[start:start + len(text)] = text
    put(top, stops[1] - 3, 'loss')
    top[stops[1] + 3] = '|'
    lower[stops[1] + 3] = '|'
    if loss_power:
        put(lower, stops[1] - len(f'x10^{loss_power}') + 1, f'x10^{loss_power}')
    if use_rss and rss_power:
        label = f'x10^{rss_power}'
        put(lower, (starts[2] + stops[2] - len(label) + 1) // 2, label)
    left = starts[prefix] + 1
    middle = stops[prefix + nmetric - 1] - 1
    right = stops[-1]
    for a, b, label in ((left, middle, 'T_'), (starts[prefix + nmetric] + 1, right, 'S_ib')):
        lower[a:b + 1] = list('-' * (b - a + 1))
        lower[a] = lower[b] = '|'
        text = ' ' + label + ' '
        put(lower, (a + b - len(text) + 1) // 2, text)
    top[left:right + 1] = list('-' * (right - left + 1))
    top[left] = top[right] = '|'
    label = ' performance metrics '
    put(top, math.ceil((left + right - len(label) + 1) / 2), label)
    return ''.join(top).rstrip() + '\n' + ''.join(lower).rstrip() + '\n'

def loss_name(c, named=False):
    """Loss name."""
    form = c.get('fdc_formulation', 1) if c['loss'] == 6 else 1
    return {1: 'ΣSAR', 2: 'ΣGLS', 3: 'Σ1-NSE', 4: 'Σ1-KGE', 5: 'ΣHuber', 6: 'Σ' + FDC[form][1], 7: 'sum(1-JKGE)' if named else 'Σ1-JKGE'}[c['loss']]

def coefficient(c, name, count):
    """Coefficient."""
    weight = c.get('_observation_weights', {}).get(name, 1 / count)
    scaling = c.get('_observation_scales', {}).get(name, 1.0)
    value = weight / scaling
    if '_observation_n_total' in c:
        value *= c['_observation_n_total'] / c.get('_observation_counts', {}).get(name, c['_observation_n_total'])
    return (value, weight, scaling)

def compact(value, width, total=False):
    """Compact."""
    value = number(value)
    decimals = max(0, 4 - math.floor(math.log10(abs(value)))) if value and math.isfinite(value) else 4
    text = fixed(value, decimals if total else 3)
    if len(text) > width:
        text = f'{value:.1e}'
    if len(text) > width:
        text = f'{value:.0e}'
    return text

class Formatter:

    def __init__(self):
        """Init."""
        self.loss_power = self.rss_power = 0
        self.loss_scale = self.rss_scale = 1

    def format(self, event, c):
        """Format."""
        scenarios = event.get('scenarios') or [s for s in ('tt', 'te', 'et', 'ee') if event.get('diagnostics', {}).get(s, {}).get('observations', {}).get(c.get('observations', ['Q'])[0], {}).get('n', 0) > 0]
        if not scenarios:
            scenarios = ['tt']
        i = event['iteration']
        header = i == 1 or (i - 1) % 25 == 0
        leader = f"it {i:4d}   l {event.get('network_parameters', 0):6,d}   CPU {event.get('iteration_seconds', event.get('seconds', 0)):5.1f} s\n"
        names = c.get('observations', ['Q'])
        diagnostics = event.get('diagnostics', {})
        if names != ['Q']:
            return self.named(event, c, scenarios, leader, header)
        stats = [diagnostics.get(s, {}).get('observations', {}).get('Q', {}) for s in scenarios]
        fdc = FDC[c.get('fdc_formulation', 1) if c['loss'] == 6 else 1][0]
        jkge = c['loss'] == 7 and any((math.isfinite(number(r.get('jkge'))) for r in stats))
        use_rss = c['loss'] != 2
        if header:
            self.loss_power, self.loss_scale = scale([r.get('loss') for r in stats])
            self.rss_power, self.rss_scale = scale([r.get('rss') for r in stats])
        metrics = ['nse', 'kge'] + (['jkge'] if jkge else []) + [fdc]
        headings = ['scen', loss_name(c)] + (['ΣGLS'] if use_rss else []) + metrics + metrics
        widths = [6] + [7] * (len(headings) - 2) + [6]
        prefix = 2 + int(use_rss)
        before = prefix + len(metrics) + 1
        title = line(headings, widths, before)
        separator = '-' * (len(title) - 1) + '\n'
        block = '\n' + separator + leader + groups(widths, prefix, len(metrics), self.loss_power, self.rss_power, use_rss) + title + separator if header else '\n' + leader
        for scenario, r in zip(scenarios, stats):
            row = [scenario, fixed(number(r.get('loss')) / self.loss_scale, 3)]
            if use_rss:
                row.append(fixed(number(r.get('rss')) / self.rss_scale, 3))
            row += [metric(r.get(key.lower())) for key in metrics]
            row += [sib(r.get('sib_' + key.lower(), r.get('sib') if key == 'nse' else None)) for key in metrics]
            block += line(row, widths, before)
        return block

    def named(self, event, c, scenarios, leader, header):
        """Named."""
        names = c['observations']
        title = f"{'scen':<4} {'L_tot':>8}"
        group = ' ' * len(title)
        for name in names:
            columns = f"{'L':>7} {'NSE':>7} {'KGE':>7} {'JKGE':>7}"
            pad = len(columns) - len(name) - 2
            left = pad // 2
            group += '   ' + '-' * left + ' ' + name + ' ' + '-' * (pad - left)
            title += ' | ' + columns
        sep = '-' * len(title) + '\n'
        block = ''
        if header:
            detail = f"data types: {{{','.join(names)}}}, loss function: {loss_name(c, True)} (loss = {c['loss']})\n"
            formula = 'total loss: L_tot = '
            for k, name in enumerate(names):
                value, weight, scaling = coefficient(c, name, len(names))
                detail += f'  c_{name} = {value:.4g} (w = {weight:.4g}, scale = {scaling:.7g})\n'
                formula += (' + ' if k else '') + f'c_{name} x L_{name}'
            block = '\n' + sep + detail + formula + '\n' + sep + leader + group + '\n' + title + '\n' + sep
        else:
            block = leader
        for scenario in scenarios:
            summary = event.get('diagnostics', {}).get(scenario, {})
            row = f"{scenario:<4} {compact(summary.get('loss'), 8, True):>8}"
            for name in names:
                data = summary.get('observations', {}).get(name, {})
                component = f"{number(data.get('loss')):.4g}"
                if len(component) > 7:
                    component = f"{number(data.get('loss')):.1e}"
                if len(component) > 7:
                    component = f"{number(data.get('loss')):.0e}"
                row += ' | ' + f'{component:>7} ' + ' '.join((f'{compact(data.get(key), 7):>7}' for key in ('nse', 'kge', 'jkge')))
            block += row + '\n'
        return block + sep
