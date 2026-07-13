import os
import math
from collections import namedtuple

# Rich result of a single-channel MC-TST rate-constant calculation.
RateResult = namedtuple('RateResult',
    ['k', 'kappa', 'Ea', 'Q_TS', 'Q_reactant', 'sigma', 'n_ts', 'n_reactant', 'imag', 'T'])
RateResult.__new__.__defaults__ = (float('nan'), float('nan'), None, None, None, 1, 0, 0, None, 298.15)


K_UNITS = 'cm3 molecule-1 s-1'


def format_rate(result):
    if result.k is None or (isinstance(result.k, float) and math.isnan(result.k)):
        return f"k({result.T} K) could not be computed (missing energies)"
    return (f"k({result.T} K) = {result.k:.3e} {K_UNITS}  "
            f"(kappa = {result.kappa:.2f}, Ea = {result.Ea:.2f} kcal/mol, sigma = {result.sigma})")


def _num(x, fmt):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return '-'
    return format(x, fmt)


def _formula(channel_name):
    return channel_name.split('_H')[0]


def _channel_tag(channel_name):
    i = channel_name.rfind('_H')
    return channel_name[i + 1:] if i != -1 else channel_name


def record_rate(rates_dir, channel_name, result, method=None):
    store = os.path.join(rates_dir, '.rates.tsv')
    rows = {}
    if os.path.exists(store):
        with open(store) as f:
            for line in f:
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 7:
                    rows[parts[0]] = parts[1:]
    rows[channel_name] = [method or '', str(result.sigma), _num(result.Ea, '.4f'),
                          _num(result.kappa, '.4f'), repr(result.k), _num(result.Q_TS, '.6e'),
                          str(result.T)]
    with open(store, 'w') as f:
        for ch, vals in rows.items():
            f.write('\t'.join([ch] + vals) + '\n')
    _write_summary(os.path.join(rates_dir, 'Rate_constants.txt'), rows)


def _write_summary(path, rows):
    def hnum(ch):
        tag = _channel_tag(ch)
        return int(tag[1:]) if tag[1:].isdigit() else 0
    channels = sorted(rows, key=hnum)
    if not channels:
        return
    formula = _formula(channels[0])
    method = rows[channels[0]][0]
    T = rows[channels[0]][6]
    total = 0.0
    for ch in channels:
        try:
            k = float(rows[ch][4])
            if not math.isnan(k):
                total += k
        except ValueError:
            pass
    W = 80
    head = f' {formula} + OH' + (f'   ({method})' if method else '')
    out = ['=' * W, f'{head:<{W - 15}}T = {T} K', '=' * W,
           f" {'Channel':<9}{'sigma':>7}{'Ea/kcal/mol':>14}{'kappa':>10}{'k / ' + K_UNITS:>28}",
           ' ' + '-' * (W - 2)]
    for ch in channels:
        _, sig, ea, kap, k, q, _t = rows[ch]
        try:
            kv = float(k)
            kf = 'failed' if math.isnan(kv) else f'{kv:.3e}'
        except ValueError:
            kf = 'failed'
        out.append(f" {_channel_tag(ch):<9}{sig:>7}{ea:>14}{kap:>10}{kf:>28}")
    out += [' ' + '-' * (W - 2),
            f' TOTAL  k({T} K) = {total:.3e} {K_UNITS}', '=' * W]
    with open(path, 'w') as f:
        f.write('\n'.join(out) + '\n')


def write_molecule_summary(path, molecules, title=None):
    W = 96
    out = ['=' * W, ' Molecule summary' + (f'  ({title})' if title else ''), '=' * W,
           f" {'Name':<34}{'Step':<12}{'E_elec / Ha':>16}{'ZPE / Ha':>14}{'Q':>16}",
           ' ' + '-' * (W - 2)]
    for m in molecules:
        out.append(f" {(m.name or '')[:33]:<34}{str(getattr(m, 'current_step', '') or '')[:11]:<12}"
                   f"{_num(getattr(m, 'electronic_energy', None), '.6f'):>16}"
                   f"{_num(getattr(m, 'zero_point', None), '.6f'):>14}"
                   f"{_num(getattr(m, 'Q', None), '.3e'):>16}")
    out.append('=' * W)
    with open(path, 'w') as f:
        f.write('\n'.join(out) + '\n')
