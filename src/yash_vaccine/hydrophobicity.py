"""Conservative membrane-like hydropathy screen, separate from somatic evidence.

Uniform 19-aa Kyte–Doolittle windows; exclude means strictly greater than 1.6.
Integer tenths make the cutoff exact. This is a heuristic, not a topology model.
"""
from collections import Counter
from .io import write_json, write_tsv

POLICY_ID = 'kyte_doolittle_19_gt_1_6'
WINDOW = 19
KD_TENTHS = dict(zip('ARNDCQEGHILKMFPSTWYV',
                    [18,-45,-35,-35,25,-35,-35,-4,-32,45,38,-39,19,28,-16,-8,-7,-9,-13,42]))
POLICY = dict(id=POLICY_ID, scale='Kyte–Doolittle (1982)', window_aa=WINDOW,
              weighting='uniform; no edge padding', exclude_mean_strictly_greater_than=1.6,
              coordinates='zero based, end exclusive',
              scale_source='https://web.expasy.org/protscale/pscale/Hydropath.Doolittle.html',
              method_source='https://doi.org/10.1016/0022-2836(82)90515-0',
              scope='Construct-compatibility heuristic, not proof of membrane insertion or its absence.')


def scan(sequence):
    """All complete windows, without rounding scores before the decision."""
    if not isinstance(sequence, str) or not sequence or set(sequence)-KD_TENTHS.keys():
        raise ValueError('Hydrophobicity requires a nonempty, uppercase, standard-amino-acid sequence')
    if len(sequence) < WINDOW:
        raise ValueError('Hydrophobicity requires at least 19 residues; no padded or partial windows')
    values = [KD_TENTHS[a] for a in sequence]
    return [dict(start_aa=i, end_aa=i+WINDOW, sequence=sequence[i:i+WINDOW],
                 mean_kd=sum(values[i:i+WINDOW])/(10*WINDOW),
                 above_threshold=sum(values[i:i+WINDOW]) > 16*WINDOW)
            for i in range(len(sequence)-WINDOW+1)]


def summarize(sequence):
    windows = scan(sequence)
    best = max(windows, key=lambda row: row['mean_kd'])  # first start breaks ties
    return dict(status='EXCLUDE_HYDROPHOBIC' if any(w['above_threshold'] for w in windows) else 'PASS',
                max_mean_kd=best['mean_kd'], max_start_aa=best['start_aa'],
                max_end_aa=best['end_aa'], max_window=best['sequence'],
                windows_scanned=len(windows), above_threshold_windows=sum(w['above_threshold'] for w in windows))


def screen_pool(pool, out):
    """Screen every nominated region before source quotas, retaining original evidence."""
    audit = [dict(source=r['source'], region_id=r['region_id'], event_id=r['event_id'],
                  gene=r['gene'], sequence=r['sequence'], **summarize(r['sequence'])) for r in pool]
    write_tsv(out/'hydrophobicity-screen.tsv', audit)
    write_json(out/'hydrophobicity-summary.json', dict(policy=POLICY, regions=len(audit),
        counts=dict(Counter(r['status'] for r in audit)),
        excluded=[{k:r[k] for k in ('source','region_id','gene','max_mean_kd','max_window')}
                  for r in audit if r['status'] != 'PASS'],
        window_selection_scope='Existing source-nominated 25-aa regions; no automatic replacement-window search.',
        changes_DNA_or_RNA_evidence=False))
    return {r['region_id']:r for r in audit}


def protein_screen(sequence, parts):
    """Scan the entire chain, while enforcing the gate outside routing elements.

Every unmasked window is exported. Windows wholly in SEC/MITD and windows
crossing their boundaries are distinguished. Boundary overlaps are reported for
review, never asserted to have experimentally validated topology. Antigens are
also screened individually, so touching routing cannot hide a hydrophobic antigen.
    """
    cursor = 0
    routing = []
    for part in parts:
        a,b = int(part['start_aa']),int(part['end_aa'])
        if a != cursor or sequence[a:b] != part['aa']:
            raise ValueError('Hydrophobicity component coordinates do not match the full protein')
        cursor = b
        if part['kind'] == 'antigen' and summarize(part['aa'])['status'] != 'PASS':
            raise ValueError('Hydrophobic antigen reached assembly: ' + part['name'])
        if part['kind'] in ('signal_peptide','mitd'):
            if (part['kind'] == 'signal_peptide' and a != 0) or (part['kind'] == 'mitd' and b != len(sequence)):
                raise ValueError('Only the fixed terminal routing elements may be exempted')
            routing.append((a,b,part['kind']))
    if cursor != len(sequence) or len({p[2] for p in routing}) != len(routing):
        raise ValueError('Invalid protein tiling or duplicated routing domains')
    windows = scan(sequence)
    for row in windows:
        a,b = row['start_aa'], row['end_aa']
        touched = [kind for x,y,kind in routing if a < y and b > x]
        contained = any(x <= a and b <= y for x,y,_ in routing)
        row.update(classification='INTENTIONAL_ROUTING_DOMAIN' if contained else
                   'ROUTING_BOUNDARY_REVIEW' if touched else 'ANTIGEN_LINKER_CASSETTE',
                   routing_overlap=';'.join(touched))
    flagged = [w for w in windows if w['above_threshold']]
    internal = [w for w in windows if w['classification'] == 'ANTIGEN_LINKER_CASSETTE']
    failures = [w for w in internal if w['above_threshold']]
    return windows, dict(policy=POLICY, status='PASS' if not failures else 'FAIL',
        protein_length_aa=len(sequence), windows_scanned=len(windows),
        antigen_linker_windows=len(internal), unexpected_high_windows=len(failures),
        max_antigen_linker_mean_kd=max((w['mean_kd'] for w in internal), default=None),
        high_window_classifications=dict(Counter(w['classification'] for w in flagged)),
        routing_boundary_policy='Report all unmasked scores; SEC/MITD overlaps are not excluded. This does not validate boundary topology.',
        all_antigen_regions_individually_pass=True)


def write_protein_screen(out, sequence, parts):
    windows, report = protein_screen(sequence, parts)
    write_tsv(out/'hydrophobicity-windows.tsv', windows)
    write_json(out/'hydrophobicity.json', report)
    if report['status'] != 'PASS':
        raise ValueError('Unexpected high-hydropathy windows in the assembled antigen/linker cassette')
    return report
