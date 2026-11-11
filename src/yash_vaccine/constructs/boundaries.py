"""Position-aware junction queries, including the empty-linker case."""
LENGTHS = (8, 9, 10, 11)

def junction_peptides(left, linker, right):
    left, right = left[-10:], right[:10]
    sequence = left + linker + right
    begin, end = len(left), len(left) + len(linker)
    return [sequence[i:i+k] for k in LENGTHS for i in range(len(sequence)-k+1)
            if (i < end and i+k > begin if linker else i < begin < i+k)]

def start_peptides(sequence):
    return [] if sequence.startswith('M') else [('M'+sequence)[:k] for k in LENGTHS]

def all_queries(regions, linkers):
    queries = set()
    for i, row in enumerate(regions):
        sequence = row['sequence']
        queries.update(start_peptides(sequence))
        queries.update(sequence[p:p+k] for k in LENGTHS for p in range(len(sequence)-k+1))
        for j, other in enumerate(regions):
            if i != j:
                for linker in linkers:
                    queries.update(junction_peptides(sequence, linker['sequence'], other['sequence']))
    return sorted(queries)

def assemble(regions, order, chosen_linkers):
    """Coordinates follow native Vaxrank: zero based, end exclusive."""
    sequence = ''; parts = []
    def add(kind, name, amino_acids):
        nonlocal sequence
        start = len(sequence); sequence += amino_acids
        parts.append(dict(kind=kind, name=name, aa=amino_acids,
                          start_aa=start, end_aa=len(sequence)))
    if not regions[order[0]]['sequence'].startswith('M'):
        add('start_codon', 'M', 'M')
    for p, index in enumerate(order):
        row = regions[index]
        add('antigen', row['antigen_id'], row['sequence'])
        if p < len(chosen_linkers):
            linker = chosen_linkers[p]
            add('linker', linker['name'], linker['sequence'])
    return sequence, parts

def scan_boundaries(sequence, parts):
    """Independent scan of the complete chain; no pairwise enumerator used."""
    for k in LENGTHS:
        for pos in range(len(sequence)-k+1):
            end = pos+k
            touched = [p for p in parts if pos < p['end_aa'] and end > p['start_aa']
                       and p['end_aa'] > p['start_aa']]
            direct = any(p['kind']=='linker' and p['start_aa']==p['end_aa']
                         and pos < p['start_aa'] < end for p in parts)
            is_link = direct or any(p['kind']=='linker' for p in touched)
            is_start = any(p['kind']=='start_codon' for p in touched)
            if is_link or is_start:
                yield dict(start_aa=pos, end_aa=end, peptide=sequence[pos:end],
                           category='junction' if is_link else 'start_m')
