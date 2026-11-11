"""Adapt reference metadata after validating shared coordinates and marker alleles.

The original Broad map is preserved. We do not disable Picard's dictionary
validation or modify either BAM. A block with an incompatible marker is held
out whole, preserving the map's anchor relationships.
"""
from pathlib import Path
import collections
import hashlib
import json
import os
import tempfile

import pysam
from common import save,sha


def compatible_panel(root,tumour,normal):
    root=Path(root);source=root/'resources/Homo_sapiens_assembly38.haplotype_database.txt'
    dictionaries=[]
    for path in [tumour,normal]:
        with pysam.AlignmentFile(path,'rb') as bam:
            bam.check_index();dictionaries.append([(row['SN'],row['LN']) for row in bam.header.to_dict()['SQ']])
    assert dictionaries[0]==dictionaries[1],'Tumour/normal dictionaries differ; review coordinates rather than relabelling them'
    dictionary=dictionaries[0];lengths=dict(dictionary)
    alignment=json.loads((root/'alignment/complete.json').read_text())
    text=source.read_text();header=[];comments=[];rows=[]
    for line in text.splitlines():
        if line.startswith('@'):header.append(line);continue
        if not line or line.startswith('#'):comments.append(line);continue
        fields=line.split('\t');assert len(fields)>=7,line
        chrom,pos,name,major,minor=fields[:5];anchor=fields[6] or name
        rows.append(dict(chrom=chrom,pos=int(pos),name=name,major=major,minor=minor,anchor=anchor,line=line))
    assert len({row['name'] for row in rows})==len(rows)
    assert {row['anchor'] for row in rows}<={row['name'] for row in rows}
    bad=collections.defaultdict(list)
    with pysam.FastaFile(alignment['reference']) as reference:
        assert all(reference.get_reference_length(name)==length for name,length in dictionary)
        for row in rows:
            reason=None;chrom,pos=row['chrom'],row['pos']
            if chrom not in lengths:reason='CONTIG_NOT_IN_SHARED_BAM_DICTIONARY'
            elif not 1<=pos<=lengths[chrom]:reason='MARKER_OUTSIDE_CONTIG'
            elif len(row['major'])!=1 or len(row['minor'])!=1 or set(row['major']+row['minor'])-set('ACGT'):reason='NON_CANONICAL_MARKER_ALLELES'
            elif reference.fetch(chrom,pos-1,pos).upper() not in [row['major'],row['minor']]:reason='MARKER_ALLELES_DO_NOT_INCLUDE_ANALYSIS_REFERENCE_BASE'
            if reason:bad[row['anchor']].append(dict(marker=row['name'],reason=reason))
    retained=[row for row in rows if row['anchor'] not in bad]
    assert len(retained)>=100,'Too few coordinate-verified fingerprint markers; review the reference'
    by_chrom=collections.defaultdict(list)
    for row in retained:by_chrom[row['chrom']].append(row['line'])
    # Preserve all marker names, alleles, positions, frequencies and anchors.
    # Only dictionary metadata and whole explicitly incompatible blocks change.
    lines=[line for line in header if not line.startswith('@SQ')]
    lines+=['@SQ\tSN:'+name+'\tLN:'+str(length) for name,length in dictionary]
    lines+=comments
    lines += [line for chrom,_ in dictionary for line in by_chrom[chrom]]
    body='\n'.join(lines)+'\n';digest=hashlib.sha256(body.encode()).hexdigest()
    directory=root/'resources/fingerprint-coordinate-map';directory.mkdir(parents=True,exist_ok=True)
    panel=directory/('haplotypes-'+digest[:16]+'.txt')
    if panel.exists():assert sha(panel)==digest
    else:
        descriptor,temporary=tempfile.mkstemp(prefix='.haplotypes-',dir=directory)
        with os.fdopen(descriptor,'w') as stream:stream.write(body)
        os.replace(temporary,panel)
    audit=directory/'coordinate-audit.json'
    save(audit,dict(source=str(source),source_sha256=sha(source),adapted_panel=str(panel),adapted_panel_sha256=digest,
        source_header_contigs=sum(line.startswith('@SQ') for line in header),shared_BAM_dictionary=dictionary,
        original_markers=len(rows),retained_markers=len(retained),original_blocks=len({row['anchor'] for row in rows}),
        retained_blocks=len({row['anchor'] for row in retained}),excluded_whole_blocks=dict(bad),
        all_retained_marker_reference_alleles_verified=True,marker_alleles_positions_names_frequencies_unchanged=True,
        original_panel_modified=False,BAM_headers_modified=False,native_dictionary_validation_disabled=False,
        reference=alignment['reference'],recorded_reference_sha256=alignment['reference_sha256'],
        scope='Coordinate/marker compatibility with the shared analysis reference; original Caris mask sequence is unavailable'))
    return panel,audit
