"""Per-allele affinity/percentile calibration built only from recorded predictions.

MHCflurry reports a predicted affinity in nM and a percentile rank of that
affinity among random peptides for the same allele. The nM scale is not
comparable between alleles; the percentile is. These helpers recover the
nM <-> percentile relation of the pinned model bundle from prediction tables
already committed to this repository, so no model or new inference is needed.
"""
from pathlib import Path
import bisect
import collections
import csv
import functools
import gzip
import math

ALLELES = ('HLA-A*01:01', 'HLA-A*33:03', 'HLA-B*08:01',
           'HLA-B*15:08', 'HLA-C*01:02', 'HLA-C*07:02')

# Every committed table that carries an affinity and its percentile for the
# six case alleles. All were produced by the same pinned MHCflurry 2.3 bundle.
CALIBRATION_SOURCES = (
    'data/shared-rna/splicing/scoring/affinities.tsv.gz',
    'designs/initial/construct/terminal-inference/predictions.tsv',
    'designs/initial/construct/protein/boundary-predictions.tsv',
    'designs/initial/construct/selected-inference/predictions.tsv',
    'designs/initial/selection/boundary-predictions.tsv',
    'designs/initial/inputs/mutations/rank/peptide-HLA-affinity.tsv',
    'designs/design-2/construct/protein/boundary-predictions.tsv',
    'designs/design-2/construct/selected-inference/predictions.tsv',
    'designs/design-3/construct/protein/boundary-predictions.tsv',
    'designs/design-3/construct/selected-inference/predictions.tsv',
    'designs/design-3/inputs/predictions/peptide-HLA-affinity.tsv',
    'designs/design-3/review/additional-mutation-affinities.tsv',
    'designs/design-3/review/candidate-annotations/predictions.tsv.gz',
)

PERCENTILE_FIELDS = ('affinity_percentile', 'percentile_rank')

class _Curve:
    """Sorted (key, value) pairs held as parallel lists for fast bisection."""

    __slots__ = ('keys', 'values')

    def __init__(self, pairs):
        self.keys = [pair[0] for pair in pairs]
        self.values = [pair[1] for pair in pairs]

    def __len__(self):
        return len(self.keys)


def open_text(path):
    path = Path(path)
    return gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()


def read_rows(path):
    with open_text(path) as stream:
        yield from csv.DictReader(stream, delimiter='\t')


def _percentile_of_row(row):
    for field in PERCENTILE_FIELDS:
        if row.get(field) not in (None, ''):
            return float(row[field])
    return None


class Calibration:
    """Monotone nM <-> percentile curves per allele, plus an allele-neutral scale."""

    def __init__(self, repo, sources=CALIBRATION_SOURCES, alleles=ALLELES):
        self.alleles = tuple(alleles)
        self.observed = {}          # (peptide, allele) -> percentile
        self.observed_nm = {}       # (peptide, allele) -> affinity in nM
        self.sources = tuple(sources)
        points = {allele: {} for allele in self.alleles}
        lengths = {allele: {} for allele in self.alleles}
        raw = {allele: [] for allele in self.alleles}
        for relative in self.sources:
            for row in read_rows(Path(repo) / relative):
                allele = row.get('allele')
                if allele not in points:
                    continue
                percentile = _percentile_of_row(row)
                if percentile is None or row.get('affinity_nm') in (None, ''):
                    continue
                affinity = float(row['affinity_nm'])
                points[allele][affinity] = percentile
                raw[allele].append((affinity, percentile))
                peptide = row.get('peptide', '')
                if peptide:
                    self.observed[(peptide, allele)] = percentile
                    self.observed_nm[(peptide, allele)] = affinity
                    lengths[allele].setdefault(len(peptide), {})[affinity] = percentile
        self.by_affinity = {a: _Curve(sorted(d.items())) for a, d in points.items()}
        self.by_percentile = {a: _Curve(sorted((p, nm) for nm, p in d.items()))
                              for a, d in points.items()}
        self.by_length = {a: {n: _Curve(sorted(d.items())) for n, d in lengths[a].items()}
                          for a in self.alleles}
        self.by_length_percentile = {
            a: {n: _Curve(sorted((p, nm) for nm, p in d.items()))
                for n, d in lengths[a].items()} for a in self.alleles}
        self.raw_rows = {a: sorted(v) for a, v in raw.items()}
        for allele in self.alleles:
            if len(self.by_affinity[allele]) < 1000:
                raise ValueError('too few calibration points for ' + allele)

    def inversions(self, tolerance=1e-9):
        """Raw rows where a higher affinity carries a lower percentile.

        Counted on every recorded row, not on the affinity-keyed curve, so that
        two rows sharing an affinity but disagreeing on percentile are visible
        rather than silently collapsed by the de-duplication in __init__.
        """
        total = 0
        for rows in self.raw_rows.values():
            for index in range(1, len(rows)):
                if rows[index][1] < rows[index - 1][1] - tolerance:
                    total += 1
        return total

    def raw_row_count(self):
        return {allele: len(rows) for allele, rows in self.raw_rows.items()}

    def tied_affinity_disagreements(self, tolerance=1e-9):
        """Affinities recorded more than once with a differing percentile.

        inversions() sorts by (affinity, percentile), so a disagreement inside
        one affinity can never present as an inversion there. This is the test
        that the pooled points really define a single curve: it counts, per
        allele, the affinities whose recorded percentiles are not all equal.
        """
        counts = {}
        for allele, rows in self.raw_rows.items():
            grouped = collections.defaultdict(list)
            for affinity, percentile in rows:
                grouped[affinity].append(percentile)
            counts[allele] = sum(1 for values in grouped.values()
                                 if max(values) - min(values) > tolerance)
        return counts

    def worst_tied_affinity_spread(self):
        """Largest relative percentile spread within a single affinity."""
        worst = 0.0
        for rows in self.raw_rows.values():
            grouped = collections.defaultdict(list)
            for affinity, percentile in rows:
                grouped[affinity].append(percentile)
            for values in grouped.values():
                low, high = min(values), max(values)
                if low > 0:
                    worst = max(worst, (high - low) / low)
        return worst

    @staticmethod
    def _interpolate(curve, value):
        keys, values = curve.keys, curve.values
        index = bisect.bisect_left(keys, value)
        if index == 0:
            return values[0]
        if index >= len(keys):
            return values[-1]
        low_key, high_key = keys[index - 1], keys[index]
        low, high = values[index - 1], values[index]
        span = math.log(max(high_key, 1e-12)) - math.log(max(low_key, 1e-12))
        if span == 0:
            return low
        fraction = (math.log(max(value, 1e-12)) - math.log(max(low_key, 1e-12))) / span
        logged = (math.log(max(low, 1e-12))
                  + fraction * (math.log(max(high, 1e-12)) - math.log(max(low, 1e-12))))
        return math.exp(logged)

    def percentile(self, allele, affinity_nm, peptide=None, length=None):
        """Recorded percentile when available, otherwise interpolated."""
        if peptide is not None:
            found = self.observed.get((peptide, allele))
            if found is not None:
                return found
        if length is not None and length in self.by_length[allele]:
            return self._interpolate(self.by_length[allele][length], affinity_nm)
        return self._interpolate(self.by_affinity[allele], affinity_nm)

    def affinity_at(self, allele, percentile, length=None):
        if length is not None and length in self.by_length_percentile[allele]:
            return self._interpolate(self.by_length_percentile[allele][length], percentile)
        return self._interpolate(self.by_percentile[allele], percentile)

    @functools.lru_cache(maxsize=None)
    def neutral_affinity(self, percentile):
        """Geometric mean, across the six alleles, of the nM at this percentile.

        Replacing an allele's nM with this value keeps the percentile of every
        peptide and removes only the per-allele offset of the nM scale.
        """
        logs = [math.log(max(self.affinity_at(a, percentile), 1e-12)) for a in self.alleles]
        return math.exp(sum(logs) / len(logs))
