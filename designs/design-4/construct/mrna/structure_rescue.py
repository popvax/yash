"""Bounded synonymous repair of a source heuristic, preserving every hard constraint."""
import re
from collections import Counter

import RNA
from Bio.Seq import Seq


def longest(structure):
    return max((len(match.group()) for match in re.finditer(r"[()]+", structure)), default=0)


def rescue(coding, protein, allowed, frequencies, forbidden, protected, target=14, budget=120):
    initial = coding
    structure, _ = RNA.fold(coding)
    initial_run = longest(structure)
    evaluations, changes = 0, []
    for _ in range(3):
        structure, _ = RNA.fold(coding)
        current_run = longest(structure)
        if current_run <= target or evaluations >= budget:
            break
        pairing = RNA.ptable(structure)
        bases = set()
        for match in re.finditer(r"[()]+", structure):
            if len(match.group()) <= target:
                continue
            for position in range(match.start(), match.end()):
                bases.add(position)
                if pairing[position + 1]:
                    bases.add(pairing[position + 1] - 1)
        candidates = sorted({
            position // 3 for position in bases
            if not any(3 * (position // 3) + base in protected for base in range(3))
        })
        best, best_run = coding, current_run
        for position in candidates:
            old = coding[3 * position:3 * position + 3]
            alternatives = sorted(
                (codon for codon in allowed[protein[position]] if codon != old),
                key=lambda codon: (-frequencies[protein[position]][codon], codon),
            )
            for codon in alternatives:
                candidate = coding[:3 * position] + codon + coding[3 * position + 3:]
                if any(motif in candidate for motif in forbidden):
                    continue
                words = Counter(candidate[start:start + 12] for start in range(len(candidate) - 11))
                if any(count > 1 for count in words.values()):
                    continue
                folded, _ = RNA.fold(candidate)
                evaluations += 1
                run = longest(folded)
                if run < best_run:
                    best, best_run = candidate, run
                if run <= target or evaluations >= budget:
                    break
            if best_run <= target or evaluations >= budget:
                break
        if best == coding:
            break
        changes.extend(index for index, (old, new) in enumerate(zip(coding, best)) if old != new)
        coding = best
    final_run = longest(RNA.fold(coding)[0])
    if str(Seq(coding).translate()) != protein or any(coding[index] != initial[index] for index in protected):
        raise ValueError("Structure repair changed the protein or a frozen nucleotide")
    return coding, {
        "initial_paired_bracket_run_nt": initial_run,
        "final_paired_bracket_run_nt": final_run,
        "target_nt": target, "target_reached": final_run <= target,
        "candidate_folds_evaluated": evaluations, "budget": budget,
        "changed_nt_positions": changes,
        "scope": "A coding-only dot-bracket-run heuristic, not measured dsRNA, stability, expression or safety.",
    }
