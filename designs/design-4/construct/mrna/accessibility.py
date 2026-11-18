"""Bounded synonymous search using ensemble accessibility, not expression prediction."""
import random
import re
from collections import Counter

import RNA
from Bio.Seq import Seq


def ensemble_accessibility(sequence, aug_start):
    compound = RNA.fold_compound(sequence)
    _, energy = compound.mfe()
    compound.exp_params_rescale(energy)
    compound.pf()
    pairs = compound.bpp()
    left, right = max(0, aug_start - 10), min(len(sequence), aug_start + 20)
    unpaired = []
    for index in range(left + 1, right + 1):
        paired = sum(pairs[other][index] for other in range(1, index))
        paired += sum(pairs[index][other] for other in range(index + 1, len(sequence) + 1))
        probability = 1.0 - paired
        if not -1e-6 <= probability <= 1 + 1e-6:
            raise ValueError("Invalid ensemble marginal unpaired probability")
        unpaired.append(min(1.0, max(0.0, probability)))
    atg = unpaired[aug_start - left:aug_start - left + 3]
    return {
        "window_start_0based": left, "window_end_exclusive": right,
        "AUG_start_0based": aug_start,
        "mean_marginal_unpaired_probability": sum(unpaired) / len(unpaired),
        "AUG_mean_marginal_unpaired_probability": sum(atg) / len(atg),
        "per_base_unpaired_probabilities": unpaired,
        "context_length_nt": len(sequence),
        "interpretation": "Mean marginal probabilities, not a joint probability that the initiation window is open, or a probability of translation.",
    }


def refine(coding, protein, allowed, forbidden, protected, utr5, suffix, seed, trials):
    generator = random.Random(seed + 2026)
    baseline = ensemble_accessibility(utr5 + coding + suffix, len(utr5))
    best, best_accessibility = coding, baseline
    initial_structure, _ = RNA.fold(coding)
    maximum_run = max((len(m.group()) for m in re.finditer(r"[()]+", initial_structure)), default=0)
    mutable = [
        position for position in range(6, min(45, len(protein)))
        if len(allowed[protein[position]]) > 1
        and not any(3 * position + base in protected for base in range(3))
    ]
    attempted = []
    for trial in range(trials):
        candidate = best
        if mutable:
            for position in generator.sample(mutable, min(len(mutable), 1 + trial % 3)):
                options = [codon for codon in allowed[protein[position]]
                           if codon != candidate[3 * position:3 * position + 3]]
                codon = generator.choice(options)
                candidate = candidate[:3 * position] + codon + candidate[3 * position + 3:]
        repeated = Counter(candidate[start:start + 12] for start in range(len(candidate) - 11))
        if any(motif in candidate for motif in forbidden) or any(count > 1 for count in repeated.values()):
            attempted.append({"trial": trial, "accepted": False, "reason": "motif_or_repeat"})
            continue
        structure, _ = RNA.fold(candidate)
        run = max((len(m.group()) for m in re.finditer(r"[()]+", structure)), default=0)
        if run > maximum_run:
            attempted.append({"trial": trial, "accepted": False, "reason": "paired_bracket_run_worsened", "run": run})
            continue
        accessibility = ensemble_accessibility(utr5 + candidate + suffix, len(utr5))
        score = accessibility["mean_marginal_unpaired_probability"]
        aug_not_worse = (
            accessibility["AUG_mean_marginal_unpaired_probability"] + 1e-9
            >= baseline["AUG_mean_marginal_unpaired_probability"]
        )
        changed = score > best_accessibility["mean_marginal_unpaired_probability"] + 1e-9 and aug_not_worse
        attempted.append({
            "trial": trial, "accepted": True, "selected_at_this_step": changed,
            "mean_marginal_unpaired_probability": score, "paired_bracket_run_nt": run,
            "AUG_not_worse_than_baseline": aug_not_worse,
        })
        if changed:
            best, best_accessibility = candidate, accessibility
    if str(Seq(best).translate()) != protein or any(best[i] != coding[i] for i in protected):
        raise ValueError("Accessibility refinement altered translation or a frozen nucleotide")
    return best, {
        "baseline": baseline, "selected": best_accessibility, "trials": attempted,
        "scope": "ViennaRNA partition-function ensemble of the known RNA segment before unspecified PolyA; standard unmodified RNA, 37 C.",
        "global_optimum_claimed": False, "expression_or_stability_validated": False,
        "AUG_marginal_nonworsening_constraint": True,
    }
