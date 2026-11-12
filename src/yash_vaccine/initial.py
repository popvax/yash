"""Reproduce every stage of Design #1 using included evidence."""
from yash_vaccine.paths import resolve_input
import importlib.metadata
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from .io import sha256, write_json


def reproduce(repo, out=None):
    repo = Path(repo).resolve()
    reference = repo / "designs/initial"
    manifest = json.loads((resolve_input(reference / "candidate.json")).read_text())
    from .provenance import RecordedSource
    recorded = RecordedSource(repo)
    recorded.verify_files({str((reference / name).resolve().relative_to(repo)): {"sha256": digest}
                           for name, digest in manifest["artifact_sha256"].items()})
    for line in (resolve_input(reference / "construct/mrna/requirements.txt")).read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        package, expected = line.split("==")
        try:
            installed = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            installed = None
        if installed != expected:
            raise ValueError("Complete reproduction needs the pinned RNA dependencies. "
                             "Install with: python -m pip install -e '.[mrna]'")

    def execute(destination):
        from .designs import reproduce as reproduce_recipe

        destination = Path(destination).resolve()
        if resolve_input(destination).exists():
            raise ValueError("Choose a new reproduction output directory: " + str(destination))
        destination.mkdir(parents=True)
        # This stage rebuilds selection and the optional internal-cassette checkpoint.
        selection = reproduce_recipe(repo, reference / "selection", destination / "selection")
        construct = destination / "construct"
        construct.mkdir()

        def run(script, *arguments):
            result = subprocess.run([sys.executable, str(script), *map(str, arguments)],
                                    cwd=repo, text=True, capture_output=True)
            if result.returncode:
                raise ValueError("Reproduction failed in " + str(script)
                                 + "\n" + result.stdout + result.stderr)

        # The protein receipts hash their generating code. Replay that exact
        # recorded runner in a temporary input tree, rather than substituting
        # current path-maintenance code into historical provenance fields.
        with tempfile.TemporaryDirectory(prefix="yash-initial-recorded-source-") as temporary:
            source = Path(temporary)
            needed = {}
            for name, digest in manifest["artifact_sha256"].items():
                include = (name.startswith("construct/")
                           and not name.startswith(("construct/protein/", "construct/mrna/")))
                include |= name in {
                    "selection/selected.tsv", "selection/verification.json", "selection/design.json",
                    "../../config/case.json", "../../config/linkers.json",
                    "../../data/checkpoints/manifest.json",
                    "../../data/checkpoints/provenance/construct-inference.json",
                    "../../data/checkpoints/provenance/construct-model-files.json",
                    "../../src/yash_vaccine/hydrophobicity.py",
                }
                if include:
                    needed[str((reference / name).resolve().relative_to(repo))] = {"sha256": digest}
            cache_manifest = json.loads(resolve_input(repo / "data/checkpoints/manifest.json").read_text())
            needed["data/checkpoints/construct-predictions.npz"] = cache_manifest["construct-predictions.npz"]
            recorded.materialize(source, needed)
            runner = source / "designs/initial/construct"
            run(runner / "run.py", "--output", construct / "protein")
            run(runner / "verify.py", "--protein", construct / "protein",
                "--output", construct / "final-check.json")

        # Run the RNA builder against the freshly reproduced protein and check.
        mrna = construct / "mrna"
        mrna.mkdir()
        for name in ("run.py", "genbank.py", "annotation-sources.json", "backbone.json",
                     "human-codon-usage.json", "requirements.txt", "optimizer.py"):
            shutil.copy2(resolve_input(reference / "construct/mrna" / name), mrna / name)
        shutil.copy2(resolve_input(reference / "construct/architecture.json"), construct / "architecture.json")
        run(mrna / "run.py", "--output", mrna)

        # Compare the freshly generated outputs, including detailed checks.
        checked = []
        for name, digest in manifest["artifact_sha256"].items():
            generated = name.startswith("construct/protein/") or name == "construct/final-check.json"
            generated |= name in manifest["nucleotide_outputs"]
            if generated:
                if sha256(destination / name) != digest:
                    raise ValueError("Reproduced Design #1 artifact differs: " + name)
                checked.append(name)
        if not checked:
            raise ValueError("No complete-design outputs were verified")
        junction = json.loads((resolve_input(construct / "final-check.json")).read_text())
        nucleotide = json.loads((resolve_input(mrna / "verification.json")).read_text())
        solver = json.loads((resolve_input(construct / "protein/construct.json")).read_text())["solver"]
        if (junction["status"] != "FRESH_FULL_JUNCTION_CHECK_PASSED"
                or not nucleotide["all_constraints_pass"]
                or nucleotide["polyA_status"] != "unspecified"
                or nucleotide["complete_construct_verified"]
                or not solver["all_phases_proven_optimal"]):
            raise ValueError("Incomplete verification of the final design")
        result = dict(status="VERIFIED", name="initial", selected_regions=selection["selected_regions"],
                      protein_length_aa=junction["protein_length_aa"],
                      rna_length_nt=nucleotide["lengths"]["transcribed_rna_nt"],
                      template_length_bp=nucleotide["lengths"]["supplied_template_context_bp"],
                      polyA_status=nucleotide["polyA_status"],
                      complete_RNA_template_sequence=False,
                      known_rna_segments_nt=nucleotide["lengths"]["known_rna_segments_nt"],
                      all_recorded_artifacts_identical=True, all_solver_phases_optimal=True,
                      hydrophobicity_pass=junction['hydrophobicity']['status']=='PASS',
                      stages=["selection", "hydrophobicity_screen", "joint_SEC_antigens_MITD", "junction_verification", "RNA_encoding"],
                      junction_prediction_pairs=junction["unique_peptide_hla_pairs"],
                      model_inference_rerun=False, cloud_jobs_started=False,
                      protein_runner="Verified recorded source in a temporary input tree",
                      replay_code_sha256={"src/yash_vaccine/initial.py": sha256(Path(__file__))},
                      source_snapshot_verification=recorded.verification_scope(),
                      output=str(destination) if out else None)
        write_json(destination / "verification.json", result)
        return result

    if out is not None:
        return execute(out)
    with tempfile.TemporaryDirectory(prefix="yash-complete-initial-") as temporary:
        return execute(Path(temporary) / "initial")
