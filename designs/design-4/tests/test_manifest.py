import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "pipeline")]

from dataio import save, sha
from verify_results import verify


class BundleManifestTests(unittest.TestCase):
    def test_deleting_manifest_cannot_disable_the_integrity_gate(self):
        with tempfile.TemporaryDirectory() as temporary:
            with patch("verify_results.HERE", Path(temporary).resolve()):
                with self.assertRaisesRegex(ValueError, "Missing bundle manifest"):
                    verify()

    def test_narrative_drift_is_detected_before_scientific_validation(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            narrative = root / "README.md"
            narrative.write_text("Research only\n")
            save(root / "manifest.json", {"files_sha256": {"README.md": sha(narrative)}})
            narrative.write_text("Unsupported altered claim\n")
            with patch("verify_results.HERE", root), patch("package_review.HERE", root):
                with self.assertRaisesRegex(ValueError, "Bundle checksum changed: README.md"):
                    verify()

    def test_manifest_cannot_silently_omit_a_published_document(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            (root / "methods.md").write_text("Research limitations\n")
            save(root / "manifest.json", {"files_sha256": {}})
            with patch("verify_results.HERE", root), patch("package_review.HERE", root):
                with self.assertRaisesRegex(ValueError, "complete published bundle"):
                    verify()


if __name__ == "__main__":
    unittest.main()
