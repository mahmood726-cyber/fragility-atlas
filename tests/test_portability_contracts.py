import json
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = REPO_ROOT / "e156-submission" / "config.json"
VALIDATION_PATHS = [
    REPO_ROOT / "tests" / "export_for_r_validation.py",
    REPO_ROOT / "tests" / "validate_vs_metafor.R",
    REPO_ROOT / "tests" / "validate_proper.R",
    REPO_ROOT / "tests" / "validate_all_estimators.R",
    REPO_ROOT / "tests" / "check_cd001431.R",
]


class PortabilityContracts(unittest.TestCase):
    def test_submission_config_uses_repo_relative_root(self) -> None:
        payload = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))

        self.assertEqual(payload["path"], "..")
        resolved_root = (CONFIG_PATH.parent / payload["path"]).resolve()
        self.assertEqual(resolved_root, REPO_ROOT.resolve())

    def test_helper_scripts_use_repo_relative_defaults(self) -> None:
        populate = (REPO_ROOT / "src" / "populate_manuscript.py").read_text(encoding="utf-8")
        embed = (REPO_ROOT / "src" / "embed_dashboard_data.py").read_text(encoding="utf-8")

        self.assertIn("DEFAULT_OUTPUT_DIR = REPO_ROOT / 'data' / 'output'", populate)
        self.assertIn("DEFAULT_MANUSCRIPT_PATH = REPO_ROOT / 'manuscript_bmj.md'", populate)
        self.assertNotIn(r"C:\\FragilityAtlas", populate)
        self.assertIn("DEFAULT_OUTPUT_DIR = REPO_ROOT / 'data' / 'output'", embed)
        self.assertIn("DEFAULT_DASHBOARD_PATH = REPO_ROOT / 'dashboard' / 'index.html'", embed)
        self.assertNotIn(r"C:\\FragilityAtlas", embed)

    def test_validation_scripts_are_repo_local_and_env_driven(self) -> None:
        for path in VALIDATION_PATHS:
            text = path.read_text(encoding="utf-8")
            self.assertNotIn(r"C:\FragilityAtlas", text, path.as_posix())
            self.assertNotIn("C:/FragilityAtlas", text, path.as_posix())
            self.assertNotIn(r"C:\Models\Pairwise70\data", text, path.as_posix())
            self.assertNotIn("C:/Models/Pairwise70/data", text, path.as_posix())


if __name__ == "__main__":
    unittest.main()
