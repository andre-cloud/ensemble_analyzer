import pytest
from unittest.mock import patch, MagicMock
from ensemble_analyzer.cli.check_setup import (
    check_python_dependencies, check_orca, check_gaussian,
    check_nwchem, check_tblite, check_mlip, check_models_dir, main
)


class TestCheckSetup:

    @patch("ensemble_analyzer.cli.check_setup.importlib.util.find_spec")
    def test_python_deps_all_pass(self, mock_find_spec):
        mock_find_spec.return_value = MagicMock()
        assert check_python_dependencies() is True

    @patch("ensemble_analyzer.cli.check_setup.importlib.util.find_spec")
    def test_python_deps_missing(self, mock_find_spec):
        mock_find_spec.side_effect = [None] + [MagicMock()] * 8 + [MagicMock()]
        assert check_python_dependencies() is False

    @patch("ensemble_analyzer.cli.check_setup.importlib.util.find_spec")
    def test_python_deps_missing_package(self, mock_find_spec):
        mock_find_spec.side_effect = [MagicMock()] * 9 + [None]
        assert check_python_dependencies() is False

    @patch("ensemble_analyzer.cli.check_setup.shutil.which")
    def test_orca_found(self, mock_which):
        mock_which.return_value = "/usr/bin/orca"
        with patch.dict("os.environ", {"ORCAVERSION": "5.0.3"}):
            assert check_orca() is True

    @patch("ensemble_analyzer.cli.check_setup.shutil.which")
    def test_orca_missing_exe(self, mock_which):
        mock_which.return_value = None
        assert check_orca() is False

    @patch("ensemble_analyzer.cli.check_setup.shutil.which")
    def test_orca_missing_env(self, mock_which):
        mock_which.return_value = "/usr/bin/orca"
        with patch.dict("os.environ", {}, clear=True):
            assert check_orca() is False

    @patch("ensemble_analyzer.cli.check_setup.shutil.which")
    def test_gaussian_found_g16(self, mock_which):
        mock_which.side_effect = lambda x: "/usr/bin/g16" if x == "g16" else None
        check_gaussian()

    @patch("ensemble_analyzer.cli.check_setup.shutil.which")
    def test_gaussian_found_g09(self, mock_which):
        mock_which.side_effect = lambda x: "/usr/bin/g09" if x == "g09" else None
        check_gaussian()

    @patch("ensemble_analyzer.cli.check_setup.shutil.which")
    def test_gaussian_not_found(self, mock_which):
        mock_which.return_value = None
        check_gaussian()

    @patch("ensemble_analyzer.cli.check_setup.shutil.which")
    def test_nwchem_found(self, mock_which):
        mock_which.side_effect = lambda x: "/usr/bin/nwchem" if x == "nwchem" else None
        check_nwchem()

    @patch("ensemble_analyzer.cli.check_setup.shutil.which")
    def test_nwchem_not_found(self, mock_which):
        mock_which.return_value = None
        check_nwchem()

    @patch("ensemble_analyzer.cli.check_setup.importlib.util.find_spec")
    def test_tblite_found(self, mock_find_spec):
        mock_find_spec.return_value = MagicMock()
        check_tblite()

    @patch("ensemble_analyzer.cli.check_setup.importlib.util.find_spec")
    def test_tblite_not_found(self, mock_find_spec):
        mock_find_spec.return_value = None
        check_tblite()

    @patch("ensemble_analyzer.cli.check_setup.importlib.util.find_spec")
    def test_ml_all_missing(self, mock_find_spec):
        mock_find_spec.return_value = None
        check_mlip()

    @patch("ensemble_analyzer.cli.check_setup.importlib.util.find_spec")
    def test_ml_aimnet_torch(self, mock_find_spec):
        def side_effect(name):
            found = {"aimnet": MagicMock(), "torch": MagicMock()}
            return found.get(name)
        mock_find_spec.side_effect = side_effect
        check_mlip()

    @patch("ensemble_analyzer.cli.check_setup.importlib.util.find_spec")
    def test_ml_torch_without_aimnet(self, mock_find_spec):
        def side_effect(name):
            found = {"torch": MagicMock()}
            return found.get(name)
        mock_find_spec.side_effect = side_effect
        check_mlip()

    @patch("ensemble_analyzer.cli.check_setup.importlib.util.find_spec")
    def test_ml_aimnet_without_torch(self, mock_find_spec):
        def side_effect(name):
            found = {"aimnet": MagicMock()}
            return found.get(name)
        mock_find_spec.side_effect = side_effect
        check_mlip()

    @patch("ensemble_analyzer.cli.check_setup.importlib.util.find_spec")
    def test_ml_uma_found(self, mock_find_spec):
        def side_effect(name):
            found = {"fairchem": MagicMock()}
            return found.get(name)
        mock_find_spec.side_effect = side_effect
        check_mlip()

    @patch("ensemble_analyzer.cli.check_setup.importlib.util.find_spec")
    def test_ml_mace_found(self, mock_find_spec):
        def side_effect(name):
            found = {"mace": MagicMock()}
            return found.get(name)
        mock_find_spec.side_effect = side_effect
        check_mlip()

    @patch("ensemble_analyzer.cli.check_setup.importlib.util.find_spec")
    def test_ml_all_found(self, mock_find_spec):
        mock_find_spec.return_value = MagicMock()
        check_mlip()

    def test_models_dir_set_and_exists(self):
        with patch("ensemble_analyzer.cli.check_setup.os.environ.get", return_value="/tmp/models"):
            with patch("ensemble_analyzer.cli.check_setup.Path.is_dir", return_value=True):
                check_models_dir()

    def test_models_dir_set_not_exists(self):
        with patch("ensemble_analyzer.cli.check_setup.os.environ.get", return_value="/tmp/models"):
            with patch("ensemble_analyzer.cli.check_setup.Path.is_dir", return_value=False):
                check_models_dir()

    def test_models_dir_not_set(self):
        with patch("ensemble_analyzer.cli.check_setup.os.environ.get", return_value=None):
            check_models_dir()

    @patch("ensemble_analyzer.cli.check_setup.check_python_dependencies")
    @patch("ensemble_analyzer.cli.check_setup.check_orca")
    @patch("ensemble_analyzer.cli.check_setup.check_gaussian")
    @patch("ensemble_analyzer.cli.check_setup.check_nwchem")
    @patch("ensemble_analyzer.cli.check_setup.check_tblite")
    @patch("ensemble_analyzer.cli.check_setup.check_mlip")
    @patch("ensemble_analyzer.cli.check_setup.check_models_dir")
    def test_main_success(self, mock_models, mock_ml, mock_tblite, mock_nwchem, mock_gauss, mock_orca, mock_deps):
        mock_deps.return_value = True
        mock_orca.return_value = True
        with pytest.raises(SystemExit) as exc:
            main()
        assert exc.value.code == 0

    @patch("ensemble_analyzer.cli.check_setup.check_python_dependencies")
    @patch("ensemble_analyzer.cli.check_setup.check_orca")
    @patch("ensemble_analyzer.cli.check_setup.check_gaussian")
    @patch("ensemble_analyzer.cli.check_setup.check_nwchem")
    @patch("ensemble_analyzer.cli.check_setup.check_tblite")
    @patch("ensemble_analyzer.cli.check_setup.check_mlip")
    @patch("ensemble_analyzer.cli.check_setup.check_models_dir")
    def test_main_failure(self, mock_models, mock_ml, mock_tblite, mock_nwchem, mock_gauss, mock_orca, mock_deps):
        mock_deps.return_value = False
        mock_orca.return_value = False
        with pytest.raises(SystemExit) as exc:
            main()
        assert exc.value.code == 1
