import warnings

import pytest

from enan_calculators._model_tasks import MODEL_TASKS, check_model_task, tasks_for_model


def _no_warning(fn, *args):
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        fn(*args)


class TestTasksForModel:
    def test_uma_1p2_has_oc22_oc25(self):
        assert {"oc22", "oc25"} <= set(tasks_for_model("fairchem", "uma-s-1p2.pt"))
        assert "oc22" not in tasks_for_model("fairchem", "uma-s-1p1.pt")

    def test_single_task_finetunes(self):
        assert tasks_for_model("fairchem", "uma-s-1p2p1-ft-r2scan.pt") == ["omat"]
        assert tasks_for_model("fairchem", "uma-s-1p2p1-ft-d4rev.pt") == ["omol"]

    def test_uma_alias_and_path(self):
        assert tasks_for_model("uma", "/models/uma/uma-s-1p1.pt") == tasks_for_model("fairchem", "uma-s-1p1.pt")

    def test_mace_cache_prefixed_name(self):
        assert tasks_for_model("mace", "mace-mp-mh-1-mace-mh-1.model") == tasks_for_model("mace", "mace-mh-1.model")

    def test_unknown_model(self):
        assert tasks_for_model("mace", "my_finetune.model") is None

    def test_uma_is_fairchem(self):
        assert MODEL_TASKS["uma"] is MODEL_TASKS["fairchem"]


class TestCheckModelTask:
    def test_supported_task_passes(self):
        _no_warning(check_model_task, "uma", "uma-s-1p2.pt", "omol")
        _no_warning(check_model_task, "mace", "mace-mh-1.model", "spice_wB97M")

    def test_unsupported_task_raises(self):
        with pytest.raises(ValueError, match="Supported: omat"):
            check_model_task("fairchem", "uma-s-1p2p1-ft-r2scan.pt", "omol")
        with pytest.raises(ValueError, match="does not support task 'rgd1_b3lyp'"):
            check_model_task("mace", "mace-mh-1.model", "rgd1_b3lyp")

    def test_mace_head_case_insensitive(self):
        _no_warning(check_model_task, "mace", "MACE-OFF23_small.model", "default")
        _no_warning(check_model_task, "mace", "EGRET_1.model", "Default")

    def test_default_task_left_to_backend(self):
        _no_warning(check_model_task, "mace", "mace-mh-0.model", None)

    def test_unknown_model_warns_no_sanity_check(self):
        with pytest.warns(UserWarning, match="no sanity check"):
            check_model_task("fairchem", "my_finetune.pt", "omat")
