import numpy as np
from incremental_slam.config import RunConfig
from incremental_slam.geometry.angles import normalize_angle
from incremental_slam.linearization.pose_gnss import linearize_pose_gnss
from incremental_slam.solvers.cholmod_manager import cholmod_available


def test_normalize_angle():
    x = normalize_angle(np.array([0.0, np.pi, -np.pi]))
    assert x.shape == (3,)


def test_linearize_pose_gnss():
    e, A = linearize_pose_gnss(np.array([1.0, 2.0, 0.5]), np.array([1.5, 2.5]))
    assert e.shape == (2,)
    assert A.shape == (2, 3)


def test_cholmod_available_flag_is_bool():
    assert isinstance(cholmod_available(), bool)


def test_runconfig_defaults_to_auto_backend():
    cfg = RunConfig()
    assert cfg.linear_backend == 'auto'
