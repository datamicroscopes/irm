from microscopes.irm import model, runner
from microscopes.irm.definition import model_definition
from microscopes.models import (
    bb,
    bbnc,
    nich,
    niw,
)
from microscopes.common.rng import rng
from microscopes.common.relation.dataview import numpy_dataview
from microscopes.irm.testutil import toy_dataset


def test_runner_default_kernel_config():
    defn = model_definition([10, 10], [((0, 0), bb), ((0, 1), nich)])
    views = map(numpy_dataview, toy_dataset(defn))
    kc = runner.default_kernel_config(defn)
    prng = rng()
    latent = model.initialize(defn, views, prng)
    r = runner.runner(defn, views, latent, kc, r=prng)
    r.run(10)
