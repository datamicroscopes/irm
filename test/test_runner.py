from microscopes.irm import model, runner
from microscopes.irm.definition import model_definition
from microscopes.models import (
    bb,
    bbnc,
    nich,
)
from microscopes.common.rng import rng
from microscopes.common.relation.dataview import numpy_dataview
from microscopes.irm.testutil import toy_dataset

import itertools as it


def _test_runner_simple(defn, kc_fn):
    views = map(numpy_dataview, toy_dataset(defn))
    kc = kc_fn(defn)
    prng = rng()
    latent = model.initialize(defn, views, prng)
    r = runner.runner(defn, views, latent, kc, r=prng)
    r.run(10)


def test_runner_default_kernel_config():
    defn = model_definition([10, 10], [((0, 0), bb), ((0, 1), nich)])
    kc_fn = runner.default_kernel_config
    _test_runner_simple(defn, kc_fn)


def test_runner_default_kernel_config_nonconj():
    defn = model_definition([10, 10], [((0, 0), bbnc), ((0, 1), nich)])
    kc_fn = runner.default_kernel_config
    _test_runner_simple(defn, kc_fn)


def test_runner_default_kernel_config_with_cluster():
    defn = model_definition([10, 10], [((0, 0), bb), ((0, 1), nich)])

    def kc_fn(defn):
        return list(it.chain(
            runner.default_assign_kernel_config(defn),
            runner.default_relation_hp_kernel_config(defn),
            runner.default_cluster_hp_kernel_config(defn)))
    _test_runner_simple(defn, kc_fn)
