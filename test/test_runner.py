from microscopes.irm import model, runner
from microscopes.irm.definition import model_definition
from microscopes.models import (
    bb,
    bbnc,
    nich,
)
from microscopes.common.rng import rng
from microscopes.common.relation.dataview import numpy_dataview
from microscopes.kernels import parallel
from microscopes.irm.testutil import toy_dataset, data_with_posterior
from microscopes.common.testutil import (
    assert_discrete_dist_approx,
    permutation_iter,
    permutation_canonical,
)

import itertools as it
import multiprocessing as mp

from nose.plugins.attrib import attr


def _test_runner_simple(defn, kc_fn):
    views = map(numpy_dataview, toy_dataset(defn))
    kc = kc_fn(defn)
    prng = rng()
    latent = model.initialize(defn, views, prng)
    r = runner.runner(defn, views, latent, kc)
    r.run(prng, 10)


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


def test_runner_default_kernel_config_convergence():
    domains = [4]
    defn = model_definition(domains, [((0, 0), bb)])
    prng = rng()
    relations, posterior = data_with_posterior(defn, prng)
    views = map(numpy_dataview, relations)
    latent = model.initialize(defn, views, prng)
    r = runner.runner(defn, views, latent, [('assign', range(len(domains)))])

    r.run(r=prng, niters=1000)  # burnin
    product_assignments = tuple(map(list, map(permutation_iter, domains)))
    idmap = {C: i for i, C in enumerate(it.product(*product_assignments))}

    def sample_fn():
        r.run(r=prng, niters=10)
        new_latent = r.get_latent()
        key = tuple(tuple(permutation_canonical(new_latent.assignments(i)))
                    for i in xrange(len(domains)))
        return idmap[key]

    assert_discrete_dist_approx(sample_fn, posterior, ntries=100)


@attr('uses_mp')
def test_runner_multiprocessing():
    defn = model_definition([10, 10], [((0, 0), bb), ((0, 1), nich)])
    views = map(numpy_dataview, toy_dataset(defn))
    kc = runner.default_kernel_config(defn)
    prng = rng()
    latents = [model.initialize(defn, views, prng)
               for _ in xrange(mp.cpu_count())]
    runners = [runner.runner(defn, views, latent, kc) for latent in latents]
    r = parallel.runner(runners)
    # check it is restartable
    r.run(r=prng, niters=10)
    r.run(r=prng, niters=10)


@attr('uses_mp')
def test_runner_multiprocessing_convergence():
    domains = [4]
    defn = model_definition(domains, [((0, 0), bb)])
    prng = rng()
    relations, posterior = data_with_posterior(defn, prng)
    views = map(numpy_dataview, relations)
    latents = [model.initialize(defn, views, prng)
               for _ in xrange(mp.cpu_count())]
    kc = [('assign', range(len(domains)))]
    runners = [runner.runner(defn, views, latent, kc) for latent in latents]
    r = parallel.runner(runners)
    r.run(r=prng, niters=10000)  # burnin
    product_assignments = tuple(map(list, map(permutation_iter, domains)))
    idmap = {C: i for i, C in enumerate(it.product(*product_assignments))}

    def sample_iter():
        r.run(r=prng, niters=10)
        for latent in r.get_latents():
            key = tuple(tuple(permutation_canonical(latent.assignments(i)))
                        for i in xrange(len(domains)))
            yield idmap[key]

    ref = [None]

    def sample_fn():
        if ref[0] is None:
            ref[0] = sample_iter()
        try:
            return next(ref[0])
        except StopIteration:
            ref[0] = None
        return sample_fn()

    assert_discrete_dist_approx(sample_fn, posterior, ntries=100, kl_places=2)
