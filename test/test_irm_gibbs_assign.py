from nose.plugins.attrib import attr

from microscopes.irm.definition import model_definition as irm_definition
from microscopes.irm.model import (
    initialize as irm_initialize,
    bind as irm_bind,
)
from microscopes.mixture.definition import model_definition as mm_definition
from microscopes.mixture.model import (
    initialize as mm_initialize,
    bind as mm_bind,
)
from microscopes.common.rng import rng
from microscopes.common.relation.dataview import (
    numpy_dataview as relation_numpy_dataview,
    sparse_2d_dataview as sparse_relation_dataview,
)
from microscopes.common.recarray.dataview import (
    numpy_dataview as rec_numpy_dataview,
)
from microscopes.models import bb, bbnc, gp
from microscopes.kernels.gibbs import assign, assign_resample
from microscopes.kernels.slice import theta

import numpy as np
import numpy.ma as ma

import itertools as it
import operator as op
import time

from scipy.sparse import coo_matrix

from nose.tools import assert_almost_equals
#from nose.plugins.attrib import attr
from microscopes.common.testutil import (
    assert_1d_lists_almost_equals,
    assert_discrete_dist_approx,
    permutation_iter,
    permutation_canonical,
    scores_to_probs,
)


def _tocsr(raw):
    n, m = raw.shape

    def indices():
        for i, j in it.product(range(n), range(m)):
            if not raw.mask[i, j]:
                yield i, j
    data = [raw[i, j] for i, j in indices()]
    i = list(map(op.itemgetter(0), indices()))
    j = list(map(op.itemgetter(1), indices()))
    return coo_matrix((data, (i, j)), shape=raw.shape).tocsr()


def test_compare_to_mixture_model():
    r = rng()

    N, D = 4, 5

    Y = np.random.uniform(size=(N, D)) > 0.8
    Y_rec = np.array([tuple(y) for y in Y], dtype=[('', bool)] * D)

    mm_view = rec_numpy_dataview(Y_rec)
    irm_view = relation_numpy_dataview(Y)

    mm_def = mm_definition(N, [bb] * D)
    irm_def = irm_definition([N, D], [((0, 1), bb)])

    perms = list(permutation_iter(N))
    assignment = perms[np.random.randint(0, len(perms))]

    mm_s = mm_initialize(mm_def, mm_view, r=r, assignment=assignment)
    irm_s = irm_initialize(irm_def,
                           [irm_view],
                           r=r,
                           domain_assignments=[
                               assignment,
                               range(D),
                           ])

    def assert_suff_stats_equal():
        assert set(mm_s.groups()) == set(irm_s.groups(0))
        assert irm_s.groups(1) == range(D)
        groups = mm_s.groups()
        for g in groups:
            for i in xrange(D):
                a = mm_s.get_suffstats(g, i)
                b = irm_s.get_suffstats(0, [g, i])
                if b is None:
                    b = {'heads': 0L, 'tails': 0L}
                assert a['heads'] == b['heads'] and a['tails'] == b['tails']

    assert_suff_stats_equal()
    assert_almost_equals(
        mm_s.score_assignment(), irm_s.score_assignment(0), places=3)

    bound_mm_s = mm_bind(mm_s, mm_view)
    bound_irm_s = irm_bind(irm_s, 0, [irm_view])

    # XXX: doesn't really have to be true, just is true of impl
    assert not bound_mm_s.empty_groups()
    assert not bound_irm_s.empty_groups()

    bound_mm_s.create_group(r)
    bound_irm_s.create_group(r)

    gid_a = bound_mm_s.remove_value(0, r)
    gid_b = bound_irm_s.remove_value(0, r)

    assert gid_a == gid_b
    assert_suff_stats_equal()

    x0, y0 = bound_mm_s.score_value(0, r)
    x1, y1 = bound_irm_s.score_value(0, r)
    assert x0 == x1  # XXX: not really a requirement

    # XXX: should really normalize and then check
    for a, b in zip(y0, y1):
        assert_almost_equals(a, b, places=2)


def test_dense_vs_sparse():
    # XXX: really belongs in irm test cases, but kernels has a nice cluster
    # enumeration iterator

    r = rng()

    n = 5
    raw = ma.array(
        np.random.choice(np.arange(20), size=(n, n)),
        mask=np.random.choice([False, True], size=(n, n)))

    dense = [relation_numpy_dataview(raw)]
    sparse = [sparse_relation_dataview(_tocsr(raw))]

    domains = [n]
    relations = [((0, 0), gp)]
    defn = irm_definition(domains, relations)

    def score_fn(data):
        def f(assignments):
            s = irm_initialize(defn, data, r=r, domain_assignments=assignments)
            assign = sum(s.score_assignment(i)
                         for i in xrange(len(assignments)))
            likelihood = s.score_likelihood(r)
            return assign + likelihood
        return f

    product_assignments = tuple(map(list, map(permutation_iter, domains)))

    dense_posterior = scores_to_probs(
        np.array(map(score_fn(dense), it.product(*product_assignments))))
    sparse_posterior = scores_to_probs(
        np.array(map(score_fn(sparse), it.product(*product_assignments))))

    assert_1d_lists_almost_equals(dense_posterior, sparse_posterior, places=3)


def _test_convergence(domains,
                      data,
                      reg_relations,
                      brute_relations,
                      kernel,
                      burnin_niters=10000,
                      skip=10,
                      ntries=50,
                      nsamples=1000,
                      places=2):
    r = rng()

    reg_defn = irm_definition(domains, reg_relations)
    brute_defn = irm_definition(domains, brute_relations)

    def score_fn(assignments):
        s = irm_initialize(
            brute_defn, data, r=r,
            domain_assignments=assignments)
        assign = sum(s.score_assignment(i) for i in xrange(len(assignments)))
        likelihood = s.score_likelihood(r)
        return assign + likelihood
    product_assignments = tuple(map(list, map(permutation_iter, domains)))
    posterior = scores_to_probs(
        np.array(map(score_fn, it.product(*product_assignments))))

    s = irm_initialize(reg_defn, data, r=r)
    bounded_states = [irm_bind(s, i, data) for i in xrange(len(domains))]

    # burnin
    start = time.time()
    last = start
    for i in xrange(burnin_niters):
        for bs in bounded_states:
            kernel(bs, r)
        if not ((i + 1) % 1000):
            print 'burning finished iteration', (i + 1), \
                'in', (time.time() - last), 'seconds'
            last = time.time()
    print 'finished burnin of', burnin_niters, \
        'iters in', (time.time() - start), 'seconds'

    idmap = {C: i for i, C in enumerate(it.product(*product_assignments))}
    #print idmap

    def sample_fn():
        for _ in xrange(skip):
            for bs in bounded_states:
                kernel(bs, r)
        key = tuple(tuple(permutation_canonical(bs.assignments()))
                    for bs in bounded_states)
        return idmap[key]

    assert_discrete_dist_approx(
        sample_fn, posterior,
        ntries=ntries, nsamples=nsamples,
        kl_places=places)


@attr('slow')
def test_one_binary():
    # 1 domain, 1 binary relation
    domains = [4]

    def mk_relations(model):
        return [((0, 0), model)]

    relsize = (domains[0], domains[0])
    data = [relation_numpy_dataview(
        ma.array(
            np.random.choice([False, True], size=relsize),
            mask=np.random.choice([False, True], size=relsize)))]
    _test_convergence(
        domains, data, mk_relations(bb), mk_relations(bb), assign)


@attr('slow')
def test_one_binary_sparse():
    # 1 domain, 1 binary relation
    domains = [4]

    def mk_relations(model):
        return [((0, 0), model)]

    relsize = (domains[0], domains[0])
    raw = ma.array(
        np.random.choice([False, True], size=relsize),
        mask=np.random.choice([False, True], size=relsize))
    data = [sparse_relation_dataview(_tocsr(raw))]
    _test_convergence(
        domains, data, mk_relations(bb), mk_relations(bb), assign)


@attr('slow')
def test_one_binary_nonconj_kernel():
    # 1 domain, 1 binary relation
    domains = [4]

    def mk_relations(model):
        return [((0, 0), model)]

    relsize = (domains[0], domains[0])
    data = [relation_numpy_dataview(
        ma.array(
            np.random.choice([False, True], size=relsize),
            mask=np.random.choice([False, True], size=relsize)))]
    kernel = lambda s, r: assign_resample(s, 10, r)
    _test_convergence(
        domains, data, mk_relations(bb), mk_relations(bb), kernel)


@attr('slow')
def test_two_binary():
    # 1 domain, 2 binary relations
    domains = [4]

    def mk_relations(model):
        return [((0, 0), model), ((0, 0), model)]

    relsize = (domains[0], domains[0])
    data = [
        relation_numpy_dataview(
            ma.array(
                np.random.choice([False, True], size=relsize),
                mask=np.random.choice([False, True], size=relsize))),
        relation_numpy_dataview(
            ma.array(
                np.random.choice([False, True], size=relsize),
                mask=np.random.choice([False, True], size=relsize))),
    ]
    _test_convergence(
        domains, data, mk_relations(bb), mk_relations(bb), assign)


@attr('slow')
def test_one_binary_one_ternary():
    # 1 domain, 1 binary, 1 ternary
    domains = [4]

    def mk_relations(model):
        return [((0, 0), model), ((0, 0, 0), model)]

    relsize2 = (domains[0], domains[0])
    relsize3 = (domains[0], domains[0], domains[0])
    data = [
        relation_numpy_dataview(
            ma.array(
                np.random.choice([False, True], size=relsize2),
                mask=np.random.choice([False, True], size=relsize2))),
        relation_numpy_dataview(
            ma.array(
                np.random.choice(
                    [False, True], size=relsize3),
                mask=np.random.choice([False, True], size=relsize3))),
    ]
    _test_convergence(
        domains, data, mk_relations(bb), mk_relations(bb), assign)


@attr('slow')
def test_one_binary_nonconj():
    # 1 domain, 1 binary relation, nonconj
    domains = [3]

    def mk_relations(model):
        return [((0, 0), model)]

    relsize = (domains[0], domains[0])
    data = [relation_numpy_dataview(
        ma.array(
            np.random.choice([False, True], size=relsize),
            mask=np.random.random(size=relsize) > 0.8))]

    def mkparam():
        return {'p': 0.05}
    params = {0: mkparam()}

    def kernel(s, r):
        assign_resample(s, 10, r)
        theta(s, r, tparams=params)
    _test_convergence(
        domains, data, mk_relations(bbnc), mk_relations(bb), kernel)


@attr('slow')
def test_two_domain_two_binary():
    # 1 domain, 2 binary relations
    domains = [3, 4]

    def mk_relations(model):
        return [((0, 0), model), ((1, 0), model)]

    relsize00 = (domains[0], domains[0])
    relsize10 = (domains[1], domains[0])
    data = [
        relation_numpy_dataview(
            ma.array(
                np.random.choice([False, True], size=relsize00),
                mask=np.random.choice([False, True], size=relsize00))),
        relation_numpy_dataview(
            ma.array(
                np.random.choice([False, True], size=relsize10),
                mask=np.random.choice([False, True], size=relsize10))),
    ]
    _test_convergence(
        domains, data, mk_relations(bb), mk_relations(bb), assign)
