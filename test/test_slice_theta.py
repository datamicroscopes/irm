from microscopes.kernels.slice import theta
from microscopes.irm.model import initialize, bind
from microscopes.common.relation.dataview import numpy_dataview
from microscopes.common.rng import rng
from microscopes.models import bbnc
from microscopes.irm.definition import model_definition

from microscopes.common.testutil import assert_1d_cont_dist_approx_sps
from scipy.stats import beta

import numpy as np

#from nose.plugins.attrib import attr


def test_slice_theta_irm():
    N = 10
    defn = model_definition([N], [((0, 0), bbnc)])
    data = np.random.random(size=(N, N)) < 0.8
    view = numpy_dataview(data)
    r = rng()
    prior = {'alpha': 1.0, 'beta': 9.0}

    s = initialize(
        defn,
        [view],
        r=r,
        cluster_hps=[{'alpha': 2.0}],
        relation_hps=[prior],
        domain_assignments=[[0] * N])

    bs = bind(s, 0, [view])

    params = {0: {'p': 0.05}}

    heads = len([1 for y in data.flatten() if y])
    tails = len([1 for y in data.flatten() if not y])

    alpha1 = prior['alpha'] + heads
    beta1 = prior['beta'] + tails

    def sample_fn():
        theta(bs, r, tparams=params)
        return s.get_suffstats(0, [0, 0])['p']

    rv = beta(alpha1, beta1)
    assert_1d_cont_dist_approx_sps(sample_fn, rv, nsamples=50000)
