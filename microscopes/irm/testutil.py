"""Test helpers specific to IRM

"""

import numpy as np
import itertools as it

from microscopes.common.rng import rng
from microscopes.common.relation.dataview import numpy_dataview
from microscopes.common.testutil import (
    permutation_iter,
    scores_to_probs,
)
from microscopes.irm import model


def toy_dataset(defn):

    def make(shape, model):
        # XXX(stephentu): create public getter for model_module
        module = model.py_desc()._model_module
        shared = module.Shared()
        shared.load(model.default_hyperparams())
        sampler = module.Sampler()
        sampler.init(shared)
        samples = [sampler.eval(shared) for _ in xrange(np.product(shape))]
        return np.array(samples).reshape(shape)

    return [make(defn.shape(i), model)
            for i, model in enumerate(defn.relation_models())]


def data_with_posterior(defn, r=None):
    # XXX(stephentu): should only accept conjugate models
    if r is None:
        r = rng()
    relations = toy_dataset(defn)
    views = map(numpy_dataview, relations)

    def score_fn(assignments):
        s = model.initialize(defn, views, r=r, domain_assignments=assignments)
        assign = sum(s.score_assignment(i) for i in xrange(len(assignments)))
        likelihood = s.score_likelihood(r)
        return assign + likelihood

    domains = defn.domains()
    product_assignments = tuple(map(list, map(permutation_iter, domains)))
    posterior = scores_to_probs(
        np.array(map(score_fn, it.product(*product_assignments))))

    return relations, posterior
