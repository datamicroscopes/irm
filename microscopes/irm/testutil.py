"""Test helpers specific to mixture models

"""

import numpy as np


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
