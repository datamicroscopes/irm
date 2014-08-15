"""Implements the Runner interface for IRM
"""

from microscopes.common import validator
from microscopes.common.rng import rng
from microscopes.common.relation._dataview import abstract_dataview
from microscopes.irm.definition import model_definition
from microscopes.irm.model import state, bind
from microscopes.kernels import gibbs, slice

def default_kernel_config(defn):
    validator.validate_type(defn, model_definition, 'defn')
    models = defn.relation_models()
    nonconj_models = filter(lambda x: x.name() == 'bbnc', models)

    # assignment
    if nonconj_models:
        # XXX(stephentu): limitation- one nonconj model triggers
        # all relations to be sampled with assign_resample(),
        # which is correct, but slower
        kernels = [('assign_resample', {'m': 10})]
        # XXX(stephentu): also slice sample on the parameter instantiations
    else:
        kernels = ['assign']

    # hyperparams
    hparams = {}
    for i, hp in enumerate(defn.relation_hyperpriors()):
        if not hp:
            continue
        hparams[i] = {k : (fn, 0.1) for k, fn in hp.iteritems()}

    kernels.append(('relation_hp', {'hparams': hparams}))
    return kernels

class runner(object):
    """The IRM model runner

    Parameters
    ----------

    defn : ``model_definition``
        The structural definition.

    views : list
        A list of the relation dataviews.

    latent : ``state``
        The initialization state.

    kernel_config : list
        A list of either `x` strings or `(x, y)` tuples, where `x` is a string
        containing the name of the kernel and `y` is a dict which configures
        the particular kernel. In the former case where `y` is omitted, then
        the defaults parameters for each kernel are used.

    r : ``rng``, optional

    """

    def __init__(self, defn, views, latent, kernel_config, r=None):
        validator.validate_type(defn, model_definition, 'defn')
        validator.validate_len(views, len(defn.relations()), 'views')
        for view in views:
            validator.validate_type(view, abstract_dataview)
        validator.validate_type(latent, state, 'latent')

        self._defn = defn
        self._views = views
        self._latent = latent

        self._kernel_config = []
        for kernel in kernel_config:
            if hasattr(kernel, '__iter__'):
                name, config = kernel
            else:
                name, config = kernel, {}
            if name == 'assign_fixed':
                raise ValueError("state should not use fixed kernel")
            elif name == 'assign':
                if config:
                    raise ValueError("assign has no parameters")
            elif name == 'assign_resample':
                # XXX(stephentu): all domains share the same m for now
                if config.keys() != ['m']:
                    raise ValueError("bad config found: {}".format(config))
            elif name == 'cluster_hp':
                valid_keys = set(xrange(len(defn.domains())))
                if not set(config.keys()).issubset(valid_keys):
                    raise ValueError("bad config found: {}".format(config))
                for args in config.values():
                    if args.keys() != ['cparam']:
                        raise ValueError("bad config found: {}".format(config))
            elif name == 'relation_hp':
                if config.keys() != ['hparams']:
                    raise ValueError("bad config found: {}".format(config))
                validator.validate_type(config['hparams'], dict)
                valid_keys = set(xrange(len(defn.relations())))
                if not set(config['hparams'].keys()).issubset(valid_keys):
                    raise ValueError("bad config found: {}".format(config))
            else:
                raise ValueError("bad kernel found: {}".format(name))
            self._kernel_config.append((name, config))

        if r is None:
            r = rng()
        validator.validate_type(r, rng, 'r')
        self._r = r

    def run(self, niters=10000):
        validator.validate_positive(niters, param_name='niters')
        inds = xrange(len(self._defn.domains()))
        models = [bind(self._latent, i, self._views) for i in inds]
        for _ in xrange(niters):
            for name, config in self._kernel_config:
                if name == 'assign':
                    for model in models:
                        gibbs.assign(model, self._r)
                elif name == 'assign_resample':
                    for model in models:
                        gibbs.assign_resample(model, config['m'], self._r)
                elif name == 'cluster_hp':
                    for idx, args in config.iteritems():
                        slice.hp(models[idx],
                                 self._r,
                                 cparam=config['cparam'])
                elif name == 'relation_hp':
                    slice.hp(models[0],
                             self._r,
                             hparams=config['hparams'])
                else:
                    assert False, "should not be reach"

    def get_latent(self):
        return self._latent
