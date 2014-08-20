"""Implements the Runner interface for IRM
"""

from microscopes.common import validator
from microscopes.common.rng import rng
from microscopes.common.relation._dataview import abstract_dataview
from microscopes.irm.definition import model_definition
from microscopes.irm.model import state, bind
from microscopes.kernels import gibbs, slice

import itertools as it


def default_assign_kernel_config(defn):
    """Creates a default kernel configuration for sampling the assignment
    (clustering) vector for every domain. The default kernel is currently a
    gibbs sampler.

    Parameters
    ----------
    defn : irm definition

    """
    validator.validate_type(defn, model_definition, 'defn')

    # XXX(stephentu): model_descriptors should implement
    # is_conjugate()

    def is_nonconj(x):
        return x.name() == 'bbnc'

    conj_inds, nonconj_inds = [], []
    for idx, m in enumerate(defn.relation_models()):
        lst = nonconj_inds if is_nonconj(m) else conj_inds
        lst.append(idx)

    assign_kernel = ('assign', conj_inds)
    assign_resample_kernel = (
        'assign_resample',
        {idx: {'m': 10} for idx in nonconj_inds}
    )
    theta_kernel = (
        'theta',
        {'tparams': {idx: {'p': 0.1} for idx in nonconj_inds}}
    )

    kernels = []
    if conj_inds:
        kernels.append(assign_kernel)
    if nonconj_inds:
        kernels.append(assign_resample_kernel)
        kernels.append(theta_kernel)
    return kernels


def default_relation_hp_kernel_config(defn):
    """Creates a default kernel configuration for sampling the component
    (feature) model hyper-parameters. The default kernel is currently
    a one-dimensional slice sampler.

    Parameters
    ----------
    defn : irm definition
        The hyper-priors set in the definition are used to configure the
        hyper-parameter sampling kernels.

    """
    validator.validate_type(defn, model_definition, 'defn')
    hparams = {}
    for i, hp in enumerate(defn.relation_hyperpriors()):
        if not hp:
            continue
        # XXX(stephentu): we are arbitrarily picking w=0.1
        hparams[i] = {k: (fn, 0.1) for k, fn in hp.iteritems()}
    return [('relation_hp', {'hparams': hparams})]


def default_cluster_hp_kernel_config(defn):
    """Creates a default kernel configuration for sampling the clustering
    (Chinese Restaurant Process) model hyper-parameter. The default kernel is
    currently a one-dimensional slice sampler.

    Parameters
    ----------
    defn : irm definition
        The hyper-priors set in the definition are used to configure the
        hyper-parameter sampling kernels.
    """
    validator.validate_type(defn, model_definition, 'defn')
    config = {}
    for i, hp in enumerate(defn.domain_hyperpriors()):
        if not hp:
            continue
        cparam = {k: (fn, 0.1) for k, fn in hp.iteritems()}
        config[i] = {'cparam': cparam}
    return [('cluster_hp', config)]


def default_kernel_config(defn):
    """Creates a default kernel configuration suitable for general purpose
    inference. Currently configures an assignment sampler followed by a
    component hyper-parameter sampler.

    Parameters
    ----------
    defn : irm definition

    """
    # XXX(stephentu): should the default config also include cluster_hp?
    return list(it.chain(
        default_assign_kernel_config(defn),
        default_relation_hp_kernel_config(defn)))


class runner(object):
    # XXX(stephentu): do a better job of documentating the kernel configuration
    # dicts
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
        A list of `(x, y)` tuples where `x` is a string containing the name of
        the kernel and `y` is kernel specific configuration.
    """

    def __init__(self, defn, views, latent, kernel_config):
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
            name, config = kernel

            if not hasattr(config, 'iteritems'):
                config = {c: {} for c in config}
            validator.validate_dict_like(config)

            def require_relation_keys(config):
                valid_keys = set(xrange(len(defn.relations())))
                if not set(config.keys()).issubset(valid_keys):
                    raise ValueError("bad config found: {}".format(config))

            def require_domain_keys(config):
                valid_keys = set(xrange(len(defn.domains())))
                if not set(config.keys()).issubset(valid_keys):
                    raise ValueError("bad config found: {}".format(config))

            if name == 'assign':
                require_domain_keys(config)
                for v in config.values():
                    validator.validate_dict_like(v)
                    if v:
                        msg = "assign has no config params: {}".format(v)
                        raise ValueError(msg)

            elif name == 'assign_resample':
                require_domain_keys(config)
                for v in config.values():
                    validator.validate_dict_like(v)
                    if v.keys() != ['m']:
                        raise ValueError("bad config found: {}".format(v))

            elif name == 'cluster_hp':
                require_domain_keys(config)
                for v in config.values():
                    validator.validate_dict_like(v)
                    if v.keys() != ['cparam']:
                        raise ValueError("bad config found: {}".format(v))

            elif name == 'relation_hp':
                if config.keys() != ['hparams']:
                    raise ValueError("bad config found: {}".format(config))
                validator.validate_dict_like(config['hparams'])
                require_relation_keys(config['hparams'])

            elif name == 'theta':
                if config.keys() != ['tparams']:
                    raise ValueError("bad config found: {}".format(config))
                validator.validate_dict_like(config['tparams'])
                require_relation_keys(config['tparams'])

            else:
                raise ValueError("bad kernel found: {}".format(name))

            self._kernel_config.append((name, config))

    def run(self, r, niters=10000):
        """Run the specified mixturemodel kernel for `niters`, in a single
        thread.

        Parameters
        ----------
        r : random state
        niters : int

        """
        validator.validate_type(r, rng, param_name='r')
        validator.validate_positive(niters, param_name='niters')
        inds = xrange(len(self._defn.domains()))
        models = [bind(self._latent, i, self._views) for i in inds]
        for _ in xrange(niters):
            for name, config in self._kernel_config:
                if name == 'assign':
                    for idx in config.keys():
                        gibbs.assign(models[idx], r)
                elif name == 'assign_resample':
                    for idx, v in config.iteritems():
                        gibbs.assign_resample(models[idx], v['m'], r)
                elif name == 'cluster_hp':
                    for idx, v in config.iteritems():
                        slice.hp(models[idx], r, cparam=v['cparam'])
                elif name == 'relation_hp':
                    slice.hp(models[0], r, hparams=config['hparams'])
                elif name == 'theta':
                    slice.theta(models[0], r, tparams=config['tparams'])
                else:
                    assert False, "should not be reached"

    def get_latent(self):
        return self._latent
