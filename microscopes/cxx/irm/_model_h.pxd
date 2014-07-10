from libcpp.vector cimport vector
from libcpp.set cimport set
from libc.stddef cimport size_t

from microscopes._shared_ptr_h cimport shared_ptr
from microscopes.cxx.common._entity_state_h cimport entity_based_state_object
from microscopes.cxx.common._random_fwd_h cimport rng_t
from microscopes.cxx.common._typedefs_h cimport hyperparam_bag_t, suffstats_bag_t
from microscopes.cxx.common.sparse_ndarray._dataview_h cimport dataview
from microscopes.cxx._models_h cimport model as component_model

cdef extern from "microscopes/irm/model.hpp" namespace "microscopes::irm":
    ctypedef vector[const dataview *] dataset_t
    cdef cppclass state:
        cppclass relation_t:
            relation_t()
            relation_t(const vector[size_t] &, const shared_ptr[component_model] &) except +

        state(const vector[size_t] &, const vector[relation_t] &) except +

        hyperparam_bag_t get_domain_hp(size_t) except +
        void set_domain_hp(size_t, const hyperparam_bag_t &) except +

        hyperparam_bag_t get_relation_hp(size_t) except +
        void set_relation_hp(size_t, const hyperparam_bag_t &) except +

        void random_initialize(const dataset_t &, rng_t &) except +
        void initialize(const vector[vector[set[size_t]]] &, const dataset_t &, rng_t &) except +

        float score_assignment() except +
        float score_likelihood(rng_t &) except +

    cdef cppclass bound_state(entity_based_state_object):
        bound_state(const shared_ptr[state] &,
                    size_t,
                    const vector[shared_ptr[dataview]] &) except +
