from libcpp.vector cimport vector
from libcpp.set cimport set
from libcpp.string cimport string
from libc.stddef cimport size_t
from libcpp cimport bool

from microscopes._shared_ptr_h cimport shared_ptr
from microscopes.common._entity_state_h cimport entity_based_state_object
from microscopes.common._random_fwd_h cimport rng_t
from microscopes.common._typedefs_h cimport hyperparam_bag_t, suffstats_bag_t
from microscopes.common.relation._dataview_h cimport dataview
from microscopes._models_h cimport model as c_model

cdef extern from "microscopes/irm/model.hpp" namespace "microscopes::irm":
    ctypedef vector[const dataview *] dataset_t

    cdef cppclass relation_definition:
        relation_definition()
        relation_definition(const vector[size_t] &,
                            const shared_ptr[c_model] &) except +

    cdef cppclass model_definition:
        model_definition(const vector[size_t] &,
                         const vector[relation_definition] &) except +

    cdef cppclass state:

        size_t ndomains()
        size_t nrelations()
        size_t nentities(size_t) except +
        size_t ngroups(size_t) except +
        const vector[ssize_t] & assignments(size_t) except +
        vector[size_t] groups(size_t) except +
        const set[size_t] & empty_groups(size_t) except +
        bool isactivegroup(size_t, size_t) except +
        size_t groupsize(size_t, size_t) except +

        hyperparam_bag_t get_domain_hp(size_t) except +
        void set_domain_hp(size_t, const hyperparam_bag_t &) except +

        hyperparam_bag_t get_relation_hp(size_t) except +
        void set_relation_hp(size_t, const hyperparam_bag_t &) except +

        bool get_suffstats(size_t, const vector[size_t] &, suffstats_bag_t &) except +
        # XXX(stephentu): no set_suffstats()

        # XXX(stephentu):
        #size_t create_group(size_t) except +
        #void delete_group(size_t, size_t) except +

        # XXX(stephentu):
        #  add_value()
        #  remove_value()
        #  score_value()

        float score_assignment(size_t) except +
        float score_likelihood(rng_t &) except +

        # stupid testing functions
        vector[vector[size_t]] entity_data_positions(size_t, size_t, const dataset_t &) except +

        string serialize() except +

    cdef cppclass model(entity_based_state_object):
        model(const shared_ptr[state] &,
              size_t,
              const vector[shared_ptr[dataview]] &) except +

cdef extern from "microscopes/irm/model.hpp" namespace "microscopes::irm::state":
    shared_ptr[state] \
    initialize(const model_definition &,
               const vector[hyperparam_bag_t] &,
               const vector[hyperparam_bag_t] &,
               const vector[vector[size_t]] &,
               const dataset_t &,
               rng_t &) except +

    shared_ptr[state] \
    deserialize(const model_definition &, const string &) except +
