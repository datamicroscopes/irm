from libcpp.vector cimport vector
from libc.stddef cimport size_t

from microscopes._shared_ptr_h cimport shared_ptr
from microscopes.cxx.common._entity_state_h cimport entity_based_state_object
from microscopes.cxx.common.sparse_ndarray._dataview_h cimport dataview
from microscopes.cxx._models_h cimport model as component_model

cdef extern from "microscopes/irm/model.hpp" namespace "microscopes::irm":
    cdef cppclass state:
        cppclass relation_t:
            relation_t()
            relation_t(const vector[size_t] &, const shared_ptr[component_model] &) except +
        state(const vector[size_t] &, const vector[relation_t] &) except +

    cdef cppclass bound_state(entity_based_state_object):
        bound_state(const shared_ptr[state] &,
                    size_t,
                    const vector[shared_ptr[dataview]] &) except +
