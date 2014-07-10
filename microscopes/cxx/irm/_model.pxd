from libcpp.vector cimport vector
from libc.stddef cimport size_t

from microscopes._shared_ptr_h cimport shared_ptr
from microscopes.cxx._models cimport factory
from microscopes.cxx._models_h cimport model as component_model
from microscopes.cxx.common.sparse_ndarray._dataview cimport abstract_dataview
from microscopes.cxx.irm._model_h cimport state as c_state, bound_state as c_bound_state

cdef class state:
    cdef shared_ptr[c_state] _thisptr
    cdef list _domains
    cdef list _relations
