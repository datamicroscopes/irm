# cython imports
from libcpp.vector cimport vector
from libc.stddef cimport size_t

from microscopes._shared_ptr_h cimport shared_ptr
from microscopes.cxx._models cimport _base
from microscopes.cxx.irm._model_h cimport \
    model_definition as c_model_definition, \
    relation_definition as c_relation_definition
 
cdef class model_definition:
    cdef shared_ptr[c_model_definition] _thisptr
    cdef public list _domains
    cdef public list _relations
