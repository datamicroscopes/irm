# python imports
from microscopes.models import model_descriptor

cdef class model_definition:
    def __cinit__(self, domains, relations):
        cdef vector[size_t] c_domains
        for d in domains:
            if d < 0:
                raise ValueError("negative domain size")
            c_domains.push_back(d)
        self._domains = list(domains)

        cdef vector[c_relation_definition] c_relations
        cdef vector[size_t] c_rdomains
        for rdomains, rmodel in relations:
            c_rdomains.clear()
            for d in rdomains:
                if d < 0 or d >= len(domains):
                    raise ValueError("invalid domain id")
                c_rdomains.push_back(d)
            if not isinstance(rmodel, model_descriptor):
                raise ValueError(
                    "invalid model given: {}".format(repr(rmodel)))
            c_relations.push_back(
                c_relation_definition(
                    c_rdomains,
                    (<_base>rmodel._c_descriptor).get()))
        self._relations = list(relations)

        self._thisptr.reset(new c_model_definition(c_domains, c_relations))
