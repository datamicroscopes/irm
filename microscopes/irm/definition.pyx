# cython: embedsignature=True


# python imports
from microscopes.models import model_descriptor
from microscopes.common import validator


cdef class model_definition:
    def __cinit__(self, domains, relations):
        validator.validate_nonempty(domains, "domains")
        validator.validate_nonempty(relations, "relations")

        cdef vector[size_t] c_domains
        for d in domains:
            validator.validate_positive(d)
            c_domains.push_back(d)
        self._domains = list(domains)

        cdef vector[c_relation_definition] c_relations
        cdef vector[size_t] c_rdomains
        for rdomains, rmodel in relations:
            if len(rdomains) < 2:
                raise ValueError("cannot have relation with arity < 2")
            if len(rdomains) > 4:
                # XXX(stephentu) need to fix in future
                raise ValueError(
                    "implementation limitation: relation arity must be <= 4")
            c_rdomains.clear()
            for d in rdomains:
                validator.validate_in_range(d, len(domains))
                c_rdomains.push_back(d)
            validator.validate_type(rmodel, model_descriptor)
            c_relations.push_back(
                c_relation_definition(
                    c_rdomains,
                    (<_base>rmodel._c_descriptor).get()))
        self._relations = list(relations)

        self._thisptr.reset(new c_model_definition(c_domains, c_relations))

    def domains(self):
        return self._domains

    def relations(self):
        return self._relations

    def shape(self, relation):
        dids = self._relations[relation][0]
        return tuple(self._domains[did] for did in dids)

    def __reduce__(self):
        args = (self._domains, self._relations)
        return (_reconstruct_model_definition, args)


def _reconstruct_model_definition(domains, relations):
    return model_definition(domains, relations)
