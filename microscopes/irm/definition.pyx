# cython: embedsignature=True


# python imports
from microscopes.models import model_descriptor
from microscopes.common import validator
import operator as op


# Why isn't this in the standard library?
def _compose(f, g):
    return lambda x: g(f(x))


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
        self._relations = []
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

            if hasattr(rmodel, '__len__'):
                validator.validate_len(rmodel, 2)
                validator.validate_type(rmodel[0], model_descriptor)
                validator.validate_type(rmodel[1], dict)
                rmodel, hp = rmodel
            else:
                validator.validate_type(rmodel, model_descriptor)
                rmodel, hp = rmodel, rmodel.default_hyperpriors()
            self._relations.append((rdomains, (rmodel, hp)))

            c_relations.push_back(
                c_relation_definition(
                    c_rdomains,
                    (<_base>rmodel._c_descriptor).get()))

        self._thisptr.reset(new c_model_definition(c_domains, c_relations))

    def domains(self):
        return self._domains

    def domain_hyperpriors(self):
        raise RuntimeError("not implemented")

    def relations(self):
        return map(op.itemgetter(0), self._relations)

    def relation_models(self):
        fn = _compose(op.itemgetter(1), op.itemgetter(0))
        return map(fn, self._relations)

    def relation_hyperpriors(self):
        fn = _compose(op.itemgetter(1), op.itemgetter(1))
        return map(fn, self._relations)

    def shape(self, relation):
        dids = self.relations()[relation]
        return tuple(self._domains[did] for did in dids)

    def __reduce__(self):
        args = (self._domains, self._relations)
        return (_reconstruct_model_definition, args)


def _reconstruct_model_definition(domains, relations):
    return model_definition(domains, relations)
