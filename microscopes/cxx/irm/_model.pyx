cdef class state:
    def __cinit__(self, domains, relations):
        self._domains = list(domains)
        self._relations = list(relations)

        cdef vector[size_t] cdomains
        for d in domains:
            cdomains.push_back(int(d))

        cdef vector[c_state.relation_t] crelations
        cdef vector[size_t] buf0
        for dlist, model in relations:
            buf0.clear()
            for d in dlist:
                buf0.push_back(int(d))
            _, c_m = model
            crelations.push_back(c_state.relation_t(buf0, (<factory>c_m).new_cmodel()))

        self._thisptr.reset(new c_state(cdomains, crelations))

    def models(self):
        return [m for _, m in self._relations]

def bind(state s, int domain, relations):
    cdef shared_ptr[c_entity_based_state_object] px
    cdef vector[shared_ptr[c_dataview]] crelations
    for reln in relations:
        if not isinstance(reln, abstract_dataview):
            raise ValueError("bad relation given: " + repr(reln))
        crelations.push_back((<abstract_dataview>reln)._thisptr)
    px.reset(new c_bound_state(s._thisptr, domain, crelations))
    cdef entity_based_state_object ret = entity_based_state_object(s.models())
    ret.set_entity(px)
    return ret
