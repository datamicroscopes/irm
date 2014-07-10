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

def bind(state s, abstract_dataview view):
    pass
