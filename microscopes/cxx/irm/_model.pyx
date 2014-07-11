from microscopes.io.schema_pb2 import CRP

cdef vector[shared_ptr[c_dataview]] get_crelations(relations):
    cdef vector[shared_ptr[c_dataview]] crelations
    for reln in relations:
        if not isinstance(reln, abstract_dataview):
            raise ValueError("bad relation given: " + repr(reln))
        crelations.push_back((<abstract_dataview>reln)._thisptr)
    return crelations

cdef vector[const c_dataview *] get_crelations_raw(relations):
    cdef vector[const c_dataview *] crelations
    for reln in relations:
        if not isinstance(reln, abstract_dataview):
            raise ValueError("bad relation given: " + repr(reln))
        crelations.push_back((<abstract_dataview>reln)._thisptr.get())
    return crelations

cdef vector[set[size_t]] get_cclustering(clustering):
    cdef vector[set[size_t]] cclustering
    cdef set[size_t] buf0
    for cluster in clustering:
        buf0.clear()
        for eid in cluster:
            buf0.insert(int(eid))
        cclustering.push_back(buf0)
    return cclustering

def fill(state s, clusters, relations, rng r):
    cdef vector[vector[set[size_t]]] cclusters
    for clustering in clusters:
        cclusters.push_back(get_cclustering(clustering))
    s._thisptr.get().initialize(cclusters, get_crelations_raw(relations), r._thisptr[0])

def random_initialize(state s, relations, rng r):
    s._thisptr.get().random_initialize(get_crelations_raw(relations), r._thisptr[0])

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

    def get_domain_hp(self, int domain):
        m = CRP()
        raw = str(self._thisptr.get().get_domain_hp(domain))
        m.ParseFromString(raw)
        return {'alpha':m.alpha}

    def set_domain_hp(self, int domain, dict d):
        m = CRP()
        m.alpha = float(d['alpha'])
        self._thisptr.get().set_domain_hp(domain, m.SerializeToString())

    def get_relation_hp(self, int relation):
        raw = str(self._thisptr.get().get_relation_hp(relation))
        return self._relations[relation][1][0].shared_bytes_to_dict(raw)

    def set_relation_hp(self, int relation, dict d):
        cdef hyperparam_bag_t raw = self._relations[relation][1][0].shared_dict_to_bytes(d)
        self._thisptr.get().set_relation_hp(relation, raw)

    def score_assignment(self):
        return self._thisptr.get().score_assignment()

    def score_likelihood(self, rng r):
        assert r
        return self._thisptr.get().score_likelihood(r._thisptr[0]);

    def entity_data_positions(self, int domain, int eid, relations):
        cdef vector[vector[size_t]] cret = self._thisptr.get().entity_data_positions(domain, eid, get_crelations_raw(relations))
        return [[x for x in inner] for inner in cret]

def bind(state s, int domain, relations):
    cdef shared_ptr[c_entity_based_state_object] px
    cdef vector[shared_ptr[c_dataview]] crelations = get_crelations(relations)
    px.reset(new c_bound_state(s._thisptr, domain, crelations))
    cdef entity_based_state_object ret = entity_based_state_object(s.models())
    ret.set_entity(px)
    return ret
