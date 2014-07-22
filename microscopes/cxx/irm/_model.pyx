from microscopes.io.schema_pb2 import CRP

def _check_relations(relations):
    for reln in relations:
        if not isinstance(reln, abstract_dataview):
            raise ValueError(
                "expected abstract_dataview, got {}".format(repr(reln)))

cdef vector[shared_ptr[c_dataview]] get_crelations(relations):
    cdef vector[shared_ptr[c_dataview]] crelations
    for reln in relations:
        crelations.push_back((<abstract_dataview>reln)._thisptr)
    return crelations

cdef vector[const c_dataview *] get_crelations_raw(relations):
    cdef vector[const c_dataview *] crelations
    for reln in relations:
        crelations.push_back((<abstract_dataview>reln)._thisptr.get())
    return crelations

cdef class state:
    def __cinit__(self, model_definition defn, **kwargs):
        self._defn = defn

        # note: python cannot overload __cinit__(), so we 
        # use kwargs to handle both the random initialization case and
        # the deserialize from string case
        if not (('data' in kwargs) ^ ('bytes' in kwargs)):
            raise ValueError("need exaclty one of `data' or `bytes'")

        cdef vector[hyperparam_bag_t] c_cluster_hps
        cdef vector[hyperparam_bag_t] c_relation_hps
        cdef vector[vector[size_t]] c_domain_assignments
        cdef vector[size_t] c_assignment

        if 'data' in kwargs:
            # handle the random initialization case
            data = list(kwargs['data'])
            if len(data) != len(defn._relations):
                raise ValueError("expecting {} relations, got {}".format(
                    len(defn._relations), len(data)))
            _check_relations(data)

            if 'r' not in kwargs:
                raise ValueError("need parameter `r'")
            r = kwargs['r']
            if not isinstance(r, rng):
                raise ValueError("need prng for parameter `r'")

            if 'cluster_hps' in kwargs:
                cluster_hps = list(kwargs['cluster_hps'])
            else:
                cluster_hps = [{'alpha':1.}]*len(defn._domains)

            def make_cluster_hp_bytes(cluster_hp):
                m = CRP()
                m.alpha = cluster_hp['alpha']
                return m.SerializeToString()
            for cluster_hp in cluster_hps:
                c_cluster_hps.push_back(make_cluster_hp_bytes(cluster_hp))

            if 'relation_hps' in kwargs:
                relation_hps = list(kwargs['relation_hps'])
                if len(relation_hps) != len(defn._relations):
                    raise ValueError("expecting {} models, got {}".format(
                        len(relation_hps), len(defn._relations)))
            else:
                relation_hps = [m.default_params() for _, m in defn._relations] 

            relation_hps_bytes = [
                m.py_desc().shared_dict_to_bytes(hp) \
                    for hp, (_, m) in zip(relation_hps, defn._relations)]
            for s in relation_hps_bytes:
                c_relation_hps.push_back(s)

            if 'domain_assignments' in kwargs: 
                domain_assignments = list(kwargs['domain_assignments'])
                if len(domain_assignments) != len(defn._domains):
                    raise ValueError(
                        "expecting {} domain assignments, got {}".format(
                            len(defn._domains), len(domain_assignments)))
                for did, assignment in enumerate(domain_assignments):
                    assignment = list(assignment)
                    if not len(assignment):
                        c_domain_assignments.push_back(vector[size_t]())
                    else:
                        if len(assignment) != len(defn._domains[did]):
                            raise ValueError(
                                "expecting {} assignments, got {}".format(
                                    len(defn._domains[did]),
                                    len(assignment)))
                        c_assignment.clear()
                        for i in assignment:
                            c_assignment.push_back(i)
                        c_domain_assignments.push_back(c_assignment)
            else:
                c_domain_assignments.resize(len(defn._domains))

            self._thisptr = c_initialize(
                defn._thisptr.get()[0],
                c_cluster_hps,
                c_relation_hps,
                c_domain_assignments,
                get_crelations_raw(data),
                (<rng>r)._thisptr[0])

        else:
            # handle the deserialize case
            self._thisptr = c_deserialize(
                defn._thisptr.get()[0],
                kwargs['bytes'])

        if self._thisptr.get() == NULL:
            raise RuntimeError("could not properly construct state")

    def models(self):
        return [m for _, m in self._defn._relations]

    def ndomains(self):
        return self._thisptr.get().ndomains()

    def nentities(self, int domain):
        return self._thisptr.get().nentities(domain)

    def groups(self, int domain):
        return [g for g in self._thisptr.get().groups(domain)]

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
        desc = self._defn._relations[relation][1].py_desc()
        return desc.shared_bytes_to_dict(raw)

    def set_relation_hp(self, int relation, dict d):
        desc = self._defn._relations[relation][1].py_desc()
        cdef hyperparam_bag_t raw = desc.shared_dict_to_bytes(d)
        self._thisptr.get().set_relation_hp(relation, raw)

    def get_suffstats(self, int relation, gids):
        desc = self._defn._relations[relation][1].py_desc()
        cdef suffstats_bag_t raw
        cdef vector[size_t] cgids
        for g in gids:
            cgids.push_back(int(g))
        cdef cbool found = self._thisptr.get().get_suffstats(relation, cgids, raw)
        if not found:
            return None
        else:
            return desc.group_bytes_to_dict(raw)

    def score_assignment(self, int domain):
        return self._thisptr.get().score_assignment(domain)

    def score_likelihood(self, rng r):
        assert r
        return self._thisptr.get().score_likelihood(r._thisptr[0]);

    def entity_data_positions(self, int domain, int eid, relations):
        cdef vector[vector[size_t]] cret = self._thisptr.get().entity_data_positions(domain, eid, get_crelations_raw(relations))
        return [[x for x in inner] for inner in cret]

    def serialize(self):
        return self._thisptr.get().serialize()

def bind(state s, int domain, relations):
    cdef shared_ptr[c_entity_based_state_object] px
    cdef vector[shared_ptr[c_dataview]] crelations = get_crelations(relations)
    px.reset(new c_model(s._thisptr, domain, crelations))
    cdef entity_based_state_object ret = entity_based_state_object(s.models())
    ret.set_entity(px)
    return ret

def initialize(model_definition defn, data, rng r, **kwargs):
    return state(defn=defn, data=data, r=r, **kwargs)

def deserialize(model_definition defn, bytes, **kwargs):
    return state(defn=defn, bytes=bytes, **kwargs)
