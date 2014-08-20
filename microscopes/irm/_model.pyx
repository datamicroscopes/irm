# cython: embedsignature=True


# python imports
from microscopes.common._rng import rng
from microscopes.common.relation._dataview import abstract_dataview
from microscopes.irm.definition import model_definition
from microscopes.common import validator
from microscopes.io.schema_pb2 import CRP


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

        valid_kwargs = ('data', 'bytes', 'r',
                        'cluster_hps', 'relation_hps', 'domain_assignments',)
        validator.validate_kwargs(kwargs, valid_kwargs)

        cdef vector[hyperparam_bag_t] c_cluster_hps
        cdef vector[hyperparam_bag_t] c_relation_hps
        cdef vector[vector[size_t]] c_domain_assignments
        cdef vector[size_t] c_assignment

        if 'data' in kwargs:
            # handle the random initialization case
            data = list(kwargs['data'])
            validator.validate_len(data, len(defn.relations()), "data")
            for rid, view in enumerate(data):
                validator.validate_type(view, abstract_dataview)
                expected = defn.shape(rid)
                actual = view.shape()
                if expected != actual:
                    msg = "expected view of shape {}, got shape {}"
                    raise ValueError(msg.format(expected, actual))

            if 'r' not in kwargs:
                raise ValueError("need parameter `r'")
            r = kwargs['r']
            validator.validate_type(r, rng, "r")

            if 'cluster_hps' in kwargs:
                cluster_hps = list(kwargs['cluster_hps'])
                validator.validate_len(
                    cluster_hps, len(defn.domains()), "cluster_hps")
            else:
                cluster_hps = [{'alpha': 1.}] *len(defn.domains())

            def make_cluster_hp_bytes(cluster_hp):
                m = CRP()
                m.alpha = cluster_hp['alpha']
                return m.SerializeToString()
            for cluster_hp in cluster_hps:
                c_cluster_hps.push_back(make_cluster_hp_bytes(cluster_hp))

            if 'relation_hps' in kwargs:
                relation_hps = list(kwargs['relation_hps'])
                validator.validate_len(
                    relation_hps, len(defn.relations()), "relation_hps")
            else:
                models = defn.relation_models()
                relation_hps = [m.default_hyperparams() for m in models]

            relation_hps_bytes = [
                m.py_desc().shared_dict_to_bytes(hp)
                for hp, m in zip(relation_hps, defn.relation_models())]
            for s in relation_hps_bytes:
                c_relation_hps.push_back(s)

            if 'domain_assignments' in kwargs:
                domain_assignments = list(kwargs['domain_assignments'])
                validator.validate_len(
                    domain_assignments,
                    len(defn.domains()),
                    "domain_assignments")
                for did, assignment in enumerate(domain_assignments):
                    assignment = list(assignment)
                    if not len(assignment):
                        c_domain_assignments.push_back(vector[size_t]())
                    else:
                        validator.validate_len(assignment, defn.domains()[did])
                        c_assignment.clear()
                        for i in assignment:
                            validator.validate_nonnegative(i)
                            c_assignment.push_back(i)
                        c_domain_assignments.push_back(c_assignment)
            else:
                c_domain_assignments.resize(len(defn.domains()))

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
        return self._defn.relation_models()

    def _validate_did(self, did, param_name=None):
        validator.validate_in_range(did, self.ndomains(), param_name)

    def _validate_rid(self, rid, param_name=None):
        validator.validate_in_range(rid, self.nrelations(), param_name)

    def _validate_eid(self, did, eid):
        validator.validate_in_range(eid, self.nentities(did))

    def _validate_gid(self, did, gid):
        if not self.isactivegroup(did, gid):
            raise ValueError("invalid gid")

    def ndomains(self):
        return len(self._defn.domains())

    def nrelations(self):
        return len(self._defn.relations())

    def nentities(self, int domain):
        self._validate_did(domain, "domain")
        return self._thisptr.get().nentities(domain)

    def ngroups(self, int domain):
        self._validate_did(domain, "domain")
        return self._thisptr.get().ngroups(domain)

    def assignments(self, int domain):
        self._validate_did(domain, "domain")
        return self._thisptr.get().assignments(domain)

    def groups(self, int domain):
        self._validate_did(domain, "domain")
        return [g for g in self._thisptr.get().groups(domain)]

    def empty_groups(self, int domain):
        self._validate_did(domain, "domain")
        return list(self._thisptr.get().empty_groups(domain))

    def isactivegroup(self, int domain, int gid):
        self._validate_did(domain, "domain")
        return self._thisptr.get().isactivegroup(domain, gid)

    def groupsize(self, int domain, int gid):
        self._validate_gid(domain, gid)
        return self._thisptr.get().groupsize(domain, gid)

    def get_domain_hp(self, int domain):
        self._validate_did(domain, "domain")
        m = CRP()
        raw = str(self._thisptr.get().get_domain_hp(domain))
        m.ParseFromString(raw)
        return {'alpha': m.alpha}

    def set_domain_hp(self, int domain, dict d):
        self._validate_did(domain, "domain")
        m = CRP()
        m.alpha = float(d['alpha'])
        self._thisptr.get().set_domain_hp(domain, m.SerializeToString())

    def get_relation_hp(self, int relation):
        self._validate_rid(relation, "relation")
        raw = str(self._thisptr.get().get_relation_hp(relation))
        desc = self._defn.relation_models()[relation].py_desc()
        return desc.shared_bytes_to_dict(raw)

    def set_relation_hp(self, int relation, dict d):
        self._validate_rid(relation, "relation")
        desc = self._defn.relation_models()[relation].py_desc()
        cdef hyperparam_bag_t raw = desc.shared_dict_to_bytes(d)
        self._thisptr.get().set_relation_hp(relation, raw)

    def get_suffstats(self, int relation, gids):
        self._validate_rid(relation, "relation")
        desc = self._defn.relation_models()[relation].py_desc()
        arity = len(self._defn.relations()[relation])
        validator.validate_len(gids, arity, "gids")
        cdef suffstats_bag_t raw
        cdef vector[size_t] cgids
        for g in gids:
            # XXX(stephentu): need to validate this, but it's OK for now since
            # the C++ API will just return false if not found
            cgids.push_back(int(g))
        cdef cbool found = (
            self._thisptr.get().get_suffstats(relation, cgids, raw)
        )
        if not found:
            return None
        else:
            return desc.group_bytes_to_dict(raw)

    def score_assignment(self, int domain):
        self._validate_did(domain, "domain")
        return self._thisptr.get().score_assignment(domain)

    def score_likelihood(self, rng r):
        validator.validate_not_none(r)
        return self._thisptr.get().score_likelihood(r._thisptr[0])

    # XXX(stephentu): this is used for debugging and should be removed
    def entity_data_positions(self, int domain, int eid, relations):
        self._validate_eid(domain, eid)
        cdef vector[vector[size_t]] cret = (
            self._thisptr.get().entity_data_positions(
                domain, eid, get_crelations_raw(relations))
        )
        return [[x for x in inner] for inner in cret]

    def serialize(self):
        return self._thisptr.get().serialize()

    def __reduce__(self):
        return (_reconstruct_state, (self._defn, self.serialize()))

    # XXX(stephentu): expose more methods


def bind(state s, int domain, relations):
    s._validate_did(domain, "domain")
    validator.validate_len(relations, s.nrelations(), "relations")
    cdef shared_ptr[c_entity_based_state_object] px
    cdef vector[shared_ptr[c_dataview]] crelations = get_crelations(relations)
    px.reset(new c_model(s._thisptr, domain, crelations))
    cdef entity_based_state_object ret = entity_based_state_object(s.models())
    ret.set_entity(px)
    ret._refs = relations
    return ret


def initialize(model_definition defn, data, rng r, **kwargs):
    return state(defn=defn, data=data, r=r, **kwargs)


def deserialize(model_definition defn, bytes):
    return state(defn=defn, bytes=bytes)


def _reconstruct_state(defn, bytes):
    return deserialize(defn, bytes)
