import pickle
import copy
import itertools as it

from microscopes.irm.definition import model_definition
from microscopes.irm import model
from microscopes.irm.testutil import toy_dataset
from microscopes.models import bb
from microscopes.common.rng import rng
from microscopes.common.relation.dataview import numpy_dataview

from nose.tools import (
    assert_equals,
    assert_almost_equals,
    assert_is_none,
    assert_is_not,
)
from distributions.tests.util import assert_close


def _assert_structure_equals(defn, s1, s2, views, r):
    assert_equals(s1.ndomains(), s2.ndomains())
    assert_equals(s1.nrelations(), s2.nrelations())
    for did in xrange(s1.ndomains()):
        assert_equals(s1.nentities(did), s2.nentities(did))
        assert_equals(s1.ngroups(did), s2.ngroups(did))
        assert_equals(s1.assignments(did),
                      s2.assignments(did))
        assert_equals(set(s1.groups(did)),
                      set(s2.groups(did)))
        assert_close(s1.get_domain_hp(did),
                     s2.get_domain_hp(did))
        assert_almost_equals(s1.score_assignment(did),
                             s2.score_assignment(did))
    for rid in xrange(s1.nrelations()):
        assert_close(s1.get_relation_hp(rid),
                     s2.get_relation_hp(rid))
        dids = defn.relations()[rid]
        groups = [s1.groups(did) for did in dids]
        for gids in it.product(*groups):
            ss1 = s1.get_suffstats(rid, gids)
            ss2 = s2.get_suffstats(rid, gids)
            if ss1 is None:
                assert_is_none(ss2)
            else:
                assert_close(ss1, ss2)
    assert_almost_equals(s1.score_likelihood(r),
                         s2.score_likelihood(r))
    before = list(s1.assignments(0))
    bound = model.bind(s1, 0, views)
    gid = bound.remove_value(0, r)
    assert_equals(s1.assignments(0)[0], -1)
    assert_equals(before, s2.assignments(0))
    bound.add_value(gid, 0, r)  # restore


def test_state_pickle():
    defn = model_definition([5], [((0, 0), bb)])
    r = rng()
    relations = toy_dataset(defn)
    views = map(numpy_dataview, relations)
    s1 = model.initialize(defn, views, r)
    s2 = pickle.loads(pickle.dumps(s1))
    _assert_structure_equals(defn, s1, s2, views, r)


def test_state_copy():
    defn = model_definition([5], [((0, 0), bb)])
    r = rng()
    relations = toy_dataset(defn)
    views = map(numpy_dataview, relations)
    s1 = model.initialize(defn, views, r)
    s2 = copy.copy(s1)
    assert_is_not(s1, s2)
    _assert_structure_equals(defn, s1, s2, views, r)

    s2 = copy.deepcopy(s1)
    assert_is_not(s1, s2)
    _assert_structure_equals(defn, s1, s2, views, r)
