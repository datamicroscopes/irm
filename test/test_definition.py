from microscopes.irm.definition import model_definition
from microscopes.models import bb, niw

from nose.tools import (
    assert_equals,
    assert_list_equal,
)
import pickle


def test_model_definition_pickle():
    defn = model_definition([10, 12], [((0, 0), bb), ((1, 0), niw(3))])
    bstr = pickle.dumps(defn)
    defn1 = pickle.loads(bstr)
    assert_list_equal(defn.domains(), defn1.domains())
    assert_equals(len(defn.relations()), len(defn1.relations()))
    zipped_relations = zip(defn.relations(), defn1.relations())
    for (dids, model), (dids1, model1) in zipped_relations:
        assert_list_equal(list(dids), list(dids1))
        assert_equals(model.name(), model1.name())
