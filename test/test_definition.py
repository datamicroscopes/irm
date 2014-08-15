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
    assert_equals(defn.relations(), defn1.relations())
    zipped_models = zip(defn.relation_models(), defn1.relation_models())
    for model, model1 in zipped_models:
        assert_equals(model.name(), model1.name())
    # XXX(stephentu): check hyperpriors
