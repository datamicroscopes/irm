from microscopes.cxx.irm.model import initialize, bind
from microscopes.cxx.common.relation.dataview import \
    numpy_dataview, sparse_2d_dataview
from microscopes.models import bb
from microscopes.irm.definition import model_definition
from microscopes.cxx.common.rng import rng

import numpy as np
import numpy.ma as ma
import itertools as it
import operator as op

from scipy.sparse import coo_matrix

def test_simple():
    domains = [5, 6]

    relations = [((0, 1), bb)]

    raw_data = [ma.array(
                np.random.choice([False,True], size=(domains[0], domains[1])),
                mask=np.random.choice([False,True], size=(domains[0], domains[1])))]

    def csr(raw):
        n, m = raw.shape
        def indices():
            for i, j in it.product(range(n), range(m)):
                if not raw.mask[i, j]:
                    yield i, j
        data = [raw[i, j] for i, j in indices()]
        i = list(map(op.itemgetter(0), indices()))
        j = list(map(op.itemgetter(1), indices()))
        return coo_matrix((data, (i, j)), shape=raw.shape).tocsr()

    defn = model_definition(domains, relations)
    data = map(numpy_dataview, raw_data)
    sparse_data = map(sparse_2d_dataview, map(csr, raw_data))

    r = rng()

    s = initialize(defn, data, r=r)
    assert s and bind(s, 0, data) and bind(s, 1, data)

    s1 = initialize(defn, sparse_data, r=r)
    assert s1 and bind(s1, 0, sparse_data) and bind(s1, 1, sparse_data)

    def entity_data_positions(domain, eid):
        def f(domains, reln):
          for pos0 in xrange(reln.shape[0]):
              for pos1 in xrange(reln.shape[1]):
                  if reln.mask[pos0, pos1]:
                      continue
                  if (domains[0] == domain and pos0 == eid) or \
                     (domains[1] == domain and pos1 == eid):
                       yield [pos0, pos1]
        return list(it.chain.from_iterable(
            f(domains, reln) for (domains, _), reln in zip(relations, raw_data)))

    def test(s):
        for did, nentities in enumerate(domains):
            for eid in xrange(nentities):
                a = entity_data_positions(did, eid)
                b = s.entity_data_positions(did, eid, data)
                assert sorted(a) == sorted(b)

    test(s)
    test(s1)
