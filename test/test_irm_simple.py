from microscopes.cxx.irm.model import state, bind
from microscopes.cxx.common.sparse_ndarray.dataview import numpy_dataview
from microscopes.cxx.models import bb

import numpy as np
import numpy.ma as ma
import itertools as it

def test_simple():
    domains = [5, 6]

    relations = [ ((0, 1), bb) ]

    raw_data = [ma.array(
                np.random.choice([False,True], size=(domains[0], domains[1])),
                mask=np.random.choice([False,True], size=(domains[0], domains[1])))]

    data = map(numpy_dataview, raw_data)

    s = state(domains, relations)
    bound_s0 = bind(s, 0, data)
    bound_s1 = bind(s, 1, data)

    assert s and bound_s0 and bound_s1

    def entity_data_positions(domain, eid):
        def f(domains, reln):
          for pos0 in xrange(reln.shape[0]):
              for pos1 in xrange(reln.shape[1]):
                  if reln.mask[pos0, pos1]:
                      continue
                  if (domains[0] == domain and pos0 == eid) or \
                     (domains[1] == domain and pos1 == eid):
                       yield [pos0, pos1]
        return list(it.chain.from_iterable(f(domains, reln) for (domains, _), reln in zip(relations, raw_data)))

    for did, nentities in enumerate(domains):
        for eid in xrange(nentities):
            a = entity_data_positions(did, eid)
            b = s.entity_data_positions(did, eid, data)
            print a, b
            assert sorted(a) == sorted(b)
