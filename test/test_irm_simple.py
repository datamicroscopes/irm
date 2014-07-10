from microscopes.cxx.irm.model import state, bind
from microscopes.cxx.common.sparse_ndarray.dataview import numpy_dataview
from microscopes.cxx.models import bb

import numpy as np
import numpy.ma as ma

def test_simple():
    domains = [5, 6]

    relations = [ ((0, 1), bb) ]

    data = [
        numpy_dataview(
              ma.array(
                np.random.choice([False,True], size=(domains[0], domains[1])),
                mask=np.random.choice([False,True], size=(domains[0], domains[1]))))
    ]

    s = state(domains, relations)
    bound_s0 = bind(s, 0, data)
    bound_s1 = bind(s, 1, data)

    assert s and bound_s0 and bound_s1
