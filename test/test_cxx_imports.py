def test_import_irm():
    from microscopes.cxx.irm.model import state, bind
    assert state and bind
