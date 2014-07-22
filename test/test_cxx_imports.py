def test_import_irm():
    from microscopes.cxx.irm.model import state, bind
    from microscopes.irm.definition import model_definition
    assert state and bind and model_definition
