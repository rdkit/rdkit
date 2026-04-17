import pytest

from rdkit import Chem

PROP_KEY = "test_prop"


@pytest.mark.parametrize("prop_value, expected_val",
                         [("value", "value"),
                          (0, 0),
                          (42.0, 42.0),
                          (True, True)])
def test_get_prop_current_behavior(prop_value, expected_val):
    """
    This test checks the current behavior of GetProp. We don't want to change
    current behavior to avoid breaking existing code.
    """
    m = Chem.MolFromSmiles("CC")
    if isinstance(prop_value, bool):
        m.SetBoolProp(PROP_KEY, prop_value)
    elif isinstance(prop_value, int):
        m.SetIntProp(PROP_KEY, prop_value)
    elif isinstance(prop_value, float):
        m.SetDoubleProp(PROP_KEY, prop_value)
    else:
        m.SetProp(PROP_KEY, prop_value)

    assert m.GetProp(PROP_KEY, True) == expected_val


@pytest.mark.parametrize("autoConvert", [False, True])
def test_get_prop_current_behavior_not_set(autoConvert):
    """
    This test checks the current behavior of GetProp when the property is not
    set. We expect a KeyError to be raised.
    """
    m = Chem.MolFromSmiles("CC")
    with pytest.raises(KeyError):
        m.GetProp(PROP_KEY, autoConvert)


@pytest.mark.parametrize("prop_value, default_val, expected_val",
                         [("Not Set", "fallback", "fallback"),
                          ("Not Set", None, None),
                          ("Not Set", 0, 0),
                          ("Not Set", 42.0, 42.0),
                          ("Not Set", True, True),
                          ("value", "fallback", "value"),
                          (0, "fallback", 0),
                          (42.0, "fallback", 42.0),
                          (True, "fallback", True)])
def test_get_prop_with_default(prop_value, default_val, expected_val):
    """
    This test describes the default_val would work. This test uses a keyword.
    """
    m = Chem.MolFromSmiles("CC")
    if prop_value != "Not Set":
        if isinstance(prop_value, bool):
            m.SetBoolProp(PROP_KEY, prop_value)
        elif isinstance(prop_value, int):
            m.SetIntProp(PROP_KEY, prop_value)
        elif isinstance(prop_value, float):
            m.SetDoubleProp(PROP_KEY, prop_value)
        else:
            m.SetProp(PROP_KEY, prop_value)

    assert m.GetPropIfPresent(PROP_KEY, default_val=default_val) == expected_val


@pytest.mark.parametrize("prop_value, default_val, expected_val",
                         [("Not Set", "fallback", "fallback"),
                          ("Not Set", None, None),
                          ("Not Set", 0, 0),
                          ("Not Set", 42.0, 42.0),
                          ("Not Set", True, True),
                          ("value", "fallback", "value"),
                          (0, "fallback", 0),
                          (42.0, "fallback", 42.0),
                          (True, "fallback", True)])
def test_get_prop_with_default(prop_value, default_val, expected_val):
    """
    This test describes the default_val would work. This test does not use
    a keyword to check if the default_val would work without using a keyword.
    """
    m = Chem.MolFromSmiles("CC")
    if prop_value != "Not Set":
        if isinstance(prop_value, bool):
            m.SetBoolProp(PROP_KEY, prop_value)
        elif isinstance(prop_value, int):
            m.SetIntProp(PROP_KEY, prop_value)
        elif isinstance(prop_value, float):
            m.SetDoubleProp(PROP_KEY, prop_value)
        else:
            m.SetProp(PROP_KEY, prop_value)

    assert m.GetPropIfPresent(PROP_KEY, default_val) == expected_val


def test_get_prop_no_default_not_set():
    """
    When a default value is not provided and the property is not set, we
    expect None to be returned instead of raising an error.
    """
    m = Chem.MolFromSmiles("CC")
    assert m.GetPropIfPresent(PROP_KEY) is None


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__]))
