import pytest

from rdkit import Chem


@pytest.mark.parametrize("auto_convert", [False, True])
@pytest.mark.parametrize("prop_value, default_val, expected_val",
                         [("Not Set", None, None),
                          ("Not Set", "fallback", "fallback"),
                          ("value", "fallback", "value"),
                          ("value", None, "value"), (None, "fallback", None),
                          (42, 0, 42), ("Not Set", 0, 0), (42.0, 0.0, 42.0),
                          ("Not Set", 0.0, 0.0)])
def test_get_prop_with_default(auto_convert, prop_value, default_val,
                               expected_val):
    m = Chem.MolFromSmiles("CC")
    prop_key = "test_prop"
    if prop_value != "Not Set":
        if isinstance(prop_value, int):
            m.SetIntProp(prop_key, prop_value)
        elif isinstance(prop_value, float):
            m.SetDoubleProp(prop_key, prop_value)
        else:
            m.SetProp(prop_key, prop_value)

    assert m.GetProp(prop_key,
                     autoConvert=auto_convert,
                     default_val=default_val) == expected_val


@pytest.mark.parametrize("auto_convert", [False, True])
def test_get_prop_no_default_not_set(auto_convert):
    m = Chem.MolFromSmiles("CC")
    with pytest.raises(KeyError):
        m.GetProp("test_prop", autoConvert=auto_convert)


@pytest.mark.parametrize("auto_convert", [False, True])
@pytest.mark.parametrize("prop_value", [None, "fallback", 0, 0.0])
def test_get_prop_no_default_set(auto_convert, prop_value):
    m = Chem.MolFromSmiles("CC")
    prop_key = "test_prop"
    if isinstance(prop_value, int):
        m.SetIntProp(prop_key, prop_value)
    elif isinstance(prop_value, float):
        m.SetDoubleProp(prop_key, prop_value)
    else:
        m.SetProp(prop_key, prop_value)
    assert m.GetProp("test_prop", autoConvert=auto_convert) == prop_value


def test_get_prop_no_keyword():
    m = Chem.MolFromSmiles("CC")
    with pytest.raises(KeyError):
        prop = m.GetProp("test_prop", True)
    true_prop = m.GetProp("test_prop", True, True)
    assert true_prop is True
