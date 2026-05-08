import pytest

from rdkit import Chem
from rdkit.Chem import AllChem


# ---------------------------------------------------------------------------
# GetProp with default
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("auto_convert", [False, True])
@pytest.mark.parametrize("default", [None, "fallback", 0, 0.0])
def test_get_prop_with_default_missing(auto_convert, default):
    """Default is returned when the property is not set."""
    m = Chem.MolFromSmiles("CC")
    assert m.GetProp("test_prop", autoConvert=auto_convert,
                     default=default) == default


@pytest.mark.parametrize("auto_convert", [False, True])
@pytest.mark.parametrize("prop_value, default, unconverted_value", [
    ("value", "fallback", "value"),
    ("value", None, "value"),
    (42, 0, '42'),
    (42.0, 0.0, '42'),
])
def test_get_prop_with_default_present(auto_convert, prop_value, default, unconverted_value):
    """The stored property value is returned, not the default."""
    m = Chem.MolFromSmiles("CC")
    if isinstance(prop_value, int):
        m.SetIntProp("test_prop", prop_value)
    elif isinstance(prop_value, float):
        m.SetDoubleProp("test_prop", prop_value)
    else:
        m.SetProp("test_prop", prop_value)
    expected = prop_value if auto_convert else unconverted_value
    assert m.GetProp("test_prop", autoConvert=auto_convert,
                     default=default) == expected


@pytest.mark.parametrize("auto_convert", [False, True])
def test_get_prop_no_default_not_set(auto_convert):
    """KeyError is raised when the property is missing and no default is given."""
    m = Chem.MolFromSmiles("CC")
    with pytest.raises(KeyError):
        m.GetProp("test_prop", autoConvert=auto_convert)


@pytest.mark.parametrize("auto_convert", [False, True])
@pytest.mark.parametrize("prop_value, unconverted_value", [
    ("hello", "hello"),
    (42, "42"),
    (42.0, "42"),
])
def test_get_prop_no_default_set(auto_convert, prop_value, unconverted_value):
    """GetProp without default returns the stored value normally."""
    m = Chem.MolFromSmiles("CC")
    if isinstance(prop_value, int):
        m.SetIntProp("test_prop", prop_value)
    elif isinstance(prop_value, float):
        m.SetDoubleProp("test_prop", prop_value)
    else:
        m.SetProp("test_prop", prop_value)
    expected = prop_value if auto_convert else unconverted_value
    assert m.GetProp("test_prop", autoConvert=auto_convert) == expected


def test_get_prop_default_positional():
    """default can be supplied positionally as the third argument."""
    m = Chem.MolFromSmiles("CC")
    with pytest.raises(KeyError):
        m.GetProp("test_prop", True)
    default_val = m.GetProp("test_prop", True, True)
    assert default_val is True


# ---------------------------------------------------------------------------
# GetProp with default on Atom, Bond, Conformer, SubstanceGroup
# ---------------------------------------------------------------------------

def test_get_prop_default_on_atom():
    m = Chem.MolFromSmiles("CC")
    atom = m.GetAtomWithIdx(0)
    assert atom.GetProp("missing", default="x") == "x"
    with pytest.raises(KeyError):
        atom.GetProp("missing")
    atom.SetProp("p", "v")
    assert atom.GetProp("p", default="x") == "v"


def test_get_prop_default_on_bond():
    m = Chem.MolFromSmiles("CC")
    bond = m.GetBondWithIdx(0)
    assert bond.GetProp("missing", default="x") == "x"
    with pytest.raises(KeyError):
        bond.GetProp("missing")
    bond.SetProp("p", "v")
    assert bond.GetProp("p", default="x") == "v"


def test_get_prop_default_on_conformer():
    m = Chem.MolFromSmiles("CC")
    AllChem.EmbedMolecule(m, randomSeed=42)
    conf = m.GetConformer()
    assert conf.GetProp("missing", default="x") == "x"
    with pytest.raises(KeyError):
        conf.GetProp("missing")
    conf.SetProp("p", "v")
    assert conf.GetProp("p", default="x") == "v"


def test_get_prop_default_on_substance_group():
    m = Chem.MolFromSmiles("CC")
    sg = Chem.CreateMolSubstanceGroup(m, "SRU")
    assert sg.GetProp("missing", default="x") == "x"
    with pytest.raises(KeyError):
        sg.GetProp("missing")
    sg.SetProp("p", "v")
    assert sg.GetProp("p", default="x") == "v"


# ---------------------------------------------------------------------------
# Typed getters with default — Mol
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("set_fn, get_fn, val, default", [
    ("SetIntProp", "GetIntProp", 42, 0),
    ("SetDoubleProp", "GetDoubleProp", 3.14, 0.0),
    ("SetBoolProp", "GetBoolProp", True, False),
    ("SetUnsignedProp", "GetUnsignedProp", 7, 0),
])
def test_typed_getter_with_default_present(set_fn, get_fn, val, default):
    """Stored typed value is returned even when a default is provided."""
    m = Chem.MolFromSmiles("CC")
    getattr(m, set_fn)("p", val)
    assert getattr(m, get_fn)("p", default=default) == val


@pytest.mark.parametrize("get_fn, default", [
    ("GetIntProp", 0),
    ("GetDoubleProp", 0.0),
    ("GetBoolProp", False),
    ("GetUnsignedProp", 0),
])
def test_typed_getter_with_default_missing(get_fn, default):
    """Default is returned when the property is not set."""
    m = Chem.MolFromSmiles("CC")
    assert getattr(m, get_fn)("missing", default=default) == default


@pytest.mark.parametrize("get_fn", ["GetIntProp", "GetDoubleProp", "GetBoolProp", "GetUnsignedProp"])
def test_typed_getter_no_default_raises(get_fn):
    """Without a default, KeyError is raised for a missing property."""
    m = Chem.MolFromSmiles("CC")
    with pytest.raises(KeyError):
        getattr(m, get_fn)("missing")


def test_typed_getter_default_wrong_type_raises_value_error():
    """GetIntProp(default=x) raises ValueError when the key exists but holds the wrong type."""
    m = Chem.MolFromSmiles("CC")
    m.SetDoubleProp("dbl_key", 3.14)
    with pytest.raises(ValueError):
        m.GetIntProp("dbl_key", default=0)


if __name__ == '__main__':
    import sys
    sys.exit(pytest.main([__file__]))
