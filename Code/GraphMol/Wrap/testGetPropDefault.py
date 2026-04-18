import pytest

from rdkit import Chem


def test_get_prop_key_error():
  m = Chem.MolFromSmiles("CC")
  with pytest.raises(KeyError):
    m.GetProp("missing_prop")


def test_get_prop_default_none():
  m = Chem.MolFromSmiles("CC")
  assert m.GetProp("missing_prop", default_val=None) is None


def test_get_prop_default_other():
  m = Chem.MolFromSmiles("CC")
  assert m.GetProp("missing_prop", default_val="fallback") == "fallback"


def test_get_prop_present():
  m = Chem.MolFromSmiles("CC")
  m.SetProp("present_prop", "value")
  assert m.GetProp("present_prop", default_val="fallback") == "value"


def test_get_prop_default_auto_convert():
  m = Chem.MolFromSmiles("CC")
  assert m.GetProp("missing_prop", autoConvert=True, default_val="fallback") == "fallback"


def test_get_prop_default_auto_convert_none():
  m = Chem.MolFromSmiles("CC")
  assert m.GetProp("missing_prop", autoConvert=True, default_val=None) is None


def test_get_prop_present_auto_convert():
  m = Chem.MolFromSmiles("CC")
  m.SetIntProp("int_prop", 42)
  assert m.GetProp("int_prop", autoConvert=True, default_val=0) == 42


def test_get_prop_auto_convert_key_error_when_missing():
  m = Chem.MolFromSmiles("CC")
  with pytest.raises(KeyError):
    m.GetProp("missing_prop", autoConvert=True)
