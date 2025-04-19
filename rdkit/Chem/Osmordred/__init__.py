try:
    from  rdkit.Chem.rdOsmordred import *
except ImportError as e:
    raise ImportError("Osmordred not installed") from e
