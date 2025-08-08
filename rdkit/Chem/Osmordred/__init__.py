# Copyright (c) 2020-2023 Greg Landrum and other RDKit contributors
#  All rights reserved.
#
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
"""
Osmordred descriptors for RDKit
"""

from rdkit import RDConfig
import os
import sys

# Add the rdOsmordred module to the path
sys.path.insert(0, os.path.join(RDConfig.RDContribDir, 'Osmordred'))

try:
    from rdkit.Chem import rdOsmordred
except ImportError:
    try:
        import rdOsmordred
    except ImportError:
        raise ImportError("rdOsmordred module not found. Please ensure Osmordred is properly built.")

import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
from typing import List, Optional, Tuple, Union
import warnings

def CalcOsmordred(smiles: str,  names: bool = False,
                  mynames: Optional[List[str]] = None) -> Union[np.ndarray, Tuple[np.ndarray, List[str]]]:
    """
    Calculate Osmordred descriptors for a single SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        names (bool): Whether to return descriptor names along with values
        mynames (List[str], optional): Custom descriptor names list
        
    Returns:
        Union[np.ndarray, Tuple[np.ndarray, List[str]]]: 
            - If names=False: numpy array of descriptor values
            - If names=True: tuple of (descriptor_values, descriptor_names)
    """

    doExEstate = True

    from rdkit import Chem
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    results = []
    descriptor_names = []

    # Define descriptors with names
    descriptors = [
        ("ABCIndex", rdOsmordred.CalcABCIndex),
        ("AcidBase", rdOsmordred.CalcAcidBase),
        ("AdjacencyMatrix", lambda mol: rdOsmordred.CalcAdjacencyMatrix(mol)),
        ("Aromatic", rdOsmordred.CalcAromatic),
        ("AtomCount", lambda mol: rdOsmordred.CalcAtomCount(mol)),
        ("Autocorrelation", rdOsmordred.CalcAutocorrelation),
        ("BCUT", rdOsmordred.CalcBCUT),
        ("BalabanJ", rdOsmordred.CalcBalabanJ),
        ("BaryszMatrix", rdOsmordred.CalcBaryszMatrix),
        ("BertzCT", rdOsmordred.CalcBertzCT),
        ("BondCount", rdOsmordred.CalcBondCount),
        ("RNCGRPCG", rdOsmordred.CalcRNCGRPCG),
        ("CarbonTypes", lambda mol: rdOsmordred.CalcCarbonTypes(mol)),
        ("Chi", rdOsmordred.CalcChi),
        ("Constitutional", rdOsmordred.CalcConstitutional),
        ("DetourMatrix", rdOsmordred.CalcDetourMatrix),
        ("DistanceMatrix", lambda mol: rdOsmordred.CalcDistanceMatrix(mol)),
        ("EState", lambda mol: rdOsmordred.CalcEState(mol, doExEstate)),
        ("EccentricConnectivityIndex", rdOsmordred.CalcEccentricConnectivityIndex),
        ("ExtendedTopochemicalAtom", rdOsmordred.CalcExtendedTopochemicalAtom),
        ("FragmentComplexity", rdOsmordred.CalcFragmentComplexity),
        ("Framework", rdOsmordred.CalcFramework),
        ("HydrogenBond", rdOsmordred.CalcHydrogenBond),
    ]

    # Always include LogS and InformationContent with maxradius=5
    descriptors.append(("LogS", rdOsmordred.CalcLogS))
    descriptors.append(("InformationContentv2", lambda mol: rdOsmordred.CalcInformationContent(mol, 5)))

    additional_descriptors = [
        ("KappaShapeIndex", rdOsmordred.CalcKappaShapeIndex),
        ("Lipinski", rdOsmordred.CalcLipinski),
        ("McGowanVolume", rdOsmordred.CalcMcGowanVolume),
        ("MoeType", rdOsmordred.CalcMoeType),
        ("MolecularDistanceEdge", rdOsmordred.CalcMolecularDistanceEdge),
        ("MolecularId", rdOsmordred.CalcMolecularId),
        ("PathCount", rdOsmordred.CalcPathCount),
        ("Polarizability", rdOsmordred.CalcPolarizability),
        ("RingCount", rdOsmordred.CalcRingCount),
        ("RotatableBond", rdOsmordred.CalcRotatableBond),
        ("SLogP", rdOsmordred.CalcSLogP),
        ("TopoPSA", rdOsmordred.CalcTopoPSA),
        ("TopologicalCharge", rdOsmordred.CalcTopologicalCharge),
        ("TopologicalIndex", rdOsmordred.CalcTopologicalIndex),
        ("VdwVolumeABC", rdOsmordred.CalcVdwVolumeABC),
        ("VertexAdjacencyInformation", rdOsmordred.CalcVertexAdjacencyInformation),
        ("WalkCount", rdOsmordred.CalcWalkCount),
        ("Weight", rdOsmordred.CalcWeight),
        ("WienerIndex", rdOsmordred.CalcWienerIndex),
        ("ZagrebIndex", rdOsmordred.CalcZagrebIndex),
    ]

    descriptors.extend(additional_descriptors)
    
    # Always include extended descriptors
    extended_descriptors = [
        ("Pol", rdOsmordred.CalcPol),
        ("MR", rdOsmordred.CalcMR),
        ("Flexibility", rdOsmordred.CalcFlexibility),
        ("Schultz", rdOsmordred.CalcSchultz),
        ("AlphaKappaShapeIndex", rdOsmordred.CalcAlphaKappaShapeIndex),
        ("HEState", rdOsmordred.CalcHEState),
        ("BEState", rdOsmordred.CalcBEState),
        ("Abrahams", rdOsmordred.CalcAbrahams),
        ("ANMat", rdOsmordred.CalcANMat),
        ("ASMat", rdOsmordred.CalcASMat),
        ("AZMat", rdOsmordred.CalcAZMat),
        ("DSMat", rdOsmordred.CalcDSMat),
        ("DN2Mat", rdOsmordred.CalcDN2Mat),
        ("Frags", rdOsmordred.CalcFrags),
        ("AddFeatures", rdOsmordred.CalcAddFeatures),
    ]
    descriptors.extend(extended_descriptors)

    for name, func in descriptors:
        try:
            value = func(mol)
            value = np.atleast_1d(np.array(value))
            results.append(value)
        except Exception as e:
            # Handle errors gracefully
            if mynames:
                arraylength = sum(1 for c in mynames if c.startswith(name))
            else:
                # Estimate array length based on function type
                arraylength = 1  # Default to 1 if we can't determine
            warnings.warn(f"Error computing {name} for SMILES: {smiles}. Error: {str(e)}")
            results.append(np.full((arraylength,), np.nan))
        
        if names:
            descriptor_names.extend([f"{name}_{i+1}" for i in range(len(results[-1]))])

    if names:
        return np.concatenate(results), descriptor_names
    return np.concatenate(results)


def Calculate(smiles_list: List[str], ids: Optional[List] = None, n_jobs: int = 4,
     names: bool = False, mynames: Optional[List[str]] = None) -> pd.DataFrame:
    """
    Compute molecular descriptors for a list of SMILES while keeping track of IDs.

    Args:
        smiles_list (List[str]): List of SMILES strings.
        ids (List, optional): List of unique identifiers (same length as smiles_list).
        n_jobs (int, optional): Number of parallel processes (default: 4).
        names (bool, optional): Whether to include names.
        mynames (List[str], optional): Custom names list.

    Returns:
        pd.DataFrame: DataFrame with results, indexed by ID with SMILES as first column.
    """
    if ids is None:
        ids = list(range(len(smiles_list)))

    if len(smiles_list) != len(ids):
        raise ValueError("smiles_list and ids must have the same length")

    results = []
    
    with ProcessPoolExecutor(max_workers=n_jobs) as executor:
        # Submit tasks with their indices
        futures = {
            executor.submit(CalcOsmordred, smi, names=False, mynames=mynames): (smi, mol_id)
            for smi, mol_id in zip(smiles_list, ids)
        }
        
        for future in tqdm(as_completed(futures), total=len(futures), desc="Processing molecules"):
            smi, mol_id = futures[future]
            try:
                result = future.result()
                if result is not None:
                    results.append((smi, mol_id, result))
                else:
                    warnings.warn(f"Failed to process molecule {mol_id} with SMILES: {smi}")
            except Exception as e:
                warnings.warn(f"Error processing molecule {mol_id} with SMILES {smi}: {e}")

    # Sort results by the original ID order
    results.sort(key=lambda x: ids.index(x[1]))

    # Convert to DataFrame
    df_results = pd.DataFrame([res[2] for res in results], 
                             index=pd.Index([res[1] for res in results], name="ID"))
    
    if mynames and len(mynames) > 0:
        df_results.columns = mynames
    
    df_results.insert(0, "cansmi", [res[0] for res in results])

    return df_results


def GetDescriptorNames() -> List[str]:
    """
    Get the list of descriptor names for a given version.
    
    Args:
        version (int): Version of descriptors (1 or 2, default: 2)
        
    Returns:
        List[str]: List of descriptor names
    """
    # version is ignored; always return v2 set
    _, names = CalcOsmordred('CCO', names=True)
    return names


# Convenience function for single molecule processing
def CalcOsmordredFromMol(mol,  names: bool = False,
                         mynames: Optional[List[str]] = None) -> Union[np.ndarray, Tuple[np.ndarray, List[str]]]:
    """
    Calculate Osmordred descriptors from an RDKit molecule object.
    
    Args:
        mol: RDKit molecule object
        version (int): Version of descriptors to use (1 or 2, default: 2)
        names (bool): Whether to return descriptor names along with values
        mynames (List[str], optional): Custom descriptor names list
        
    Returns:
        Union[np.ndarray, Tuple[np.ndarray, List[str]]]: 
            - If names=False: numpy array of descriptor values
            - If names=True: tuple of (descriptor_values, descriptor_names)
    """
    from rdkit import Chem
    smiles = Chem.MolToSmiles(mol)
    return CalcOsmordred(smiles, names, mynames)
