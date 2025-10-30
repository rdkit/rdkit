//
//  Copyright (C) 2025 and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
%{
#include <GraphMol/Subset.h>
%}

%ignore RDKit::copyMolSubset;
%include <GraphMol/Subset.h>

%{
RDKit::ROMol *copyMolSubsetHelper(const RDKit::ROMol& mol,
			   const std::vector<unsigned int> &atoms,
			   const std::vector<unsigned int> &bonds,
			   const RDKit::SubsetOptions &options=RDKit::SubsetOptions()
			   )
{
  return copyMolSubset(mol, atoms, bonds, options).release();
}
 

RDKit::ROMol *copyMolSubsetHelper(const RDKit::ROMol& mol,
		     const std::vector<unsigned int> &atoms,
		     const std::vector<unsigned int> &bonds,
		     RDKit::SubsetInfo &subsetInfo,
		     const RDKit::SubsetOptions &options=RDKit::SubsetOptions()
		     ) {
  return copyMolSubset(mol, atoms, bonds, subsetInfo, options).release();
}

RDKit::ROMol *copyMolSubsetHelper(const RDKit::ROMol& mol,
				  const std::vector<unsigned int>& path,
				  const RDKit::SubsetOptions &options = RDKit::SubsetOptions()) {
  return copyMolSubset(mol, path, options).release();
}


RDKit::ROMol *copyMolSubsetHelper(const RDKit::ROMol& mol,
				  const std::vector<unsigned int>& path,
				  RDKit::SubsetInfo &subsetInfo,
				  const RDKit::SubsetOptions &options = RDKit::SubsetOptions()) {
  return copyMolSubset(mol, path, subsetInfo, options).release();
}
  
%}

%rename("copyMolSubset") copyMolSubsetHelper;

RDKit::ROMol *copyMolSubsetHelper(const RDKit::ROMol& mol,
				  const std::vector<unsigned int> &atoms,
				  const std::vector<unsigned int> &bonds,
				  const RDKit::SubsetOptions &options=RDKit::SubsetOptions()
				  );
RDKit::ROMol *copyMolSubsetHelper(const RDKit::ROMol& mol,
				  const std::vector<unsigned int> &atoms,
				  const std::vector<unsigned int> &bonds,
				  RDKit::SubsetInfo &subsetInfo,
				  const RDKit::SubsetOptions &options=SRDKit::ubsetOptions()
				  );
RDKit::ROMol *copyMolSubsetHelper(const RDKit::ROMol& mol,
				  const std::vector<unsigned int>& path,
				  const RDKit::SubsetOptions &options = SubsetOptions());

RDKit::ROMol *copyMolSubsetHelper(const RDKit::ROMol& mol,
				  const std::vector<unsigned int>& path,
				  RDKit::SubsetInfo &subsetInfo,
				  const RDKit::SubsetOptions &options = RDKit::SubsetOptions());
