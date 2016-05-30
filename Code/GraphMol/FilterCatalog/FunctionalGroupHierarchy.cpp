//
//  Copyright (c) 2016, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include "FilterMatchers.h"
#include "FilterCatalog.h"
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/shared_ptr.hpp>
#ifdef RDK_THREADSAFE_SSS
#include <boost/thread/once.hpp>
#endif
#include <boost/algorithm/string.hpp>
#include <RDGeneral/BoostEndInclude.h>

namespace RDKit {

namespace {
struct FuncData_t {
  int level;
  const char *name;
  const char *smarts;
  const char *label;
  const char *removalReaction;
};

// Exploit the fact that the current hierarchy is a three level tree
//  0 -> gets added to root node
//  1 -> gets added to last 0 node
//  2 -> gets added to last 1 node
const int MAX_DEPTH = 3;

const FuncData_t FuncDataArray[] = {
    {0, "AcidChloride", "C(=O)Cl", "Acid Chloride", 0},
    {1, "AcidChloride.Aromatic", "[$(C-!@[a])](=O)(Cl)", "Aromatic", 0},
    {1, "AcidChloride.Aliphatic", "[$(C-!@[A;!Cl])](=O)(Cl)", "Aliphatic", 0},

    {0, "CarboxylicAcid", "C(=O)[O;H,-]", "Carboxylic acid", 0},
    {1, "CarboxylicAcid.Aromatic", "[$(C-!@[a])](=O)([O;H,-])", "Aromatic", 0},
    {1, "CarboxylicAcid.Aliphatic", "[$(C-!@[A;!O])](=O)([O;H,-])", "Aliphatic",
     0},
    {1, "CarboxylicAcid.AlphaAmino",
     "[$(C-[C;!$(C=[!#6])]-[N;!H0;!$(N-[!#6;!#1]);!$(N-C=[O,N,S])])](=O)([O;H,-"
     "])",
     "alpha Amino Acid", 0},

    {0, "SulfonylChloride", "[$(S-!@[#6])](=O)(=O)(Cl)", "Sulfonyl Chloride",
     0},
    {1, "SulfonylChloride.Aromatic", "[$(S-!@c)](=O)(=O)(Cl)", "Aromatic", 0},
    {1, "SulfonylChloride.Aliphatic", "[$(S-!@C)](=O)(=O)(Cl)", "Aliphatic", 0},

    {0, "Amine", "[N;$(N-[#6]);!$(N-[!#6;!#1]);!$(N-C=[O,N,S])]", "Amine", 0},
    {1, "Amine.Primary", "[N;H2;D1;$(N-!@[#6]);!$(N-C=[O,N,S])]", "Primary", 0},
    {2, "Amine.Primary.Aromatic", "[N;H2;D1;$(N-!@c);!$(N-C=[O,N,S])]",
     "Primary aromatic", 0},
    {2, "Amine.Primary.Aliphatic", "[N;H2;D1;$(N-!@C);!$(N-C=[O,N,S])]",
     "Primary aliphatic", 0},
    {1, "Amine.Secondary", "[N;H1;D2;$(N(-[#6])-[#6]);!$(N-C=[O,N,S])]",
     "Secondary", 0},
    {2, "Amine.Secondary.Aromatic", "[N;H1;D2;$(N(-[c])-[#6]);!$(N-C=[O,N,S])]",
     "Secondary aromatic", 0},
    {2, "Amine.Secondary.Aliphatic", "[N;H1;D2;$(N(-C)-C);!$(N-C=[O,N,S])]",
     "Secondary aliphatic", 0},
    {1, "Amine.Tertiary", "[N;H0;D3;$(N(-[#6])(-[#6])-[#6]);!$(N-C=[O,N,S])]",
     "Tertiary", 0},
    {2, "Amine.Tertiary.Aromatic",
     "[N;H0;D3;$(N(-[c])(-[#6])-[#6]);$(N-C=[O,N,S])]", "Tertiary aromatic", 0},
    {2, "Amine.Tertiary.Aliphatic", "[N;H0;D3;$(N(-C)(-C)-C);!$(N-C=[O,N,S])]",
     "Tertiary aliphatic", 0},
    {1, "Amine.Aromatic", "[N;$(N-c);!$(N-[!#6;!#1]);!$(N-C=[O,N,S])]",
     "Aromatic", 0},
    {1, "Amine.Aliphatic", "[N;!$(N-c);$(N-C);!$(N-[!#6;!#1]);!$(N-C=[O,N,S])]",
     "Aliphatic", 0},
    {1, "Amine.Cyclic", "[N;R;$(N-[#6]);!$(N-[!#6;!#1]);!$(N-C=[O,N,S])]",
     "Cyclic", 0},

    {0, "BoronicAcid", "[$(B-!@[#6])](O)(O)", "Boronic Acid",
     "[#6:1]-!@[B:2]([O:3])[O:4]>>[#6:1][X].[B:2]([O:3])[O:4]"},
    {1, "BoronicAcid.Aromatic", "[$(B-!@c)](O)(O)", "Aromatic",
     "[c:1]-!@[B:2]([O:3])[O:4]>>[c:1][X].[B:2]([O:3])[O:4]"},
    {1, "BoronicAcid.Aliphatic", "[$(B-!@C)](O)(O)", "Aliphatic",
     "[C:1]-!@[B:2]([O:3])[O:4]>>[C:1][X].[B:2]([O:3])[O:4]"},

    {0, "Isocyanate", "[$(N-!@[#6])](=!@C=!@O)", "Isocyanate", 0},
    {1, "Isocyanate.Aromatic", "[$(N-!@c)](=!@C=!@O)", "Aromatic", 0},
    {1, "Isocyanate.Aliphatic", "[$(N-!@C)](=!@C=!@O)", "Aliphatic", 0},

    {0, "Alcohol", "[O;H1;$(O-!@[#6;!$(C=!@[O,N,S])])]", "Alcohol", 0},
    {1, "Alcohol.Aromatic", "[O;H1;$(O-!@c)]", "Aromatic", 0},
    {1, "Alcohol.Aliphatic", "[O;H1;$(O-!@[C;!$(C=!@[O,N,S])])]", "Aliphatic",
     0},

    {0, "Aldehyde", "[CH;D2;!$(C-[!#6;!#1])]=O", "Aldehyde", 0},
    {1, "Aldehyde.Aromatic", "[CH;D2;$(C-!@[a])](=O)", "Aromatic", 0},
    {1, "Aldehyde.Aliphatic", "[CH;D2;$(C-!@C)](=O)", "Aliphatic", 0},

    {0, "Halogen",
     "[$([F,Cl,Br,I]-!@[#6]);!$([F,Cl,Br,I]-!@C-!@[F,Cl,Br,I]);!$([F,Cl,Br,I]-["
     "C,S](=[O,S,N]))]",
     "Halogen", 0},
    {1, "Halogen.Aromatic", "[F,Cl,Br,I;$(*-!@c)]", "Aromatic", 0},
    {1, "Halogen.Aliphatic",
     "[$([F,Cl,Br,I]-!@C);!$([F,Cl,Br,I]-!@C-!@[F,Cl,Br,I])]", "Aliphatic", 0},
    {1, "Halogen.NotFluorine",
     "[$([Cl,Br,I]-!@[#6]);!$([Cl,Br,I]-!@C-!@[F,Cl,Br,I]);!$([Cl,Br,I]-[C,S](="
     "[O,S,N]))]",
     "Not Fluorine", 0},
    {2, "Halogen.NotFluorine.Aliphatic",
     "[$([Cl,Br,I]-!@C);!$([Cl,Br,I]-!@C-!@[F,Cl,Br,I]);!$([Cl,Br,I]-[C,S](=[O,"
     "S,N]))]",
     "Aliphatic Not Fluorine", 0},
    {2, "Halogen.NotFluorine.Aromatic", "[$([Cl,Br,I]-!@c)]",
     "Aromatic Not Fluorine", 0},
    {1, "Halogen.Bromine",
     "[Br;$([Br]-!@[#6]);!$([Br]-!@C-!@[F,Cl,Br,I]);!$([Br]-[C,S](=[O,S,N]))]",
     "Bromine", 0},
    {2, "Halogen.Bromine.Aliphatic",
     "[Br;$(Br-!@C);!$(Br-!@C-!@[F,Cl,Br,I]);!$(Br-[C,S](=[O,S,N]))]",
     "Aliphatic Bromine", 0},
    {2, "Halogen.Bromine.Aromatic", "[Br;$(Br-!@c)]", "Aromatic Bromine", 0},
    {2, "Halogen.Bromine.BromoKetone", "[Br;$(Br-[CH2]-C(=O)-[#6])]",
     "Bromoketone", 0},

    {0, "Azide", "[N;H0;$(N-[#6]);D2]=[N;D2]=[N;D1]", "Azide", 0},
    {1, "Azide.Aromatic", "[N;H0;$(N-c);D2]=[N;D2]=[N;D1]", "Aromatic Azide",
     0},
    {1, "Azide.Aliphatic", "[N;H0;$(N-C);D2]=[N;D2]=[N;D1]", "Aliphatic Azide",
     0},

    {0, "Nitro", "[N;H0;$(N-[#6]);D3](=[O;D1])~[O;D1]", "Nitro", 0},
    {1, "Nitro.Aromatic", "[N;H0;$(N-c);D3](=[O;D1])~[O;D1]", "Aromatic Nitro",
     0},
    {1, "Nitro.Aliphatic", "[N;H0;$(N-C);D3](=[O;D1])~[O;D1]",
     "Aliphatic Nitro", 0},

    {0, "TerminalAlkyne", "[C;$(C#[CH])]", "Terminal Alkyne", 0}};

const unsigned int NUM_FUNCS =
    static_cast<unsigned int>(sizeof(FuncDataArray) / sizeof(FuncData_t));

FilterCatalog &hierarchy_get() {
  static FilterCatalog fgroup;
  return fgroup;
}

std::map<std::string, ROMOL_SPTR> &flatten_get() {
  static std::map<std::string, ROMOL_SPTR> flattenedHierarchy;
  return flattenedHierarchy;
}

std::map<std::string, ROMOL_SPTR> &flatten_normalized_get() {
  static std::map<std::string, ROMOL_SPTR> flattenedHierarchy;
  return flattenedHierarchy;
}

void hierarchy_create() {
  FilterCatalog &fgroupHierarchy = hierarchy_get();
  std::map<std::string, ROMOL_SPTR> &flattenedHierarchy = flatten_get();
  std::map<std::string, ROMOL_SPTR> &flattenedHierarchyNorm =
      flatten_normalized_get();

  std::vector<FilterHierarchyMatcher *> toplevel;
  FilterHierarchyMatcher *stack[MAX_DEPTH];

  for (size_t i = 0; i < NUM_FUNCS; ++i) {
    // Make a new node
    ROMOL_SPTR pattern(SmartsToMol(FuncDataArray[i].smarts, 0, true));
    pattern->setProp("Label", FuncDataArray[i].label);
    if (FuncDataArray[i].removalReaction) {
      pattern->setProp("RemovalReaction", FuncDataArray[i].removalReaction);
    }
    std::string key(FuncDataArray[i].name);
    flattenedHierarchy[key] = pattern;
    boost::to_lower(key);
    flattenedHierarchyNorm[key] = pattern;

    FilterHierarchyMatcher node(SmartsMatcher(FuncDataArray[i].name, pattern));

    if (FuncDataArray[i].level == 0) {
      toplevel.push_back(new FilterHierarchyMatcher(node));
      stack[0] = toplevel.back();

    } else {
      PRECONDITION(
          FuncDataArray[i].level < MAX_DEPTH,
          std::string(
              "Invalid Depth in Built in Functional Group Hierarchy: ") +
              FuncDataArray[i].name);
      boost::shared_ptr<FilterHierarchyMatcher> real_node =
          stack[FuncDataArray[i].level - 1]->addChild(node);

      stack[FuncDataArray[i].level] = real_node.get();
    }
  }

  // add the top levels to the filter catalog and give them ownership
  //  to the filter catalog
  for (size_t i = 0; i < toplevel.size(); ++i) {
    fgroupHierarchy.addEntry(new FilterCatalogEntry(
        toplevel[i]->getName(),
        boost::shared_ptr<FilterMatcherBase>(toplevel[i])));
  }
}
}

const FilterCatalog &GetFunctionalGroupHierarchy() {
#ifdef RDK_THREADSAFE_SSS
#ifdef BOOST_THREAD_PROVIDES_ONCE_CXX11
  static boost::once_flag flag;
#else
  static boost::once_flag flag = BOOST_ONCE_INIT;
#endif
  boost::call_once(&hierarchy_create, flag);
#else
  static bool loaded = false;
  if (!loaded) {
    hierarchy_create();
    loaded = true;
  }
#endif
  return hierarchy_get();
}

const std::map<std::string, ROMOL_SPTR> &GetFlattenedFunctionalGroupHierarchy(
    bool normalize) {
  GetFunctionalGroupHierarchy();
  if (normalize) return flatten_normalized_get();
  return flatten_get();
}
}
