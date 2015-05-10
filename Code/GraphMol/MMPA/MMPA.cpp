//
//  Copyright (C) 2015 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <map>
#include <vector>
#include <algorithm>
#include <math.h>
#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>
//#include <iostream>
//#include <sstream>
//#include <boost/property_tree/ptree.hpp>

#include "../SmilesParse/SmilesParse.h"
#include "../SmilesParse/SmilesWrite.h"
#include "../Substruct/SubstructMatch.h"

//#include "../SmilesParse/SmartsWrite.h"

#include "MMPA.h"

namespace RDKit {
  namespace MMPA {

    typedef std::vector< std::pair<unsigned, unsigned> >  BondVector_t; //pair of BeginAtomIdx, EndAtomIdx

    static inline
    void convertMatchingToBondVect(std::vector<BondVector_t >& matching_bonds,
                             const std::vector<MatchVectType>& matching_atoms, const ROMol& mol){
        for(std::vector<MatchVectType>::const_iterator m=matching_atoms.begin(); m!=matching_atoms.end(); ++m){
            matching_bonds.push_back(BondVector_t());
            BondVector_t& mb = matching_bonds.back();    //current match
            // assume patern is only one bond pattern
            unsigned a1 = (unsigned) (*m)[0].second; // mol atom 1 index
            unsigned a2 = (unsigned) (*m)[1].second; // mol atom 2 index
            mb.push_back(std::pair<unsigned, unsigned>(a1, a2));
//            // loop throght all matched bonds
//            for(MatchVectType::const_iterator mai=m->begin(); mai!=m->end(); ++mai){
//            }
        }
    }

    static inline
    void appendBonds(BondVector_t& bonds,
               const BondVector_t& matching_bonds){
        for(BondVector_t::const_iterator b=matching_bonds.begin(); b!=matching_bonds.end(); ++b)
            bonds.push_back(*b);
    }

    static
    void addResult(std::vector< std::pair<ROMOL_SPTR,ROMOL_SPTR> >& res,
                   const ROMol& mol, const BondVector_t& bonds_selected) {
        RWMol em(mol);
        // loop through the bonds to delete. == deleteBonds()
        unsigned isotope = 0;
        std::map<unsigned, unsigned> isotope_track;
        for(size_t i=0; i < bonds_selected.size(); i++ ){
            isotope += 1;
            // remove the bond
            em.removeBond(bonds_selected[i].first, bonds_selected[i].second);

            // now add attachement points
            unsigned newAtomA = em.addAtom(new Atom(0), true, true);
            em.addBond(bonds_selected[i].first, newAtomA, Bond::SINGLE);
            unsigned newAtomB = em.addAtom(new Atom(0), true, true);
            em.addBond(bonds_selected[i].second, newAtomB, Bond::SINGLE);

            // keep track of where to put isotopes
            isotope_track[newAtomA] = isotope;
            isotope_track[newAtomB] = isotope;
        }
            
        // should be able to get away without sanitising mol as the existing valencies/atoms not changed
        ROMol modifiedMol(em);
        // Split ions:
        //#canonical smiles can be different with and without the isotopes
        //#hence to keep track of duplicates use fragmented_smi_noIsotopes
        std::string fragmented_smi_noIsotopes = MolToSmiles(modifiedMol, true);

        bool valid = true;
        std::vector<std::string> fragments;

        //#check if it's a valid triple cut
        if(isotope == 3){
            valid = false;
            boost::split(fragments, fragmented_smi_noIsotopes, boost::is_any_of("."), boost::token_compress_on);
            for(std::vector<std::string>::const_iterator f=fragments.begin(); f!=fragments.end(); f++){
                if( boost::regex_search(*f, boost::regex("\\*.*\\*.*\\*"), boost::match_any)){
                    valid = true;
                    break;
                }
            }
        }

        RWMol *m1=NULL, *m2=NULL;   // core & side_chains molecules
        if(!valid)
            return;
// DEBUG PRINT
std::cout<<"isotope="<< isotope <<"\n";
        if(isotope == 1){
            m2 = new RWMol(em); // output = '%s,%s,,%s.%s'
// * DEBUG PRINT
//? replace any [*] with [:1]
//?            fragmented_smi_noIsotopes = re.sub('\[\*\]', '[*:1]', fragmented_smi_noIsotopes)
            boost::split(fragments, fragmented_smi_noIsotopes, boost::is_any_of("."), boost::token_compress_on);
            RWMol* s1 = SmilesToMol(fragments[0]);
            RWMol* s2 = SmilesToMol(fragments[1]);

            //#need to cansmi again as smiles can be different
            std::cout<<"smi,id,,"<< MolToSmiles(*s1, true)<<"."<<MolToSmiles(*s2, true)<<"\n";
// * /
        }
        else if(isotope >= 2){
/*
            #add the isotope labels
            for key in isotope_track:
                #to add isotope lables
                modifiedMol.GetAtomWithIdx(key).SetIsotope(isotope_track[key])
            fragmented_smi = Chem.MolToSmiles(modifiedMol,isomericSmiles=True)

            //#change the isotopes into labels - currently can't add SMARTS or labels to mol
            fragmented_smi = re.sub('\[1\*\]', '[*:1]', fragmented_smi)
            fragmented_smi = re.sub('\[2\*\]', '[*:2]', fragmented_smi)
            fragmented_smi = re.sub('\[3\*\]', '[*:3]', fragmented_smi)

            fragments = fragmented_smi.split(".")

            #identify core/side chains and cansmi them
            core,side_chains = find_correct(fragments)

            #now change the labels on sidechains and core
            #to get the new labels, cansmi the dot-disconnected side chains
            #the first fragment in the side chains has attachment label 1, 2nd: 2, 3rd: 3
            #then change the labels accordingly in the core

            #this is required by the indexing script, as the side-chains are "keys" in the index
            #this ensures the side-chains always have the same numbering

            isotope_track = {}
            side_chain_fragments = side_chains.split(".")

            for s in xrange( len(side_chain_fragments) ):
                matchObj = re.search( '\[\*\:([123])\]', side_chain_fragments[s] )
                if matchObj:
                    #add to isotope_track with key: old_isotope, value:
                    isotope_track[matchObj.group(1)] = str(s+1)

            #change the labels if required
            if(isotope_track['1'] != '1'):
                core = re.sub('\[\*\:1\]', '[*:XX' + isotope_track['1'] + 'XX]' , core)
                side_chains = re.sub('\[\*\:1\]', '[*:XX' + isotope_track['1'] + 'XX]' , side_chains)
            if(isotope_track['2'] != '2'):
                core = re.sub('\[\*\:2\]', '[*:XX' + isotope_track['2'] + 'XX]' , core)
                side_chains = re.sub('\[\*\:2\]', '[*:XX' + isotope_track['2'] + 'XX]' , side_chains)

            if(isotope == 3):
                if(isotope_track['3'] != '3'):
                    core = re.sub('\[\*\:3\]', '[*:XX' + isotope_track['3'] + 'XX]' , core)
                    side_chains = re.sub('\[\*\:3\]', '[*:XX' + isotope_track['3'] + 'XX]' , side_chains)

            #now remove the XX
            core = re.sub('XX', '' , core)
            side_chains = re.sub('XX', '' , side_chains)

            output = '%s,%s,%s,%s' % (smi,id,core,side_chains)
            if( (output in out) == False):
                out.add(output)
*/
        }

        res.push_back(std::pair<ROMOL_SPTR,ROMOL_SPTR>(ROMOL_SPTR(m1), ROMOL_SPTR(m2)));
    }

    static
    void processCuts(size_t i, size_t maxCuts, BondVector_t& bonds_selected, const std::vector<BondVector_t>& matching_bonds,
                     const ROMol& mol, std::vector< std::pair<ROMOL_SPTR,ROMOL_SPTR> >& res){
        for(size_t x=i; x < matching_bonds.size(); x++ ){
            appendBonds(bonds_selected, matching_bonds[x]);
            addResult(res, mol, bonds_selected);
            for(i++; i < maxCuts; i++ ){
                processCuts (i, maxCuts, bonds_selected, matching_bonds, mol, res);
            }
            bonds_selected.pop_back();
        }
    }

    bool fragmentMol(const ROMol& mol,
                     std::vector< std::pair<ROMOL_SPTR,ROMOL_SPTR> >& res,
                     unsigned int maxCuts,
                     const std::string& pattern) {
        res.clear();
        std::auto_ptr<const ROMol> smarts((const ROMol*)SmartsToMol(pattern));
        std::vector<MatchVectType> matching_atoms; //one bond per match ! with default pattern
        unsigned int total = SubstructMatch(mol, *smarts, matching_atoms);
        if(0==total){
//tmp            res.push_back(std::pair<ROMOL_SPTR,ROMOL_SPTR>(NULL, NULL)); //print "mol,id,,"
            return false;
        }

        std::vector<BondVector_t> matching_bonds; // List of matched query's bonds
        convertMatchingToBondVect(matching_bonds, matching_atoms, mol);

        //loop to generate every single, double and triple cut in the molecule
        BondVector_t bonds_selected;
        processCuts (0, maxCuts, bonds_selected, matching_bonds, mol, res);
        return true;
    }
  }
}   // namespace RDKit

