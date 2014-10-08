//
//  Copyright (c) 2007-2014, Novartis Institutes for BioMedical Research Inc.
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
//       products derived from this software without specific prior written permission.
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

#ifndef __RD_REACTION_H_17Aug2006__
#define __RD_REACTION_H_17Aug2006__

#include <GraphMol/RDKitBase.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <vector>

namespace RDKit{
  class ReactionPickler;
     
  //! used to indicate an error in the chemical reaction engine
  class ChemicalReactionException : public std::exception {
  public:
    //! construct with an error message
    explicit ChemicalReactionException(const char *msg) : _msg(msg) {};
    //! construct with an error message
    explicit ChemicalReactionException(const std::string msg) : _msg(msg) {};
    //! get the error message
    const char *message () const { return _msg.c_str(); };
    ~ChemicalReactionException () throw () {};
  private:
    std::string _msg;
  };
   
  //! This is a class for storing and applying general chemical reactions.
  /*!
     basic usage will be something like:
     
     \verbatim
     ChemicalReaction rxn;
     rxn.addReactantTemplate(r1);
     rxn.addReactantTemplate(r2);
     rxn.addProductTemplate(p1);
     rxn.initReactantMatchers();
     
     MOL_SPTR_VECT prods;
     for(MOL_SPTR_VECT::const_iterator r1It=reactantSet1.begin();
         r1It!=reactantSet1.end();++r1It;){
       for(MOL_SPTR_VECT::const_iterator r2It=reactantSet2.begin();
	   r2It!=reactantSet2.end();++r2It;){
	 MOL_SPTR_VECT rVect(2);
         rVect[0] = *r1It;
         rVect[1] = *r2It;
             
         std::vector<MOL_SPTR_VECT> lprods;
         lprods = rxn.runReactants(rVect);
         for(std::vector<MOL_SPTR_VECT>::const_iterator lpIt=lprods.begin();
            lpIt!=lprods.end();++lpIt){
            // we know this is a single-product reaction:
            prods.push_back((*lpIt)[0]);
         }
       }     
     }
     \endverbatim     

     NOTES:
       - to allow more control over the reaction, it is possible to flag reactant
         atoms as being protected by setting the "_protected" property on those
         atoms. Here's an example:
         \verbatim
            std::string smi="[O:1]>>[N:1]";
            ChemicalReaction *rxn = RxnSmartsToChemicalReaction(smi); 
            rxn->initReactantMatchers();
        
            MOL_SPTR_VECT reacts;
            reacts.clear();
            smi = "OCO";
            ROMol *mol = SmilesToMol(smi);
            reacts.push_back(ROMOL_SPTR(mol));
            std::vector<MOL_SPTR_VECT> prods;
            prods = rxn->runReactants(reacts);
            // here prods has two entries, because there are two Os in the
            // reactant. 

            reacts[0]->getAtomWithIdx(0)->setProp("_protected",1);
            prods = rxn->runReactants(reacts);
            // here prods only has one entry, the reaction at atom 0
            // has been blocked by the _protected property
         \endverbatim

  */
  class ChemicalReaction {
    friend class ReactionPickler;

  public:
    ChemicalReaction() : df_needsInit(true), df_implicitProperties(false) {};
    ChemicalReaction(const ChemicalReaction &other){
        df_needsInit=other.df_needsInit;
        df_implicitProperties=other.df_implicitProperties;
        m_reactantTemplates=other.m_reactantTemplates;
        m_productTemplates=other.m_productTemplates;
        m_agentTemplates=other.m_agentTemplates;
    }
    //! construct a reaction from a pickle string
    ChemicalReaction(const std::string &binStr);

    //! Adds a new reactant template
    /*!
      \return the number of reactants

    */
    unsigned int addReactantTemplate(ROMOL_SPTR mol){
      this->df_needsInit = true;
      this->m_reactantTemplates.push_back(mol);
      return this->m_reactantTemplates.size();
    }

    //! Adds a new agent template
    /*!
      \return the number of agent

    */
    unsigned int addAgentTemplate(ROMOL_SPTR mol){
      this->m_agentTemplates.push_back(mol);
      return this->m_agentTemplates.size();
    }

    //! Adds a new product template
    /*!
      \return the number of products

    */
    unsigned int addProductTemplate(ROMOL_SPTR mol){
      this->m_productTemplates.push_back(mol);
      return this->m_productTemplates.size();
    }
    
    //! Removes the reactant templates from a reaction if atom mapping ratio is below a given threshold
    /*! By default the removed reactant templates were attached to the agent templates.
        An alternative will be to provide a pointer to a molecule vector where these reactants should be saved.
    */
    void removeUnmappedReactantTemplates(
      double thresholdUnmappedAtoms = 0.2,
      bool moveToAgentTemplates = true,
      MOL_SPTR_VECT *targetVector = NULL);

    //! Removes the product templates from a reaction if its atom mapping ratio is below a given threshold
    /*! By default the removed products templates were attached to the agent templates.
        An alternative will be to provide a pointer to a molecule vector where these products should be saved.
    */
    void removeUnmappedProductTemplates(
      double thresholdUnmappedAtoms =0.2,
      bool moveToAgentTemplates = true,
      MOL_SPTR_VECT *targetVector = NULL);
      
    //! Runs the reaction on a set of reactants
    /*!
     
      \param reactants: the reactants to be used. The length of this must be equal to
                        this->getNumReactantTemplates()
                         
      \return a vector of vectors of products. Each subvector will be 
              this->getNumProductTemplates() long.          
      
      We return a vector of vectors of products because each individual template may 
      map multiple times onto its reactant. This leads to multiple possible result
      sets.
    */
    std::vector<MOL_SPTR_VECT> runReactants(const MOL_SPTR_VECT reactants) const;

    MOL_SPTR_VECT::const_iterator beginReactantTemplates() const {
        return this->m_reactantTemplates.begin();    
    }
    MOL_SPTR_VECT::const_iterator endReactantTemplates() const {
        return this->m_reactantTemplates.end();    
    }

    MOL_SPTR_VECT::const_iterator beginProductTemplates() const {
        return this->m_productTemplates.begin();    
    }
    MOL_SPTR_VECT::const_iterator endProductTemplates() const {
        return this->m_productTemplates.end();    
    }

    MOL_SPTR_VECT::const_iterator beginAgentTemplates() const {
        return this->m_agentTemplates.begin();
    }
    MOL_SPTR_VECT::const_iterator endAgentTemplates() const {
        return this->m_agentTemplates.end();
    }

    MOL_SPTR_VECT::iterator beginReactantTemplates() {
        return this->m_reactantTemplates.begin();    
    }
    MOL_SPTR_VECT::iterator endReactantTemplates() {
        return this->m_reactantTemplates.end();    
    }

    MOL_SPTR_VECT::iterator beginProductTemplates() {
        return this->m_productTemplates.begin();    
    }
    MOL_SPTR_VECT::iterator endProductTemplates() {
        return this->m_productTemplates.end();    
    }

    MOL_SPTR_VECT::iterator beginAgentTemplates() {
        return this->m_agentTemplates.begin();
    }
    MOL_SPTR_VECT::iterator endAgentTemplates() {
        return this->m_agentTemplates.end();
    }
    unsigned int getNumReactantTemplates() const { return this->m_reactantTemplates.size(); };
    unsigned int getNumProductTemplates() const { return this->m_productTemplates.size(); };
    unsigned int getNumAgentTemplates() const { return this->m_agentTemplates.size(); };

    //! initializes our internal reactant-matching datastructures.
    /*! 
        This must be called after adding reactants and before calling
        runReactants.
    */
    void initReactantMatchers();

    bool isInitialized() const { return !df_needsInit; };
    
    //! validates the reactants and products to make sure the reaction seems "reasonable"
    /*! 
        \return   true if the reaction validates without errors (warnings do not stop
                  validation)
         
        \param numWarnings: used to return the number of validation warnings
        \param numErrors:   used to return the number of validation errors
        
        \param silent: If this bool is true, no messages will be logged during the validation. 
                       By default, validation problems are reported to the warning and error 
                       logs depending on their severity.
                       
    */
    bool validate(unsigned int &numWarnings,unsigned int &numErrors,bool silent=false) const;
        

    //! returns whether or not the reaction uses implicit
    //! properties on the product atoms
    /*!

      This toggles whether or not unspecified atomic properties in the
      products are considered to be implicit and should be copied from
      the actual reactants. This is necessary due to a semantic difference
      between the "reaction SMARTS" approach and the MDL RXN
      approach:
        In "reaction SMARTS", this reaction:
          [C:1]-[Br:2].[O-:3]>>[C:1]-[O:3].[Br-:2]
        applied to [CH4+]Br should yield [CH4+]O
        Something similar drawn in an rxn file, and applied to
        [CH4+]Br should yield [CH3]O. 
        In rxn there is no charge on the product C because nothing is
        specified in the rxn file; in "SMARTS" the charge from the
        actual reactants is not *removed* because no charge is
        specified in the reaction.

    */
    bool getImplicitPropertiesFlag() const { return df_implicitProperties; };
    //! sets the implicit properties flag. See the documentation for
    //! getImplicitProertiesFlag() for a discussion of what this means.
    void setImplicitPropertiesFlag(bool val) { df_implicitProperties=val; };

  private:
    bool df_needsInit;
    bool df_implicitProperties;
    MOL_SPTR_VECT m_reactantTemplates,m_productTemplates,m_agentTemplates;
    ChemicalReaction &operator=(const ChemicalReaction &); // disable assignment
  };

  //! tests whether or not the molecule has a substructure match
  //! to any of the reaction's reactants
  //! the \c which argument is used to return which of the reactants
  //! the molecule matches. If there's no match, it is equal to the number
  //! of reactants on return
  bool isMoleculeReactantOfReaction(const ChemicalReaction &rxn,const ROMol &mol,
                                      unsigned int &which);
  //! \overload
  bool isMoleculeReactantOfReaction(const ChemicalReaction &rxn,const ROMol &mol);
  
  //! tests whether or not the molecule has a substructure match
  //! to any of the reaction's products
  //! the \c which argument is used to return which of the products
  //! the molecule matches. If there's no match, it is equal to the number
  //! of products on return
  bool isMoleculeProductOfReaction(const ChemicalReaction &rxn,const ROMol &mol,
                                   unsigned int &which);
  //! \overload
  bool isMoleculeProductOfReaction(const ChemicalReaction &rxn,const ROMol &mol);

  //! tests whether or not the molecule has a substructure match
  //! to any of the reaction's agents
  //! the \c which argument is used to return which of the agents
  //! the molecule matches. If there's no match, it is equal to the number
  //! of agents on return
  bool isMoleculeAgentOfReaction(const ChemicalReaction &rxn,const ROMol &mol,
                                   unsigned int &which);
  //! \overload
  bool isMoleculeAgentOfReaction(const ChemicalReaction &rxn,const ROMol &mol);

  //! returns indices of the atoms in each reactant that are changed
  //! in the reaction
  /*!
    \param rxn the reaction we are interested in

    \param mappedAtomsOnly if set, atoms that are not mapped will not be included in
         the list of changed atoms (otherwise they are automatically included)

     How are changed atoms recognized?
         1) Atoms whose degree changes 
         2) Atoms whose bonding pattern changes
         3) unmapped atoms (unless the mappedAtomsOnly flag is set)
         4) Atoms connected to unmapped atoms
         5) Atoms whose atomic number changes (unless the
            corresponding product atom is a dummy)
         6) Atoms with more than one atomic number query (unless the
            corresponding product atom is a dummy)

     Note that the atomic number of a query atom depends on how it's constructed.
       When coming from SMARTS: if the first query is an atomic label/number that
          sets the atomic number, otherwise it's zero.
          For example [O;$(OC)] is atomic number 8 while [$(OC);O] is atomic
          number 0.
       When coming from RXN: the atomic number of the atom in the rxn file sets
          the value.
   */
  VECT_INT_VECT getReactingAtoms(const ChemicalReaction &rxn,bool mappedAtomsOnly=false);

  //! add the recursive queries to the reactants of a reaction
  /*!
    This does its work using RDKit::addRecursiveQueries()

      \param rxn the reaction we are interested in
      \param queries        - the dictionary of named queries to add
      \param propName       - the atom property to use to get query names
      optional:
      \param reactantLabels - to store pairs of (atom index, query string)
                              per reactant
  
      NOTES:
        - existing query information, if present, will be supplemented (AND logic)
        - non-query atoms will be replaced with query atoms using only the query
          logic
        - query names can be present as comma separated lists, they will then
          be combined using OR logic.
        - throws a KeyErrorException if a particular query name is not present
          in \c queries

   */  
  void addRecursiveQueriesToReaction(ChemicalReaction &rxn,
                                     const std::map<std::string,ROMOL_SPTR> &queries,
		  std::string propName,
		  std::vector<std::vector<std::pair<unsigned int,std::string> > > *reactantLabels=NULL);

} // end of RDKit namespace

namespace RDDepict {
  //! \brief Generate 2D coordinates (a depiction) for a reaction
  /*! 

    \param rxn the reaction we are interested in

    \param spacing the spacing between components of the reaction

    \param updateProps if set, properties such as conjugation and
        hybridization will be calculated for the reactant and product
        templates before generating coordinates. This should result in
        better depictions, but can lead to errors in some cases.

    \param canonOrient canonicalize the orientation so that the long
    axes align with the x-axis etc.

    \param nFlipsPerSample - the number of rotatable bonds that are
    flipped at random for each sample

    \param nSamples - the number of samples

    \param sampleSeed - seed for the random sampling process

    \param permuteDeg4Nodes - try permuting the drawing order of bonds around
          atoms with four neighbors in order to improve the depiction

    for the other parameters see the documentation for compute2DCoords()

  */
  void compute2DCoordsForReaction(RDKit::ChemicalReaction &rxn,
                                  double spacing=2.0,
                                  bool updateProps=true,
                                  bool canonOrient=false,
                                  unsigned int nFlipsPerSample=0,
                                  unsigned int nSamples=0,
                                  int sampleSeed=0,
                                  bool permuteDeg4Nodes=false);

} // end of RDDepict namespace

#endif
