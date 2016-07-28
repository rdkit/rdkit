//
//  Copyright (C) 2016 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//#include <>
#include "StructChecker.h"

namespace RDKit {
 namespace StructureCheck {

     void createRingList(const ROMol &mol, std::vector<unsigned> &ring_list) {
     }

     void CombineRings(const ROMol &mol, std::vector<unsigned> &ring_list) {
     }

/*
struct tree_cell
{
int color;
int link;
};

bond_set_node *RingList(unsigned bonds[][2], unsigned nbonds)
* Returns a list of basis rings of the graph defined by
* bonds[0..nbonds-1][0..1]
*
* Generalized to ignore any bonds that contain nodes with number 0.
* This can be used to perform ring analysis only on a selected subset
* of the bonds of the source graph.

     {
         int i;
         unsigned int b;
         int at1, at2;
         int level1, level2;
         int trace, tmp;
         int color, old_color, new_color;
         int tmp_link, new_link;
         bond_set_node *result, *p;
         unsigned natoms;
         struct tree_cell *tree;

         result = (bond_set_node *)NULL;

         for (i = 0, natoms = 0; i<nbonds; i++)         //make space for spanning tree 
         {
             if (natoms < bonds[i][0]) natoms = bonds[i][0];
             if (natoms < bonds[i][1]) natoms = bonds[i][1];
         }
         tree = TypeAlloc(natoms + 1, struct tree_cell);

         for (i = 0; i <= natoms; i++)
         {
             tree[i].color = NO_COLOR; tree[i].link = UNLINKED;
         }

         //Add bonds to spanning tree one at a time. If a bond doesn't link 
         //to an old component, label the tree cells for both atoms with a  
         //new (unique) label. If one of the atoms corresponds to a cell    
         //that already has a label, label the cell corresponding to the    
         //other atom with the same label. If both atoms correspond to      
         //different labels, relabel the cells with the higher label to the 
         //lower one to show, that they belong to the same component. If    
         //both cells are label the same, a ring is found. Then the links   
         //of the tree cells are followed to trace back to the common parent
         //and the new ring, i.e. the set of bonds, is added to the result. 

         for (b = 0, color = NO_COLOR; b<nbonds; b++)
         {
             at1 = bonds[b][0]; at2 = bonds[b][1];
             if (at1 == 0 || at2 == 0) continue;  //ignore bonds with atom 0 
             if (tree[at1].color == NO_COLOR &&     //new component 
                 tree[at2].color == NO_COLOR)
             {
                 color++; tree[at2].link = b;
                 tree[at1].color = tree[at2].color = color;
             }
             else if (tree[at1].color == NO_COLOR)           //link first atom 
             {
                 tree[at1].color = tree[at2].color;
                 tree[at1].link = b;
             }
             else if (tree[at2].color == NO_COLOR)         //link second atom 
             {
                 tree[at2].color = tree[at1].color;
                 tree[at2].link = b;
             }
             else if (tree[at1].color != tree[at2].color)  //link two compnts. 
             {
                 new_color = tree[at1].color; old_color = tree[at2].color;
                 for (i = 0; i <= natoms; i++)
                     if (tree[i].color == old_color) tree[i].color = new_color;

                 tmp_link = tree[at2].link;     //trace the links of component 2 
                 tree[at2].link = b;            //and revert linkage             
                 while (tmp_link != UNLINKED)
                 {
                     at2 = TheOtherAtom(bonds[tmp_link], at2);
                     new_link = tmp_link;
                     tmp_link = tree[at2].link;
                     tree[at2].link = new_link;
                 }
             }
             else                            //ring found -> add it to list 
             {                               //trace back to root from both atoms 
                 for (trace = at1, level1 = 0; tree[trace].link != UNLINKED; level1++)
                     trace = TheOtherAtom(bonds[tree[trace].link], trace);
                 for (trace = at2, level2 = 0; tree[trace].link != UNLINKED; level2++)
                     trace = TheOtherAtom(bonds[tree[trace].link], trace);

                 if (level1 > level2)   //make path 1 the shorter one of the two 
                 {
                     tmp = level1; level1 = level2; level2 = tmp;
                     tmp = at1; at1 = at2; at2 = tmp;
                 }

                 p = NewBondSetNode(nbonds); p->next = result; result = p;
                 p->bond_set = ClearSet(p->bond_set); p->cardinality = 1;
                 p->bond_set = PutMember(p->bond_set, b);

                 for (i = 0; i<level2 - level1; i++) //trace back excess of long path 
                 {
                     p->bond_set = PutMember(p->bond_set, (unsigned)tree[at2].link);
                     p->cardinality++;
                     at2 = TheOtherAtom(bonds[tree[at2].link], at2);
                 }

                 while (at1 != at2)     //simultaneously trace back both paths 
                 {
                     p->bond_set = PutMember(p->bond_set, (unsigned)tree[at1].link);
                     p->cardinality++;
                     at1 = TheOtherAtom(bonds[tree[at1].link], at1);
                     p->bond_set = PutMember(p->bond_set, (unsigned)tree[at2].link);
                     p->cardinality++;
                     at2 = TheOtherAtom(bonds[tree[at2].link], at2);
                 }
             }       //else ring found 
         }       //for all bonds 

         free((char *)tree);
         return(result);
     }

     bond_set_node *SortRings(bond_set_node *list)
         /*
         * Sorts *list into descending order with respect to cardinality.
         
     {
         bond_set_node *p1, *p2;
         bit_set_t *set;
         int size;

         for (p1 = list; p1 != (bond_set_node *)NULL; p1 = p1->next)
         {
             for (p2 = p1->next; p2 != (bond_set_node *)NULL; p2 = p2->next)
                 if (p1->cardinality > p2->cardinality)  //out of order -> exchange 
                 {
                     size = p1->cardinality;
                     p1->cardinality = p2->cardinality;
                     p2->cardinality = size;
                     set = p1->bond_set;
                     p1->bond_set = p2->bond_set;
                     p2->bond_set = set;
                 }
         }

         return(list);
     }

     bond_set_node *CombineRings(bond_set_node *list)
// Combines pairs of rings until selfconsistency to get a list of smaller basis rings.
     {
         bond_set_node *p1, *p2;
         int size;
         bit_set_t *set, *tmp;
         int changed;
         int ntoggle;

         if (list == (bond_set_node *)NULL) return(list);

         srand(1);         //make sure that algorithm works reproducibly 
         ntoggle = 0;      //safeguard against infinite looping 
         set = NewSet(MaxMember(list->bond_set));  //assumes allocated size is 
         do                                        //the same for all rings!!! 
         {
             changed = FALSE;
             list = SortRings(list);  //loop over all pairs of different rings 
             for (p1 = list; p1 != (bond_set_node *)NULL; p1 = p1->next)
                 for (p2 = p1->next; p2 != (bond_set_node *)NULL; p2 = p2->next)
                 {
                     // check first if there is any overlap
                     if (IntersectionIsEmpty(p1->bond_set, p2->bond_set)) continue;;
                     set = SetExclusiveUnion(CopySet(set, p1->bond_set), p2->bond_set);
                     size = Cardinality(set);
                     if (size > 0 &&
                         (size <= p1->cardinality || size <= p2->cardinality))
                     {
                         if (p1->cardinality > p2->cardinality)
                         {
                             if (p1->cardinality > size || (rand() / 10) % 2)
                             {
                                 if (p1->cardinality == size)
                                     ntoggle++;
                                 else
                                     ntoggle = 0;
                                 changed = TRUE; p1->cardinality = size;
                                 tmp = p1->bond_set; p1->bond_set = set; set = tmp;
                             }
                         }
                         else
                         {
                             if (p2->cardinality > size || (rand() / 10) % 2)
                             {
                                 if (p2->cardinality == size)
                                     ntoggle++;
                                 else
                                     ntoggle = 0;
                                 changed = TRUE; p2->cardinality = size;
                                 tmp = p2->bond_set; p2->bond_set = set; set = tmp;
                             }
                         }
                     }
                 }
             if (ntoggle > 4) changed = FALSE;          // limit to 4 toggles 
         } while (changed);

         DisposeSet(set);
         return(list);
     }
*/     

 }// namespace StructureCheck
} // namespace RDKit
