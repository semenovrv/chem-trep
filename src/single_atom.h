/*-*-C++-*-

**********************************************************************
Copyright (C) 2007,2008 by Sergei V. Trepalin sergey_trepalin@chemical-block.com
Copyright (C) 2007,2008 by Andrei Gakh andrei.gakh@nnsa.doe.gov

This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************
*/
/*
  Diagram is generated using templates, which are stored in SD file templates.sdf
  The SD file is usual SD file, which contain chemical structures and might contain data.
  Only chemical structures are used. Subgraph isomorphisme search is executed and coordinates
  of atoms are determined from templates. See Molecules, 11, 129-141 (2006) for algorithm decription.
  Structures in SD file are converted in next manner:
  1. All atoms, except explicit hydrogens, are replaced with generic ANY_ATOM (matched with any atom in subgraph isomorphisme search)
  2. All bonds are replaces with generic ANY_BOND, which can be matched with any bond in molecule
  3. All hydrogen are removed, but they are used for search-query and structure atom matching is believed fo be
     sucessfukk if chemical structure contains more or equal number of hydrogens, than query. Using explicitly-defined hydrogens
	 on query enables ones to remove substitutors attachment for atom, which are sterically hidden on templates
  if the file will not be found, predefined templates will be used
*/
/*
	Class TSingleAtom was separated from simple_molecule.h.
*/

#ifndef _HDR_SINGLE_ATOM_
#define _HDR_SINGLE_ATOM_

#include <string>
#include <vector>
#include <iostream>
#include "mol_struct_common.h"

#define RMG_QUERY_SIZE 2

namespace MolStruct {

class  TSingleAtom  {
public:
	short int na;   /*Atom position in the Periodic Table (nucley charge)*/
	short int nv;   /*Current atom's valence*/
	short int nc;   /*-9..+9 atom's charge*/
	short int iz;   /*Isotope difference between round(AtomMass)*/
	double    rx;   /*Internal X-coordinate representation*/
	double    ry;   /*Internal Y-coordinate representation*/
	short int rl;   /*radical label*/
	short int nb;   /*Number of neightboring atoms*/
	short int currvalence;
	/*Valence, excluding non-defined hydrogen
	atoms. Equal sum of valences for each  bond type.*/
	short int special;
	/*Special=0 - no special definition
	=1 - no other atoms, except explicitly defined, can
	be connected with the atom (query structure)*/
	int ac[CONNMAX];
	/*atom's numbers in array ATOM, which are
	connected with the atom selected*/
	short int astereo;   /*=0 - no stereo, =1 - R, = 2 -S,  = 3 - unknown*/

	std::string anum;
	int enumerator;
	int fragIndex;

	std::vector<int> * atomList;
	bool rejectInList;
	
	int rmgAtom[RMG_QUERY_SIZE];
	
	//bool processAddInfo;
	short int hybridization;  //=0 - no meaning, 1-sp3, 2-sp2 3-sp
	short int aromatic; //=0-no meaning 1-aromatic 2-non-aromatic
	short int chirality; //=0-no meaning =1-alpha =2-beta
	short int coordination; //=0-no meaning >0 - number of bonds, attached to this atom
	short int rKind1;  //=0 - no meaning =1-planar =2-non planar =3 any ring
	short int rSize1;  //=0 - no meaning =1 - any ring =3... ring size
	short int rKind2; 
	short int rSize2; 
	short int rKind3; 
	short int rSize3; 
	short int rKind4;
	short int rSize4;
	short int rKind5;
	short int rSize5;
	short int rKind6;
	short int rSize6;


	int encoder();
	int chargeConversion();
	int valencyConversion();
	int allAtAtom();
	void atomCopy(TSingleAtom * source);
	void addToList(int an);
	static bool atomEquivalent(TSingleAtom * structure, TSingleAtom * query,
						int nHStr, int nHQuery, bool chargeSensitivity, bool isotopeSensitivity,
						bool compareHard, bool rejectAnyInStructure, bool isUnsaturated, bool useAddInfo=false);
   static bool atomMolDynEquivalent(TSingleAtom * structure, TSingleAtom * query);
   static bool ringConditionsEquivalent(int queryKind, int querySize, int structureKind, int structureSize);
						
	static int chargeDeltaValency(int atomNo);

	TSingleAtom*				clone(void) const;

	void setRingProperty(int ringSize, bool isPlanar);


	TSingleAtom(void);
	virtual ~TSingleAtom(void);

	/*!Output operator.*/
	friend std::ostream& operator<<(std::ostream& out, const TSingleAtom& sa);
};




} // namespace MolStruct

#endif
