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
	Class TSingleBond was separated from simple_molecule.h.
*/

#ifndef _HDR_SINGLE_BOND_
#define _HDR_SINGLE_BOND_

#include <string>
#include <iostream>
#include "mol_struct_common.h"

#define NBONDTYPES 11

namespace MolStruct {

static bool checkAromatic=true;

class  TSingleBond  {
public:


	/*!
	Bond type:
	    0 - single or triple
		1 - single
		2 - double
		3 - triple
		4 - aromathic (query)
		5 - single/double (query)
		6 - coordinational
		7 - double or triple
		8 - any (query)
		9 - up
		10 - down
		11 - either
	*/

	enum		BondType {
		Undefined,
		Single,
		Double,
		Triple,
		Aromatic,
		SingleOrDouble,
		Coordinational,
		Dummy,
		Any,
		Up,
		Down,
		Either
	};

	short int tb;

	int at[2];     /*Bond definition-atoms number in array ATOM*/
	short int db;        /*Ring/Chain conditions*/
	short int bstereo;   /* =0 - no, =1 - E , = 2 - Z, =3 E/Z*/
	
	/*	=1 - Chain (query)
		=2 - Ring  (query)
		=3 - CIS/TRANS
		=4 - CIS/TRANS+Chain;
		=5 - CIS/TRANS+Ring */
	short int special;

	int		 enumerator;

	TSingleBond*		clone() const;
	int					bondConversion();
	void					bondCopy(TSingleBond * source);
	static bool			bondEquivalent(TSingleBond * sBond, TSingleBond * qBond, bool compareHard);
	int					getValence();

	 TSingleBond(void);
	 TSingleBond(int atom0, int atom1, int type = 1);
	 virtual ~TSingleBond(void){};

	 friend std::ostream& operator<<(std::ostream& out, const TSingleBond& sb);
};


} // namespace MolStruct

#endif