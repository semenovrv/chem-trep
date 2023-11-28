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
	Class TEditedMolecule was separated from simple_molecule.h.
*/


#ifndef _HDR_EDITED_MOLECULE_
#define _HDR_EDITED_MOLECULE_

#include <string>
#include "simple_molecule.h"
#include "scalar.h"

namespace MolStruct {

class TEditedMolecule: public TSimpleMolecule {
protected:
	std::vector<adjustedlist>			queryBKNew;
	std::vector<adjustedlist> structureBKNew;
	std::vector<int>		queryQHydr;
	std::vector<int>		queryAGer;
	std::vector<int>		queryBQCounter;
	std::vector<int>		queryCurrentAssignment;
	std::vector<int>		queryEnum;
	std::vector<int>		queryInverse;

	bool						queryStereoQ,fIsQueryPrepare,queryTestSpace;

	std::vector<int>		aSTested;
	std::vector<int>		bSTested;
	std::vector<int>		bSTestedStore;
	TScalar                 queryStas;

	//bool			aEQ[NATOMSMAX][NATOMSMAX];
	//bool			bEQ[NBONDSMAX][NBONDSMAX];
	void			directBondAss(int& bnq, bool& test, bool& test1, /*const bool beq[NBONDSMAX][NQUERYMAX],
						const bool aeq[NATOMSMAX][NQUERYMAX]*/ const std::vector<bool> & beq, const std::vector<bool> & aeq, std::vector<int>& bqcounter,
						std::vector<int>& aqtested, std::vector<int>& bstested,
						std::vector<int>& bqtested, std::vector<int>& astested, 
						const std::vector<int> ager, const std::vector<adjustedlist> & bqconn,
						const std::vector<adjustedlist> & bsconn, TSimpleMolecule *smol, TScalar& stas);
	bool			allQueryPresent(const std::vector<int> qA, const std::vector<int> qB, int nA, int nB);
	bool            getEquivalentValue(int indexQuery, int indexStructure, int nQueryElements, const std::vector<bool>&matrix);
	void            setEquivalentValue(int indexQuery, int indexStructure, int nQueryElements, bool value, std::vector<bool>&matrix);

	double		distance(int atom1, int atom2) const;
public:
//	static int const		NOOTHER_MASK=1; 
//	static int const		AROMATIC_MASK=NOOTHER_MASK << 1;
//	static int const		EXACTNUMBER_MASK=AROMATIC_MASK << 1;
	bool                    compareHard;
	bool                    rejectAnyInStructure;
	std::vector<int>		queryAQTested;
	std::vector<int>		queryBQTested;
	std::vector<int>		*fIncludedList;
	double					xShiftResize,yShiftResize,scaleResize;  //Is used to store RectResize values

	int						prepareQuery(TSimpleMolecule & sMol);
	bool					fragmentSearch(TSimpleMolecule * molecule1, std::vector<int>* bondLabel, bool useFragmentNo=false, bool useAddInfo=false);
	void					findAllInsertions(TSimpleMolecule * molecule1, std::vector<std::vector<int> *> & all);

	TEditedMolecule*		extractFragment(int atomN, std::vector<int> * enumerator);
	void					addAsTemplate(TSimpleMolecule& fragmentMol, int thisAN, int smAN, int thisBN, int smBN, bool isAddition);
	int						addFragment(TSimpleMolecule & eMolecule, int naDEF, int cha, int chb,
												int chb1, std::vector<int>& list, double xOldCenter,
												double yOldCenter, double xNewCenter,
												double yNewCenter, double scale, double cFi, double sFi,
												int buttonStatus, bool clearEnumerator);

	int						makePoly(int size, int cha1, int chb1, const double& xC, const double& yC);
	void						rectResize(double XMin, double YMin, double XMax, double YMax);


	static int sproduct(TSimpleMolecule & molecule, int br, int i1, int i2);


	TEditedMolecule() : TSimpleMolecule()
	{
		fIncludedList = NULL;
		fIOPT10 = true;   //charge sensitivite
		fIOPT11 = false;   //isotopes difference
		fIOPT12 = 2;   //stereo bond change
		fIOPT13 = true;  
		compareHard=false;
		rejectAnyInStructure=false;
		fIncludedList=NULL;
		queryStereoQ=false;
		fIsQueryPrepare=false;
		queryTestSpace=false;
	};

	friend class TAtomCenteredMolecule;
};

void createSingleMolecule(TSimpleMolecule * largeMolecule, std::vector<TSimpleMolecule *> & molList);
void decomposeSingleMoleculeNoCoordinates(TSimpleMolecule * largeMolecule, std::vector<TSimpleMolecule *> & molList, bool removeHydrogens);



} // namespace MolStruct

#endif