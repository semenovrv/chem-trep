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
	Class TSimpleMolecule was separated from simple_molecule.h.
*/

#ifndef _HDR_SIMPLE_MOLECULE_
#define _HDR_SIMPLE_MOLECULE_

#include "single_atom.h"
#include "single_bond.h"

#include <string>
#include <vector>


namespace MolStruct {

const int CYCLE_MAX_SIZE = 10;


class TSimpleMolecule {
private:
	bool aromatic(int cycleSize, const std::vector<int> bondList, std::vector<int>& arom);
	void twoAtomUnitVector(int na1, int na2, double & xv, double & yv, const std::vector<int>atomDefine);
	void defC(int& currNumDef, int baseCycle, int atomClean, std::vector<int>& cycleDefine,
		std::vector<int>& atomDefine, std::vector<int>& atomCycle, std::vector<int>& cycleAddress, 
		std::vector<int>& cycleSize, std::vector<int>& dsATN, std::vector<int>& dsTP, 
		std::vector<int>& dsSC, std::vector<int>& dsNA1, std::vector<int>& dsNA2);
	void defA(int& currNumDef, int atomClean, int sPN, int baseCycle, std::vector<int>& atomDefine, const std::vector<int> listAtomClean,
		std::vector<int>& cycleDefine, std::vector<int>& cycleSize, std::vector<int>& cycleAddress, std::vector<int>& atomCycle,
		std::vector<int>& dsATN, std::vector<int>& dsTP, std::vector<int>& dsNA1, std::vector<int>& dsNA2);
	void canonizeCycle(int ringSize, std::vector<int> & bondList);
	//Chain rotate members
	double atomDistanceMetric(int an);
	bool bondsOverlapped(int bN1, int bN2, double delta);
	bool threeBondResolve(int an, int bondExcluded, double& xv, double& yv, std::vector<adjustedlist> * bkExt) const;
	bool unitVectorCoincident(int aN, double xV, double yV) const;
	void makeStandardStructure();
    void correctValencies();
    bool deleteHydrogenDecCharge(int an, int anInc);
    void tautomerAnalize();
    void correctAzulenes(const std::vector<adjustedlist> & bk);
    int  chargeConversion(int atn);
    int  valencyConversion(int atn);
    int  allAtAtom(int atn);
    int  bondConversion(int bnb);
    void processCodeList();
	bool bondConnectsCH2(int bondNo);
	std::string getSimplePolymerMolBlock() const;
protected:
	std::vector<TSingleAtom*>			fAtom;
	std::vector<TSingleBond*>			fBond;

public:
	std::ostream							*refofs;
	
	//Very bad!!! Make it private!!!
	bool						fIOPT10;   //charge sensitivite
	bool                        fIOPT11;   //isotope
	int							fIOPT12;   //stereo bond change
	bool						fIOPT13;   //semipolar bond as double
	int	  						unC1, unC2, unC3;
	std::string                 inChIKey;
	int							userData, chiral;
	unsigned int				fScreenData[5];
	
	std::vector<fragmentCode> mixtureCodeList;
	std::vector<twoSphereRecord> fragmentCodeList;
	//std::string inChIKey;


	int nAtoms() const;
	int nBonds() const;
	TSingleAtom * getAtom(int index);
	TSingleBond * getBond(int index);

	const TSingleAtom * getAtom(int index) const;
	const TSingleBond * getBond(int index) const;
	void defineAtomConn();
	void defineBondConn(std::vector<adjustedlist> & bondConnection) const;
    void calculateRMGTypes(std::vector<adjustedlist> * bC);

	void newB(const std::vector<adjustedlist> & bk, int bnum, int anum, int & total, int * e, int * e1) const;
	void singleVawe(const std::vector<adjustedlist> & bk, std::vector<int> & alreadyDefined,
		std::vector<int> & prevBond, std::vector<int> & prevAtom,
		int & nPrev, std::vector<int> & curBond, std::vector<int> & curAtom) const;
	void vaweBond(int bondN, const std::vector<adjustedlist> & bk, int & ringSize, std::vector<int> & bondList, int maxSize=0) const;
	void atomDistanceList(int atomNo, int prohibitedAtomNo, std::vector<int> & atomDistance) const;
	void allAboutCycles(int maxSize = 0);
	void zeroMolDynProperties();
	void redraw(const std::vector<int>listAtomClean, const std::vector<int>listBondClean, 
		int & atomClean, int & bondClean, int spn, int sCHA1, int sCHB1, bool iOPT7);
	void clear();
	void readConnectionMatrix(const std::vector<int>iA1, const std::vector<int>iA2, int nAtoms, int nBonds);
	void readConnectionMatrix(const std::vector<int>iA1, const std::vector<int>iA2, const std::vector<double>rx, const std::vector<double>ry, int nAtoms, int nBonds);
	void redrawMolecule();
	void getMolfile(std::ostream & data, bool writeend=false, std::vector<double> * zCoor=NULL) const;
	std::string getMolfile() const;
	void saveMolfile(const std::string fileName);
    bool readMolfile(std::istream & data);
	bool readMolfile(const std::vector<std::string> & data);
	bool readMolfile(const std::string fileName);
	int  loadAllFromStream(std::istream & data);
	void removeHydrogen(std::vector<int> * qHydr, std::vector<int> * qEnumerator);
	void atomBondChange();
	bool stereoBondChange();


	bool checkOverlapped();
	double averageBondLength() const;
	double averageAtomDistance();
	void flipSmall(int cHB);
	int listarSize() const;
	bool makeFragment(int& n, std::vector<int>& list, int aT, int aTEx) const;
	void getFragment(int & an, int & bn, TSimpleMolecule & fragment);
	void moleculeCopy(TSimpleMolecule & source);
	void moleculeCopy(const TSimpleMolecule & source);
	int correctOverlapped();
	void normalizeCoordinates(double aveBL);
	void makeEquivalentList(std::vector<int>& equivalenceList, bool isTopologyOnly) const;
	void bondEnlarge(int bN);
	void deleteAtom(int index);
	void deleteBond(int index);
	void addAtom(TSingleAtom * sa);
	void addAtom(int na, int charge, double rx, double ry);
	void addBond(TSingleBond * sb);
	void addBond(int tb, int at1, int at2);
	void addFragmentAtPos(TSimpleMolecule & otherMol, int an, double posX, double posY);
	void setCoordinatesString(const std::string value);
	void setCoordinatesStringEx(const std::string value);
	void unitVector(int aN, double& xV, double& yV) const;
	void bondUnitVector(int bn, double& xv, double& yv) const;
	double bondLength(int index) const;
	int singleAtomicDescriptor(int aNumber,int bNumber, bool useEnumerator);
	int hasOverlapped(double delta, bool findFirst);
	int getNH(int atomNo) const;
	int getExplicitH(int atomNo) const;
    int cyclesStatistic(std::vector<int>& nCyclesSize, std::vector<int>& cycleDetails, bool testEnumerator);
    int getNGauscheCarbons();
    int getNOrtho();
    void addAllHydrogens(bool forStereo);
    bool addDefaultHydrogen(int atn, bool forStereo);
    bool isCoordinatesBad(int atn);
    void assignBondFromStereodescriptor(int atomNo, int atomChirality, bool ignoreH);
	void calculateAllIndeces();
    void extractFragmentFast(const std::vector<int> aList, TSimpleMolecule & smExtracted);
    bool addMoleculeNoScale(TSimpleMolecule * molecule);
    void createFragmentList(std::vector<TSimpleMolecule *> & list);
	int  makeUniqueSymbol(const std::vector<adjustedlist> & bondConnection);
	int  makeSecondSymbol(const std::vector<adjustedlist> & bondConnection);
	int  setResourceString(const std::string value);
	int  setJMEString(const std::string value);


	TSimpleMolecule * extractFragment(int atomN, std::vector<int> * enumerator); 
	void extractFragment(int atomNo, int topoDistance, int firstAtomCharge, TSimpleMolecule & saveStructure);

//	void extractFragment(int atomNo, int topoDistance, int firstAtomCharge, TSimpleMolecule & saveStructure, std::vector<adjustedlist> * bk);

	

    int processCharges();
    int  fragmentSecond(int sphere, int att, int secAt, const std::vector<int>a, const std::vector<int>b, 
      const std::vector<adjustedlist> & bk, std::vector<int> & wSphere);

	bool formulaIdentical(TSimpleMolecule * other);
	void getHomolog(TSimpleMolecule * sMol);
	bool isHomolog(TSimpleMolecule * other);



	std::string				getAtomSymbol(int atAtom);
	double					getAtomMass(int atAtom);
	double					getMolWeight();
	std::string				getMolformula(bool isHTML=false) const;
	int						getNumberOfAtoms(int na);
	void					removeExplicitHydrogens(void);
	void					rescaleToLength(double newAvgBondLength);
    void					init();
	TSimpleMolecule&		operator=(const TSimpleMolecule& sm);

	//molecular dynamic operators
    void fillMolDynStructureProperties();

	bool isSP3Atom(int atomNo, const std::vector<adjustedlist> & bondConnection) const;
	bool isSP2Atom(int atomNo, const std::vector<adjustedlist> & bondConnection) const;
	bool isSPAtom(int atomNo, const std::vector<adjustedlist> & bondConnection);
	bool isAromaticAtom(int atomNo, const std::vector<adjustedlist> & bondConnection);
	bool isPlanarRing(int cycleSize, const std::vector<int> & bondList, const std::vector<adjustedlist> & bondConnection);
	int  topologicalRadius(int atomIndex);
	void extractSphericalEnviroment(int atomIndex, int topoDistance, TSimpleMolecule & sMol);
	void extractZeroBonds(int bondIndex, TSimpleMolecule & sMol, std::vector<adjustedlist> * bk);
	bool isAromatic(const std::vector<int> & bondList) const;
	void bondListToAtomList(const std::vector<int> & bondList, std::vector<int> & atomList) const;
	int getOppositeLabel(int atomNo) const;
	int getOppositeAttachedAtom(int atomNo) const;
	int bondNumber(int an1, int an2) const;
	bool atomInRing(int atomNo, const std::vector<int> & bondList) const;
	void cyclesCalculate(int & nTotalCycles, int & nAromFive, int & nAromSix, int & nCondenced, std::vector<std::vector<int> *> * ringList) const;
	bool hasDoubleBond(int atomNo,const  std::vector<adjustedlist> & bk) const;
    void removeFarPendant(const std::vector<std::vector<int> *> & ringList, const std::vector<int> & mainChain, int backboneDistance);
	bool extractFragment(int atomNo, int acyclicBondNo, TSimpleMolecule & smExtracted) const;
	int  topoDistance(int an1, int an2) const;
	//std::string calculateInChIKey();
	//std::string calculateInChIKeyEx() const;
	void makeStandardPolymer(const std::vector<int> & mainChain);



	TSimpleMolecule(const TSimpleMolecule& sm);
	TSimpleMolecule();

	virtual ~TSimpleMolecule()
	{
		clear();
	};

};

bool fragmentsIdentical(int fromSphere1, int fromSphere2, int sphereCount, const std::vector<int> & sphereList, const std::vector<int> & eqList);
void analyzeRotors(TSimpleMolecule & sm, std::vector<std::string> & explanation, std::vector<int> & symnumber);
bool listPresent(const std::vector<std::vector<int> *> & ringData, const std::vector<int> & tempList);
bool bondPresent(std::vector<int> & bondList, int bn);
bool incrementValues(std::vector<int>& currentValues, const std::vector<int> maxValues);
void writeSDF(const std::string & fieldName, const std::string & data, std::vector<std::string> & fileData);



} // namespace MolStruct

#endif