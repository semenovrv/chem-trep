
#ifndef _HDR_MOLUTILS_
#define _HDR_MOLUTILS_

#include <vector>
#include <string>
#include "simple_molecule.h"
#include "mol_struct_common.h"



namespace MolStruct {
	void alternateBonds(TSimpleMolecule & sm);
	//path determine.pas
	void pathCalc(const TSimpleMolecule & sm, int an1, int an2, std::vector<int> & list);
	int  topoDistance(const TSimpleMolecule & sm, int an1, int an2);
	int determineAllPathesForPrediction(TSimpleMolecule & molIn, int an1, int an2, const std::vector<std::vector<int> *> & ringList, std::vector<int> & mainChain, const std::vector<adjustedlist> & bkIn);
	int determineAllPathes(TSimpleMolecule molIn, int an1, int an2, std::vector<std::vector<int> *> & pathways);
	int determineAllPathes(const TSimpleMolecule & sm, int an1, int an2, std::vector<int>& mainChain, const std::vector<adjustedlist> & bk);
	void chainAtomList(std::vector<std::vector<int> *> & pathways, std::vector<int>& atomList);
	//radprocess.pas
	void makeEquivalentList(const TSimpleMolecule & molecule, std::vector<int> & equivalenceList, bool onlyTopological);
	int indexOf(const std::vector<int>& values, int value);
	int indexOf(const std::vector<std::string>& values, std::string value);
	bool saveTextFile(const std::string fileName, const std::vector<std::string> & textData);
	bool saveTextFile(const std::string fileName, const std::string & textData);
	bool loadTextFile(const std::string fileName, std::vector<std::string> & textData);
	std::string format(const std::string & fmtString, const std::string & argument);
	std::string format(const std::string & fmtString, const std::vector<std::string> & arguments);
	std::string formatPrecision(double value, int precision);
	void removeSpaces(std::string & data);

};//MolStruct

#endif
