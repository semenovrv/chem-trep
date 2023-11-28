#include <cmath>
#include "simple_molecule.h"
#include "edited_molecule.h"
#include "mol_struct_common.h"
//#include "mcdlutil.h"

#include <fstream>
#include <iostream>

namespace MolStruct {

	void alternateBonds(TSimpleMolecule & sm) {
		std::vector<int> aPosition, aCharge, aRad, nHydr, iA1, iA2, bondOrders;
		int nAtoms, nBonds, i,n;
		bool alternateNeed=false;
		
		for (i = 0; i < sm.nBonds(); i++) if (sm.getBond(i)->tb == 4) {
			alternateNeed = true;
			break;
		}
		if (!alternateNeed) return;
		nAtoms = sm.nAtoms();
		aPosition.reserve(nAtoms);
		aCharge.reserve(nAtoms);
		aRad.reserve(nAtoms);
		nHydr.reserve(nAtoms);
		for (i = 0; i < nAtoms; i++) {
			aPosition.push_back(sm.getAtom(i)->na);
			aCharge.push_back(0);
			aRad.push_back(0);
			nHydr.push_back(0);
		}

		nBonds = sm.nBonds();
		iA1.reserve(nBonds);
		iA2.reserve(nBonds);
		bondOrders.reserve(nBonds);
		for (i = 0; i < nBonds; i++) {
			iA1.push_back(sm.getBond(i)->at[0]);
			iA2.push_back(sm.getBond(i)->at[1]);
			n = sm.getBond(i)->tb;
			if (n == 4) n = 0;
			bondOrders.push_back(n);
		}
//		alternate(aPosition, aCharge, aRad,  nHydr, iA1, iA2, bondOrders, nAtoms, nBonds);
		for (i = 0; i < nBonds; i++) sm.getBond(i)->tb = bondOrders[i];

	}

    //pathdetermine.pas
	void pathCalc(const TSimpleMolecule & sm, int an1, int an2, std::vector<int> & list) {
		int i, j, n, k, m;
		bool test;
		std::vector<int> prevAtom, prevDistance;

		list.clear();
		prevAtom.resize(sm.nAtoms());
		prevDistance.resize(sm.nAtoms());
		for (i = 0; i < sm.nAtoms(); i++) {
			prevAtom[i] = -1;
			prevDistance[i] = 0;
		};
		prevAtom[an1] = sm.nAtoms()+1;
		prevDistance[an1] = 1;
		test = true;
		m = 0;
		while (test) {
			test = false;
			n = -1;
			m++;
			for (i = 0; i < sm.nAtoms(); i++) if ((prevAtom[i] >= 0) && (prevDistance[i] == m)) {
				for (j = 0; j < sm.getAtom(i)->nb; j++) {
					n = sm.getAtom(i)->ac[j];
					if (prevAtom[n] < 0) {
						prevAtom[n] = i;
						prevDistance[n] = m+1;
					}
					if (n == an2) break; else test = true;
				};
				if (n == an2) break;
			};
			if (n == an2) break;
		};
		if (n == an2) {
			list.push_back(n);
			test = true;
			while (test) {
				n = prevAtom[n];
				if ((n >= 0) && (n < sm.nAtoms())) list.push_back(n); else test = false;
			};
		};
		/*
		for (i = 0; i < sm.nAtoms(); i++) sm.getAtom(i)->enumerator = 0;
		n = 0;
		sm.getAtom(an1)->enumerator = 1;
		test = true;
		while (test && (list.size() == 0)) {
		    test = false;
		    n++;
		    for (i = 0; i < sm.nAtoms(); i++) if (sm.getAtom(i)->enumerator == n) {
			    for (j = 0; j < sm.getAtom(i)->nb; j++) {
				    k = sm.getAtom(i)->ac[j];
				    if (sm.getAtom(k)->enumerator == 0) {
					    sm.getAtom(k)->enumerator = n + 1;
					    test = true;
				    };
					if (k == an2) {
						list.push_back(an2);
						break;
					};
			    };
		    };
		    if (list.size() > 0) test = false;
		};// until(list.Count > 0) or (not test);



	    if (list.size()>0) {
		    k = list[0];
            n= sm.getAtom(k)->enumerator;
			test = true;
			while (test && (n > 0)) {
				test = false;
				for (j = 0; j < sm.getAtom(k)->nb; j++) {
					m = sm.getAtom(k)->ac[j];
					if (sm.getAtom(m)->enumerator == (n - 1)) {
						list.push_back(m);
						k = m;
						n--;
						test = true;
						break;
					};
				};
			};// until(not test) or (n == 0);
	    };
		*/
	};


	int  topoDistance(const TSimpleMolecule & sm, int an1, int an2) {
		std::vector<int> list;
		if ((an1 < 0) || (an2 < 0)) return 0;
	    pathCalc(sm, an1, an2, list);
        return list.size();
	};

	void removeSpaces(std::string & data) {
		while (data.length() > 0) {
			if (data.at(0) == ' ') data = data.substr(1); else
			if (data.at(data.length() - 1) == ' ') data = data.substr(0, data.length() - 1); else break;
		}
	}

	int determineAllPathesForPrediction(TSimpleMolecule & molIn, int an1, int an2, const std::vector<std::vector<int> *> & ringList, std::vector<int> & mainChain, const std::vector<adjustedlist> & bkIn) {

		//This method allows to remove spirocycles, attached to backbone
		std::vector<int> * ringBonds;
		std::vector<int> ringAtoms, usedRings;
		bool added;
		int i, j, n, countStore;

		mainChain.clear();
		pathCalc(molIn, an1, an2, mainChain);
		usedRings.resize(ringList.size());
		added = true;
		while (added) {
			added = false;
			countStore = mainChain.size();
			for (i = 0; i < ringList.size(); i++) if (usedRings[i] == 0) {
				ringBonds = ringList[i];
				molIn.bondListToAtomList(*ringBonds, ringAtoms);
				n = 0;
				for (j = 0; j < ringAtoms.size(); j++) if (findQuick(mainChain, countStore, ringAtoms[j]) >= 0) n++;
				if (n > 1) {  //skip spiro compounds-only fused rings
					added = true;
					usedRings[i] = 1;
					for (j = 0; j < ringAtoms.size(); j++) mainChain.push_back(ringAtoms[j]);
				};
			};
			if (added) quickSortIntegers(mainChain);
		}; // until not added;
	//below usedRings are used as temporary array
		usedRings.clear();
		if (mainChain.size() > 0) {
			usedRings.reserve(mainChain.size());
			usedRings.push_back(mainChain[0]);
			for (i = 1; i < mainChain.size(); i++) if (mainChain[i] != mainChain[i - 1]) usedRings.push_back(mainChain[i]);
			mainChain.resize(usedRings.size());
			for (i = 0; i < mainChain.size(); i++) mainChain[i] = usedRings[i];
			quickSortIntegers(mainChain);
		}
		return 0;
	};


	int determineAllPathes(const TSimpleMolecule molIn, int an1, int an2, std::vector<std::vector<int> *> & pathways) {
		//==0 result - OK, ==1-an1 at different fragment than an2, ==2 - an1 or an2 out of range
		TSimpleMolecule sm;
		std::vector<int> bnList;
		std::vector<int>* list0;
		int i, j, k, n1, n2, n;
		bool test;

		int result = 0;
		list0 = NULL;
		if ((an1 < 0) || (an2 < 0) || (an1 >= molIn.nAtoms()) || (an2 >= molIn.nAtoms())) {
			result = 2;
			return result;
		};

		for (i = 0; i < pathways.size(); i++) delete (pathways[i]);
		pathways.clear();
	

		sm.moleculeCopy(molIn);
		sm.defineAtomConn();
		sm.allAboutCycles();
		list0 = new std::vector<int>;
 		pathCalc(sm, an1, an2, *list0);
		if (list0->size() > 0) {
			quickSortIntegers(*list0);
			for (i = 0; i < sm.nBonds(); i++) if (sm.getBond(i)->db > 1) if ((findQuick(*list0, list0->size(), sm.getBond(i)->at[0]) >= 0)
				&& (findQuick(*list0, list0->size(), sm.getBond(i)->at[1]) >= 0)) bnList.push_back(i);
			pathways.push_back(list0);
			//another pathways
			for (i = 0; i < bnList.size(); i++) {
				n1 = sm.getBond(bnList[i])->at[0];
				n2 = sm.getBond(bnList[i])->at[1];

				test = false;
				for (j = 0; j < sm.getAtom(n1)->nb; j++) {
					if (sm.getAtom(n1)->ac[j] == n2) test = true;
					if (test) sm.getAtom(n1)->ac[j] = sm.getAtom(n1)->ac[j + 1];
				};
				sm.getAtom(n1)->nb = sm.getAtom(n1)->nb - 1;

				test = false;
				for (j = 0; j < sm.getAtom(n2)->nb; j++) {
					if (sm.getAtom(n2)->ac[j] == n1) test = true;
					if (test) sm.getAtom(n2)->ac[j] = sm.getAtom(n2)->ac[j + 1];
				};
				sm.getAtom(n2)->nb = sm.getAtom(n2)->nb - 1;

				//new path
				list0 = new (std::vector<int>);
				pathCalc(sm, an1, an2, *list0);
				//restore
				sm.getAtom(n1)->nb = sm.getAtom(n1)->nb + 1;
				sm.getAtom(n1)->ac[sm.getAtom(n1)->nb - 1] = n2;
				sm.getAtom(n2)->nb = sm.getAtom(n2)->nb + 1;
				sm.getAtom(n2)->ac[sm.getAtom(n2)->nb - 1] = n1;
				//sm.getAtom.nb[n2]= sm.getAtom[n2].NB + 1;
				//sm.getAtom.AC[n2, sm.getAtom[n2].NB]= n1;
				//checking list
				if (list0->size() > 0) {
					quickSortIntegers(*list0);
					for (j = 0; j < pathways.size(); j++) {
						//tempList = pathways[j];
						if (pathways[j]->size() == list0->size()) {
							test = true;
							for (k = 0; k < pathways[j]->size(); k++) if ((*pathways[j])[k] != (*list0)[k]) {
								test = false;
								break;
							};
							if (test) {
								delete(list0);
								list0 = NULL;
							};
						};
						if (!list0) break;
					};
					if (list0) pathways.push_back(list0);
				}
				else {
					delete(list0);
					list0 = NULL;
				};
			};
		} else {
		    delete(list0);
			result= 1;
		};
	};

	int determineAllPathes(TSimpleMolecule & sm, int an1, int an2, std::vector<int>& mainChain, const std::vector<adjustedlist> & bk) {
	//==0 result - OK, ==1-an1 at different fragment than an2, ==2 - an1 or an2 out of range

		int i, j, k, n, cCount;
		std::vector<int> tempArray;
		bool test;

		int result = 0;
		if ((an1<0) || (an2<0) || (an1>=sm.nAtoms()) || (an2>=sm.nAtoms())) {
			result = 2;
			return result;
		};
		mainChain.clear();
		pathCalc(sm, an1, an2, mainChain);
		if (mainChain.size() > 0) {
			test = true;
			while (test) {
				test = false;
				cCount = mainChain.size();
				quickSortIntegers(mainChain);
				for (i = 0; i < sm.nAtoms(); i++) if (findQuick(mainChain, cCount, i) >= 0) for (j = 0; j < bk[i].nb; j++) {
					n = bk[i].adjusted[j];
					if (sm.getBond(n)->db > 1) {
						if (sm.getBond(n)->at[0] == i) k = sm.getBond(n)->at[1]; else k = sm.getBond(n)->at[0];
						if (findQuick(mainChain, cCount, k) < 0) {
							mainChain.push_back(k);
							test = true;
						};
					};
				};
			};  // until not test;
			tempArray.reserve(mainChain.size());
			tempArray.push_back(mainChain[0]);
			for (i = 1; i < mainChain.size(); i++) if (mainChain[i] != mainChain[i - 1]) tempArray.push_back(mainChain[i]);
			mainChain.clear();
			mainChain.reserve(tempArray.size());
			for (i = 0; i < tempArray.size(); i++) mainChain.push_back(tempArray[i]);
		} else {
			result = 1;
		};
	};
	

	void chainAtomList(std::vector<std::vector<int> *> & pathways, std::vector<int>& atomList) {
		std::vector<int> list;
		int i, j;
	
		for (i = 0; i < pathways.size(); i++) for (j = 0; j < pathways[i]->size(); j++) list.push_back((*pathways[i])[j]);
	
		quickSortIntegers(list);
		atomList.clear();
		atomList.reserve(list.size());
		if (list.size() > 0) {
			atomList.push_back(list[0]);
			for (j = 1; j < list.size(); j++) if (list[j] != list[j - 1]) atomList.push_back(list[j]);
		};
	};

	//radprocess.pas
	void makeEquivalentList(const TSimpleMolecule & molecule, std::vector<int> & equivalenceList, bool onlyTopological) {
		molecule.makeEquivalentList(equivalenceList, onlyTopological);
	};


	int indexOf(const std::vector<int>& values, int value) {
		int result = -1;
		int i;
		for (i = 0; i < values.size(); i++) if (values[i] == value) {
			result = i;
			break;
		}
		return result;
	}

	int indexOf(const std::vector<std::string>& values, std::string value) {
		int result = -1;
		int i;
		for (i = 0; i < values.size(); i++) if (values[i] == value) {
			result = i;
			break;
		}
		return result;
	}


	bool saveTextFile(const std::string fileName, const std::vector<std::string> & textData) {
		bool result = false;
		std::ofstream fileData(fileName);
		std::string s;
		int i;

		if (fileData.is_open()) {
			for(i=0; i<textData.size(); i++){
				s = textData[i];
				fileData << s << std::endl;
			};
			fileData.close();
			result = true;
		};
		return result;
	};

	bool saveTextFile(const std::string fileName, const std::string & textData) {
		bool result = false;
		std::ofstream fileData(fileName);

		fileData << textData << std::endl;
		fileData.close();
		result = true;
		return result;
	};


	bool loadTextFile(const std::string fileName, std::vector<std::string> & textData) {
		bool result = false;
		std::ifstream fileData(fileName);
		std::string s;

		textData.clear();
		if (fileData.is_open()) {
			//result = readMolfile(fileMol);
			while (std::getline(fileData, s)) textData.push_back(s);
			fileData.close();
			result = true;
		};
		return result;
	};

	std::string format(const std::string & fmtString, const std::string & argument) {
		std::string result = fmtString;
		int n = result.find("{}");
		if (n != std::string::npos) {
			result = result.substr(0, n) + argument + result.substr(n + 2);
		}
		return result;
	}

	std::string format(const std::string & fmtString, const std::vector<std::string> & arguments) {
		std::string result = fmtString;
		int k = 0;
		bool test = true;
		while (test) {
			int n = result.find("{}");
			if (n == std::string::npos) test = false; else {
				result = result.substr(0, n) + arguments[k] + result.substr(n + 2);
			};
			k++;
			if (k >= arguments.size()) test = false;
		};
		return result;
	}

	std::string formatPrecision(double value, int precision) {
		double multiplicator = 1;
		int i,n;
		std::string result,exponent;

		result = std::to_string(value);
		exponent = "";
		n = result.find("e");
		if (n == std::string::npos) n = result.find("E");
		if (n != std::string::npos) {
			exponent = result.substr(n);
			result = result.substr(0, n);
		}
		else {
			if (precision > 0) for (i = 0; i < precision; i++) multiplicator = 10 * multiplicator;
			result = std::to_string(round(multiplicator*value) / multiplicator);
		};
		n = result.find(".");
		if (n != std::string::npos) {
			if (precision == 0) result = result.substr(0, n); else result = result.substr(0, n + precision + 1);
		}
		result = result + exponent;
		return result;
	}


};//Molspace