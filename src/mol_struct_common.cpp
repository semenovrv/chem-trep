/*****************************************************************************
 *
 *   Copyright (c), Kinetic Technologies Ltd. 2003-2011    All Rights Reserved.
 *
 *   Author	: Denis Shirabaikin, KINTECH, Moscow, Russia
 *
 *   Project	: MolStructure
 *
 *   $Revision:  $
 *   $Date:  $
 *   $Author:  $
 *   @(#) $HeadURL:  $
 *
 *****************************************************************************/

#include <vector>
#include <string>
#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>

#include "simple_molecule.h"
#include "single_atom.h"
#include "single_bond.h"
//#include "template_redraw.h"

using namespace std;

namespace MolStruct {


    //Hydrogen valencies. Zero dummy element is the first 
	const int hVal[NELEMMCDL] = {  
	0,1,0,0,0,3,4,3,2,1,
	0,0,0,3,4,3,2,1,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,4,3,2,1,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	2,3,2,1,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,2,0,2,1,0,1,2,0,
	0,0,0,0,0,0,0,0,0,0,//Fn
	0,0,0,0,1,1,0,0,0,1,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0};

	const int maxVal[NELEMMCDL] = {
    0,1,0,1,2,4,4,5,2,1,
	0,1,2,4,4,6,6,7,0,1,
	2,3,4,5,6,7,6,4,4,2,
	2,3,4,5,6,7,8,1,2,3,
	4,5,6,7,8,6,6,2,2,3,
	4,5,6,7,8,1,2,3,4,4,
	3,3,3,3,3,4,3,3,3,3,
	3,3,4,5,6,7,8,6,6,3,
	2,3,4,5,6,7,8,1,2,3,
	4,5,6,6,6,6,3,4,3,3,//Fm
	3,3,1,1,1,1,0,0,0,1,
	0,8,1,8,5,0,0,0,0,0,
	0,8,8,8,8,8,8,8,8,8,
	8,8,8,8,8};

    const int chargeVal[NELEMMCDL] = {  //0 - dummy
   0,-1,-1,-1,-1,-1,-1, 1, 1, 1,-1, //Ne
  -1,-1,-1,-1, 1, 1, 1,-1,-1,-1, //Ca
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1, //Zn
  -1,-1, 1, 1, 1,-1,-1,-1,-1,-1, //Zr
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1, //Sn
   1, 1, 1,-1,-1,-1,-1,-1,-1,-1, //Nd
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1, //Yb
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1, //Hg
  -1,-1, 1, 1, 1,-1,-1,-1,-1,-1, //Th
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1, //Fm
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1, //RMG 
  -1,-1,-1};


   const string aSymb[NELEMMCDL] = {"0",
    "H" ,"He","Li","Be","B" ,"C" ,"N" ,"O" ,"F" ,"Ne",
    "Na","Mg","Al","Si","P" ,"S" ,"Cl","Ar","K" ,"Ca",
    "Sc","Ti","V" ,"Cr","Mn","Fe","Co","Ni","Cu","Zn",
    "Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y" ,"Zr",
    "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn",
    "Sb","Te","I" ,"Xe","Cs","Ba","La","Ce","Pr","Nd",
    "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb",
    "Lu","Hf","Ta","W" ,"Re","Os","Ir","Pt","Au","Hg",
    "Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th",
    "Pa","U" ,"Np","Pu","Am","Cm","Bk","Cf","Es","Fm",
    "Md","No","Lr","D" ,"T" ,"G" ,"0" ,"Xx","Eh","lp"  ,
    "M" ,"X" ,"A" ,"Q" ,"*"  ,""  ,""  ,""  ,""  ,"",
     "R!H","Cbf","Cdd","Sis","Sid","Cs","Cd","Ct","Cb","CO",
	 "Os","Od","Ss","R"};

	const double aMass[NELEMMCDL]  = { 0.0,
		1.00794,    4.002602, 6.941,     9.012182,  10.811,    12.011,   14.00674,  15.9994,
		18.9984032, 20.1797,  22.989768, 24.3050,   26.981539, 28.0855,  30.973762, 32.066,
		35.4527,    39.948,   39.0983,   40.078,    44.955910, 47.88,    50.9415,   51.9961,
		54.93805,   55.847,   58.93320,  58.69,     63.546,    65.39,    69.723,    72.61,
		74.92159,   78.96,    79.904,    83.80,     85.4678,   87.62,    88.90585,  91.224,
		92.90638,   95.94,    97.9072,   101.07,    102.90550, 106.42,   107.8682,  112.411,
		114.82,     118.710,  121.75,    127.60,    126.90447, 131.29,   132.90543, 137.327,
		138.9055,   140.115,  140.90765, 144.24,    144.9127,  150.36,   151.965,   157.25,
		158.92534,  162.50,   164.93032, 167.26,    168.93421, 173.04,   174.967,   178.49,
		180.9479,   183.85,   186.207,   190.2,     192.22,    195.08,   196.96654, 200.59,
		204.3833,   207.2,    208.98037, 208.9824,  209.9871,  222.0176, 223.0197,  226.0254,
		227.0278,   232.0381, 231.0359,  238.0289,  237.0482,  244.0642, 243.0614,  247.0703,
		247.0703,   251.0796, 252.083,   257.0951,  258.10,    259.1009, 262.11,    2.0140,
		3.016045,   0.00000,  0.00000,   0.00000,   1.00794,   0.00000,  0.00000,   0.00000,
		0.00000,    0.00000,  0.00000,   0.00000,   0.00000,   0.00000,  0.00000,   0.00000,
		0.00000,    0.00000,  0.00000,   0.00000,   0.00000,   0.00000,  0.00000,   0.00000,
		0.00000,    0.00000,  0.00000,   0.00000,   0.00000,   0.00000
	};


const int possibleAromatic [NAROMMAX] = {7,8,15,16,33,34,51,52,HETERO_ATOM};
const int metals[NMETALS] = {
  3,4,11,12,13,19,20,21,22,23,24,25,26,27,28,29,
  30,31,37,38,39,40,41,42,43,44,45,46,47,48,49,50,55,56,57,58,59,60,61,62,63,
  64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,87,88,89,90,91,
  92,93,94,95,96,97,98,99,100,101,102,103};

//#define A B   ??? - 
const int lightMetals[NLIGHT_METALS] = {
  3,4,11,12,13,19,20,21,22,23,24,25,26,27,28,29,30,31,37,38};
const int heavyMetals[NHEAVY_METALS] = {
    39,40,41,42,43,44,45,46,47,48,49,50,55,56,57,58,59,60,61,62,63,
    64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,87,88,89,90,91,
    92,93,94,95,96,97,98,99,100,101,102,103};
const int halogens[NHALOGENS] = {9,17,35,53,85};
const int hetero[NHETERO] = {7,8,14,15,16,33,34,51,52,84};


const int NUMBEREXACT = 21;
const int exactAtom[NUMBEREXACT]={6,14,5,50,82,8,16,34,52,7,15,33,51,9,17,35,53,32,13,26,80};
const int alkaly[5]={3,11,19,37,55};
const int alkalyEarth[5]={4,12,20,38,56};
const int trivalent[31]={21,31,39,49,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,81,89,90,91,92,93,94,95,96,97,98,99};
const int titan[3]={22,40,72};
const int vanadium[3]={23,41,73};
const int cromium[3]={24,42,74};
const int manganeze[3]={25,43,75};
const int likeFe[2]={27,28};
const int platinum[6]={44,45,46,76,77,78};
const int copper[3]={29,47,79};
const int zink[2]={30,48};



const string fsastart="{SA:";
const string fsbstart="{SB:";

int getAtomPosition(const std::string atomicSymbol){
  int result=-1;
  int i;

  for (i=0; i<NELEMMCDL; i++) if (aSymb[i].size() == 3) if (atomicSymbol.find(aSymb[i]) == 0) {
	result=i;
	break;
  };
  for (i=0; i<NELEMMCDL; i++) if (aSymb[i].size() == 2) if (atomicSymbol.find(aSymb[i]) == 0) {
	result=i;
	break;
  };
  if (result < 0) for (i=0; i<NELEMMCDL; i++) if (aSymb[i].size() == 1) if (atomicSymbol.find(aSymb[i]) == 0) {
	result=i;
	break;
  };
  return result;
}

int encoder(int na) {

 //For given atom's number ATN in array ATOM returns a number from 1 to 32.
 //Atoms with related properties have the same number. It is required to reduce
 //the possible number of fragments in structure}
  int i;
  int result=32;

  for (i=0; i<NUMBEREXACT; i++) if (exactAtom[i] == na) {
    result=i;
    return result;
  };

  for (i=0; i<5; i++) if (na == alkaly[i]) {
	result=21; return result;
  };
  for (i=0; i<5; i++) if (na == alkalyEarth[i]) {
	result=22; return result;
  };
  for (i=0; i<31; i++) if (na == trivalent[i]) {
	result=23; return result;
  };
  for (i=0; i<3; i++) if (na == titan[i]) {
	result=24; return result;
  };
  for (i=0; i<3; i++) if (na == vanadium[i]) {
	result=25; return result;
  };
  for (i=0; i<3; i++) if (na == cromium[i]) {
	result=26; return result;
  };
  for (i=0; i<3; i++) if (na == manganeze[i]) {
	result=27; return result;
  };
  for (i=0; i<3; i++) if (na == copper[i]) {
	result=28; return result;
  };
  for (i=0; i<2; i++) if (na == likeFe[i]) {
	result=29; return result;
  };
  for (i=0; i<2; i++) if (na == zink[i]) {
	result=30; return result;
  };
  for (i=0; i<6; i++) if (na == platinum[i]) {
	result=31; return result;
  };
  return result;
};


bool isHetero(int na) {
  bool result=false;
  for (int i=0; i<NHETERO; i++) if (hetero[i] == na) {
	result=true;
	break;
  };
  return result;
}

bool isHalogen(int na) {
  bool result=false;
  for (int i=0; i<NHALOGENS; i++) if (halogens[i] == na) {
	result=true;
	break;
  };
  return result;
}


bool isMetall(int na){
  bool result=false;
  for (int i=0; i<NMETALS; i++) if (metals[i] == na) {
	result=true;
	break;
  };
  return result;
};


string intToStr(int k) {
	char temp[16];
	sprintf(temp,"%d",k);
	string line=temp;
	return line;
};

int hydrogenValency(int na) {   //Hydrogen valency
	int		result = 0;
	if(na > 0 && na < NELEMMCDL) result = hVal[na];
	return result;
};

int maxValency(int na) {        //Maximal valency of a specified element
	int		result = 8;
	if(na > 0 && na < NELEMMCDL) result = maxVal[na];
	return result;
};

std::string getAtomSymbol(int na) 
{
	if(na > 0 && na < NELEMMCDL) return aSymb[na];
	return string("");
};

double getAtomMass(int na) 
{
	if(na > 0 && na < NELEMMCDL) return aMass[na];
	return 0.0;
};

int compareStringsNumbers(string s1, string s2) {

  int result;
  int n,i;

  n=s1.length();
  if (s2.length()>n) n=s2.length();
  if (s1.length()<n) for (i=s1.length(); i<n; i++) {
    if ((s1.at(0) >= '0') && (s1.at(0) <= '9')) s1="0"+s1; else s1=s1+"0";
  };
  if (s2.length()<n) for (i=s2.length(); i<n; i++) {
    if ((s2.at(0) >= '0') && (s2.at(0) <= '9')) s2="0"+s2; else s2=s2+"0";
  };
  //lowest priority -zz

  result=s1.compare(s2);
  return result;
}

std::string getInChIKeyString(const fragmentCode & data){
  std::string result;

 // return "";
  result.assign(data.inChIKeyB);
  return result;
};  

void setInChIKeyString(const std::string & dataToSet, fragmentCode & data){
 // return;
  std::string s;

  memset(&data.inChIKeyB[0],0,INCHI_KEY_SIZE);
  s=dataToSet;
  if (s.size() >= INCHI_KEY_SIZE) {
	s=s.substr(0,INCHI_KEY_SIZE-1);
  };
  //if (dataToSet.size() >= INCHI_KEY_SIZE) throw;
  if (s.size() > 0) strcpy(&data.inChIKeyB[0],s.c_str());
};

std::string clearInChIKey(const std::string & inChI) {
  std::string result=inChI; 
  if (result.size() >= INCHI_KEY_SIZE) result=result.substr(0,INCHI_KEY_SIZE-1);
  return result;
};

int readSDFRecord(std::istream & data, std::vector<std::string> * recordData) {
	//return number of successfully read lines
	bool test = data.good();
	std::string s;
	int result=0;

	while (test) {
		getline(data,s);
		if (recordData) recordData->push_back(s);
		result++;
		if (s=="$$$$") test = false;
		if (test && (!data.good())) test = false;
	}
	return result;
}

std::string vectorToString(const std::string & delimiter, const std::vector<std::string> & data) {
	std::string result = "";
	int i;
	for (i = 0; i < data.size(); i++) result = result + data[i] + delimiter;
	return result;
}


/*
std::string intToStr(int n) {
  std::stringstream ss;
  ss << n;
  std::string result = ss.str();
  return result;
};
*/

void sortIntegers(int iLo, int iHi, std::vector<int> & data) {
	int lO, hI, mid, n;
	if ((iLo <0) || (iHi < 0)) return;
    lO = iLo;
    hI = iHi;
    mid = data[(lO + hI)/2];
	while (lO <= hI) {
		while (data[lO]<mid)  lO++;
		while(data[hI]>mid) hI--;
		if (lO <= hI) {			
			n = data[lO];
			data[lO] = data[hI];
			data[hI] = n;
			lO++;
			hI--;
		};
	}; // until Lo > Hi;
	if (hI>iLo) sortIntegers(iLo, hI, data);
	if (lO<iHi) sortIntegers(lO, iHi, data);
}

void quickSortIntegers(std::vector<int> & data) {
	int iLo, iHi;
	iLo = 0; iHi = data.size() - 1;
	sortIntegers(iLo, iHi, data);
}

void sortStrings(int iLo, int iHi, std::vector<std::string> & data, std::vector<int> * dataObject) {
	int lO, hI, n;
	std::string mid, s;

	if ((iLo <0) || (iHi < 0)) return;
	lO = iLo;
	hI = iHi;
	mid = data[(lO + hI) / 2];
	while (lO <= hI) {
		while (data[lO]<mid)  lO++;
		while (data[hI]>mid) hI--;
		if (lO <= hI) {
			s = data[lO];
			data[lO] = data[hI];
			data[hI] = s;
			if (dataObject) {
				n = (*dataObject)[lO];
				(*dataObject)[lO] = (*dataObject)[hI];
				(*dataObject)[hI] = n;
			};
			lO++;
			hI--;
		};
	}; // until Lo > Hi;
	if (hI>iLo) sortStrings(iLo, hI, data, dataObject);
	if (lO<iHi) sortStrings(lO, iHi, data, dataObject);
}

void quickSortStrings(std::vector<std::string> & data, std::vector<int> * dataObject) {
	int iLo, iHi;
	iLo = 0; iHi = data.size() - 1;
	sortStrings(iLo, iHi, data, dataObject);
}


int findQuick(const std::vector<int> & data, int count, int value) {
	int result = -1;
	int iHi, iLo, iMid, n;
	bool test;
	if (data.size() == 0) return result;
	iLo = 0; iHi = count-1;// data.size() - 1; non fullecontent of vector may be used foe quick search
	test = true;
	while (test) {
		iMid = (iLo + iHi) / 2;
		n = data[iMid];
		if (n > value) iHi = iMid - 1; else
		if (n < value) iLo = iMid + 1; else {
			result = iMid;
			test = false;
		}

		if (iLo > iHi) test = false;
	}

	return result;
};

int findQuickString(const std::vector<std::string> & data, int count, const std::string & value) {
	int result = -1;
	int iHi, iLo, iMid;
	std::string s;
	bool test;
	if (data.size() == 0) return result;
	iLo = 0; iHi = count - 1;// data.size() - 1; non fullecontent of vector may be used foe quick search
	test = true;
	while (test) {
		iMid = (iLo + iHi) / 2;
		s = data[iMid];
		if (s > value) iHi = iMid - 1; else
			if (s < value) iLo = iMid + 1; else {
				result = iMid;
				test = false;
			}

		if (iLo > iHi) test = false;
	}

	return result;
};

int findCodeIndex(const std::string & s, const std::vector<std::string> & codeList, const std::vector<std::string> & unsortedCodes) {
	int i;
	int result = findQuickString(codeList, codeList.size(), s);
	if (result < 0) for (i = 0; i < unsortedCodes.size(); i++) if (s == unsortedCodes[i]) {
		result = TEMP_SHIFT + i;
		break;
	};
	return result;
}


} // namespace MolStruct
