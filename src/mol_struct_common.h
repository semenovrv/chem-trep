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

#ifndef _HDR_MOL_STRUCT_COMMON_
#define _HDR_MOL_STRUCT_COMMON_

#include <string>

namespace MolStruct {

#define RUNDEF -1.2345678E9
#define INCHI_KEY_SIZE 28

//common constants
//static const int MAXBONDS=999;
static const int MAXFRAGS=200;
static const int MAXCHARS=1000;
static const int MAX_DEPTH=10;
static const int NELEMMAX=120;
static const int TEMP_SHIFT = 2000000000;

#define NELEMMCDL 135
#define OBMCDL 

#ifndef BUFF_SIZE
#define BUFF_SIZE 32768
#endif

//#define CONNMAX 15
#define CONNMAX 20
//#define NATOMSMAX 999
//#define NBONDSMAX 999
//#define NQUERYMAX 64
//#define NATOMSMAXSMBIG 2000
//#define NBONDSMAXSMBIG 2000



extern const int hVal[NELEMMCDL];
extern const int maxVal[NELEMMCDL];
extern const int chargeVal[NELEMMCDL];
extern const std::string aSymb[NELEMMCDL];
extern const double aMass[NELEMMCDL];
//extern const double vdDist[NELEMMCDL];
//extern const double kappa[NELEMMCDL];
//extern const int nElectrons[NELEMMCDL];

#define NMETALS 78
#define NHALOGENS 5
#define NHETERO 10
#define NAROMMAX 9
#define NLIGHT_METALS 20
#define NHEAVY_METALS 58
#define NPHYSICAL_ATOMS 110

#define DEUTERIUM_ATOM 103
#define TRITIUM_ATOM 104
#define EXPLICIT_HYDROGEN 109
#define LP_ATOM 110

#define METALL_ATOM 111
#define HALOGEN_ATOM 112
#define ANY_ATOM 113
#define HETERO_ATOM 114
#define ID_ZVEZDA 115

//RMG atoms
#define R_NOT_H_ATOM 121
#define CBF_ATOM 122
#define CDD_ATOM 123
#define SIS_ATOM 124
#define SID_ATOM 125
#define CS_ATOM 126
#define CD_ATOM 127
#define CT_ATOM 128
#define CB_ATOM 129
#define CO_ATOM 130
#define OS_ATOM 131
#define OD_ATOM 132
#define SS_ATOM 133
#define R_ATOM 134

/*
  ID_METALLATOM=211;
  ID_HALOGENATOM=212;
  ID_ANYATOM=213;
  ID_HETEROATOM=214;

*/

#define ANY_BOND 8
#define SINGLE_OR_TRIPLE_BOND 0
#define UNSATURATED_BOND 7
#define AROMATIC_BOND 4
#define SINGLE_OR_DOUBLE_BOND 5

#define NOOTHER_MASK 1
#define AROMATIC_MASK 2
#define EXACTNUMBER_MASK 4
#define LIST_MASK 8
#define NOTLIST_MASK 16
#define NONSATURATED_MASK 32



extern const int possibleAromatic [NAROMMAX];
extern const int metals[NMETALS];

extern const int lightMetals[NLIGHT_METALS];
extern const int heavyMetals[NHEAVY_METALS];

extern const int halogens[NHALOGENS];
extern const int hetero[NHETERO];


typedef struct adjustedlist{
  int nb;
  int adjusted[CONNMAX];
} adjustedlist;

typedef struct fragmentCode{
  int unc1, unc2, unc3;
  char inChIKeyB[INCHI_KEY_SIZE];
} fragmentCode;


typedef struct twoSphereRecord {
 int c1,c2;
} twoSphereRecord;


//typedef adjustedlist neigbourlist[NATOMSMAXSMBIG];



std::string		getAtomSymbol(int na);
double			getAtomMass(int na);
// Return valency by hydrogen for given atomic position in the Periodic Table
int				hydrogenValency(int na);
int				maxValency(int na);
bool			isHetero(int na);
bool			isHalogen(int na);
bool			isMetall(int na);
int             getAtomPosition(const std::string atomcSymbol);

std::string getInChIKeyString(const fragmentCode & data);  
void setInChIKeyString(const std::string & dataToSet, fragmentCode & data);
std::string clearInChIKey(const std::string & inChI);
void quickSortIntegers(std::vector<int> & data);
int findQuick(const std::vector<int> & data, int count, int value);
void quickSortStrings(std::vector<std::string> & data, std::vector<int> * dataObject=NULL);
int findQuickString(const std::vector<std::string> & data, int count, const std::string & value);
int findCodeIndex(const std::string & s, const std::vector<std::string> & codeList, const std::vector<std::string> & unsortedCodes);
int readSDFRecord(std::istream & data, std::vector<std::string> * recordData);
std::string vectorToString(const std::string & delimiter, const std::vector<std::string> & data);


struct cmpFragmentCode {
    bool operator()(const fragmentCode & a, const fragmentCode & b) const {
	  bool result=false;
	  if (a.unc1 < b.unc1) result=true; else if (a.unc1 == b.unc1) {
		if (a.unc2 < b.unc2) result=true; else if (a.unc2 == b.unc2) {
		  if (a.unc3 < b.unc3) result=true; else if (a.unc3 == b.unc3) {
			std::string sA=getInChIKeyString(a);
			std::string sB=getInChIKeyString(b);
			if (sA < sB) result=true;
		  };
		};
	  };
      return result;
    }
};


std::string intToStr(int k);

int compareStringsNumbers(std::string s1, std::string s2);
int encoder(int na);



class Rect {
  public:
    double left,top,right,bottom;
};

class Point {
  public:
    double x,y;
};

class cfIOPT {
public:
	cfIOPT(void)
	{
		fIOPT1=2;
		fIOPT2=true;
		fIOPT3=1;
		fIOPT4=false;
		fIOPT5=2;
		fIOPT7=false;
		fIOPT8=true;
		fIOPT9=true;
		fIOPT10=true;
		fIOPT11=true;
		fIOPT12=3;
		fIOPT13=true;
		fIOPT15=12;
		fIOPT16=10;
	};
  int			fIOPT1;
    /*Hydrogen Show =1 - Off =2 -On =3 - Hetero
      Determines whether or not implicitly defined hydrogens atoms
      should be shown on a screen*/
  bool	   fIOPT2;
    /*Stereo bond   =1 - On  =2 -Off
      Determines whether or not stereo bonds (BT=9,10,11-see RBOND
      type definition) should be shown as stereo or as a single.*/
  int			fIOPT3;
    /*Stereo subsriptors (R/S/Z/E) =1 - Symbol (draw string 'R','S'..
                                   =2 - Color
                                   =3 - Off (not show)*/
  bool		fIOPT4;
    //{Atom number Show =1 - Off, =2 - On}
  int       fIOPT5;
    /*Bond spacing = 1 - Small, =2 - Normal, =3 - Large. Determines
      line distance in double, triple and related bonds.*/
  bool    fIOPT7;
    /*Clean:Atoms Shift =1-On, -2-Off. Determines whether or not
      small shifts for atom positions is necessary while CLEAN
      procedure is executed, to avoid atom's overlapping.*/
  bool    fIOPT8;
   /* AutoFitScreen =1-On, =2-Off. Determines whether or not the
      structure should be rescaled, if some elements of structure
      exceed boundaries of Wn[5].*/
  bool    fIOPT9;
   /*Flip (command) =1-Global, =2-Local. If Global was installed,
     Flip command enable while structure drawing, otherwise only through icon.*/
  bool    fIOPT10;
   /*Charge Sensitivity =1-Yes, =2-No. Determines whether or not
     atoms, which differ in the only charge, should be different
     or identical while Structure/Substructure search.*/
  bool    fIOPT11;
   /*Isotop Sensitivity =1-Yes, =2-No. The same as for CHARGE.*/
  int       fIOPT12;
   /*Stereo search =1-Ignore, =2-Exact, =3-Replace. Ignore-all
     stereo bonds are treated as single. Replace-Up and Down
     stereobonds in structure may be replaced.*/
  bool    fIOPT13;
   /*Semipolar bond (A+---B-). =1-Unique, =2-Equal Double. Type
     of representation of semipolar bond while Structure/Substructure search.*/
  int       fIOPT15;
   /*Number of points on ring while COMPASS draw.*/
  int       fIOPT16;
   /*CleanBond ratio 1:IOPT[15]. It is used for Clean Bond command
    (see CLEANATOMS procedure SED5.PAS file for detailed)*/

};

}


#endif