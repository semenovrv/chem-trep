
#include <vector>
#include <string>
#include <cmath>
#include <cstdio>

#include "single_atom.h"

using namespace std; 

namespace MolStruct {

#define NEXACTATOMS 21
const int exactAtom[NEXACTATOMS]={6,14,5,50,82,8,16,34,52,7,15,33,51,9,17,35,53,32,13,26,80};
#define NALKALYATOMS 5
const int alkaly[NALKALYATOMS] ={3,11,19,37,55};
#define NALKALYEARTHATOMS 5
const int alkalyEarth[NALKALYEARTHATOMS] ={4,12,20,38,56};
#define NTRIVALENTATOMS 31
const int trivalent[NTRIVALENTATOMS] ={21,31,39,49,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,81,89,90,91,92,93,94,95,96,97,98,99};
#define NTITANATOMS 3
const int titan[NTITANATOMS] ={22,40,72};
#define NVANADIUMATOMS 3
const int vanadium[NVANADIUMATOMS] ={23,41,73};
#define NCHROMIUMATOMS 3
const int cromium[NCHROMIUMATOMS] ={24,42,74};
#define NMANGANESEATOMS 3
const int manganeze[NMANGANESEATOMS] ={25,43,75};
#define NLIKEFEATOMS 2
const int likeFe[NLIKEFEATOMS] ={27,28};
#define NPLATINUMATOMS 6
const int platinum[NPLATINUMATOMS] ={44,45,46,76,77,78};
#define NCOPPERATOMS 3
const int copper[NCOPPERATOMS] ={29,47,79};
#define NZINKATOMS 2
const int zink[NZINKATOMS] ={30,48};



TSingleAtom * TSingleAtom::clone(void) const
{
	TSingleAtom * result;
	result=new TSingleAtom();

	result->astereo=this->astereo;
	for (int i=0; i<CONNMAX; i++) result->ac[i]=this->ac[i];
	result->currvalence=this->currvalence;
	result->iz=this->iz;
	result->na=this->na;
	result->nb=this->nb;
	result->nc=this->nc;
	result->nv=this->nv;
	result->rl=this->rl;
	result->rx=this->rx;
	result->ry=this->ry;
	result->special=this->special;
	result->enumerator=this->enumerator;
	result->fragIndex=this->fragIndex;
	result->anum=this->anum;
	if (atomList) {
	  result->atomList=new vector<int>();
	  result->atomList->resize(atomList->size());
	  for (int i=0; i<atomList->size(); i++) (*(result->atomList))[i]=(*atomList)[i];
	};
	for (int i=0; i<RMG_QUERY_SIZE; i++) result->rmgAtom[i]=rmgAtom[i];

   //molecular dynamic
	result->rejectInList=this->rejectInList;
   //result->processAddInfo=this->processAddInfo;
	result->hybridization=this->hybridization; 
	result->coordination = this->coordination;
	result->aromatic=this->aromatic; 
	result->chirality=this->chirality; 
	result->rKind1=this->rKind1; 
	result->rSize1=this->rSize1; 
	result->rKind2=this->rKind2; 
	result->rSize2=this->rSize2; 
	result->rKind3=this->rKind3; 
	result->rSize3=this->rSize3;    
	result->rKind4=this->rKind4;
	result->rSize4=this->rSize4;
	result->rKind5=this->rKind5;
	result->rSize5=this->rSize5;
	result->rKind6=this->rKind6;
	result->rSize6=this->rSize6;

	return result;
};


int TSingleAtom::chargeDeltaValency(int atomNo) 
{
	int result=-1;
	if (atomNo < NELEMMCDL) result=chargeVal[atomNo];
	return result;
};

int TSingleAtom::encoder() 
{
	// For given atom's number ATN in array ATOM returns a number from 1 to 32.
	//  Atoms with related properties have the same number. It is required to reduce
	//  the possible number of fragments in structure}
	int i;

	for (i=0; i<NEXACTATOMS; i++) if (exactAtom[i] == na) return i;
	for (i=0; i<NALKALYATOMS; i++) if (na == alkaly[i]) return 21;
	for (i=0; i<NALKALYEARTHATOMS; i++) if (na == alkalyEarth[i]) return 22;
	for (i=0; i<NTRIVALENTATOMS; i++) if (na == trivalent[i]) return 23;
	for (i=0; i<NTITANATOMS; i++) if (na == titan[i]) return 24;
	for (i=0; i<NVANADIUMATOMS; i++) if (na == vanadium[i]) return 25;
	for (i=0; i<NCHROMIUMATOMS; i++) if (na == cromium[i]) return 26;
	for (i=0; i<NMANGANESEATOMS; i++) if (na == manganeze[i]) return 27;
	for (i=0; i<NCOPPERATOMS; i++) if (na == copper[i]) return 28;
	for (i=0; i<NLIKEFEATOMS; i++) if (na == likeFe[i]) return 29;
	for (i=0; i<NZINKATOMS; i++) if (na == zink[i]) return 30;
	for (i=0; i<NPLATINUMATOMS; i++) if (na == platinum[i]) return 31;
	return 32;
};

int TSingleAtom::chargeConversion() 
{

	//  3 - if radical label present.
	//  2 - if charge <0
	//  1 - if charge >0
	//  0 - charge=0

	if (nc < 0) return 2;
	if (nc > 0) return 1;
	return 0;
};

int TSingleAtom::valencyConversion() 
{

	//for atom's number ATN in array ATOM returns some connected with valency value:
	// =2-Valence of ATN is less, then usual
	// =1-Valence of ATN is more, then usual
	// =0-usual valency


	int k1,k2;
	int result=0;

	//Default hydrogen
	k1=nv;
	k1=k1-currvalence-abs(nc)-rl;
	if (k1<0) k1=0;
	k2=hVal[na];
	k2=k2-currvalence-abs(nc)-rl;
	if (k2<0) k2=0;
	if (k1 == k2) result=0; else if (k1 < k2) result=1; else result=2;
	return result;
};

int TSingleAtom::allAtAtom()
{
	//Define a digital representation of atom, they include:
	//  a) Position of atom in the Periodic System
	//  b) Its charge
	//  c) Its valency}

	int b1,b2,b3,w;

	b1=encoder();
	b2=chargeConversion();
	b3=valencyConversion();
	w=b3;
	w=w << 2;
	w=w+b2;
	w=w << 5;
	w=w+b1;
	if (rl != 0) w= ~w;
	return w;
};

void TSingleAtom::addToList(int an){
  if (! atomList) atomList=new std::vector<int>();
  atomList->push_back(an);
};


void TSingleAtom::atomCopy(TSingleAtom * source)
{
	this->astereo=source->astereo;
	for (int i=0; i<CONNMAX; i++) this->ac[i]=source->ac[i];
	this->currvalence=source->currvalence;
	this->iz=source->iz;
	this->na=source->na;
	this->nb=source->nb;
	this->nc=source->nc;
	this->nv=source->nv;
	this->rl=source->rl;
	this->rx=source->rx;
	this->ry=source->ry;
	this->special=source->special;
	if (source->atomList) {
      if (! atomList) atomList=new std::vector<int>();
	  atomList->resize(source->atomList->size());
	  for (int i=0; i<source->atomList->size(); i++) (*atomList)[i]=(*(source->atomList))[i];
	} else {
	  delete atomList;
	  atomList=NULL;
	};
	for (int i=0; i<RMG_QUERY_SIZE; i++) rmgAtom[i]=source->rmgAtom[i];

   //molecular dynamice
	this->rejectInList=source->rejectInList;
   //this->processAddInfo=source->processAddInfo;
	this->hybridization=source->hybridization;
	this->coordination = source->coordination;
	this->aromatic=source->aromatic; 
	this->chirality=source->chirality; 
	this->rKind1=source->rKind1; 
	this->rSize1=source->rSize1; 
	this->rKind2=source->rKind2; 
	this->rSize2=source->rSize2; 
	this->rKind3=source->rKind3; 
	this->rSize3=source->rSize3; 
	this->rKind4=source->rKind4;
	this->rSize4=source->rSize4;
	this->rKind5=source->rKind5;
	this->rSize5=source->rSize5;
	this->rKind6=source->rKind6;
	this->rSize6=source->rSize6;
};

bool TSingleAtom::ringConditionsEquivalent(int queryKind, int querySize, int structureKind, int structureSize) {
  bool result;
  
  result=((queryKind == structureKind) && (querySize == structureSize));
  if ((! result) && (queryKind == structureKind)) result=(querySize == 1);
  
  return result;	
}

bool TSingleAtom::atomMolDynEquivalent(TSingleAtom * structure, TSingleAtom * query){
  bool result=true;
  bool firstRingChecked,secondRingChecked,thirdRingChecked,fourthRingChecked,fifthRingChecked,sixthRingChecked;
  bool test;

  if (query->hybridization > 0) result=(query->hybridization == structure->hybridization);
  if (result && (query->coordination > 0)) {
	  result = (query->coordination == structure->nb);
	  //if (result && (query->coordination==4)) {
	  	//	  test = true;
	  //};
  };

  if (result && (query->aromatic > 0)) result=(query->aromatic == structure->aromatic);
  // bad! if (result && (query->aromatic > 0)) result = ((structure->rSize1 == 5) || (structure->rSize1 == 6) || (structure->rSize2 == 5) || (structure->rSize2 == 6));
  if (result && (query->chirality > 0)) result=(query->chirality == structure->chirality);
  if (! result) return result;
  firstRingChecked=false;
  secondRingChecked=false;
  thirdRingChecked=false;
  fourthRingChecked=false;
  fifthRingChecked=false;;
  sixthRingChecked=false;
  if (result) if ((query->rKind1 > 0) || (query->rSize1 > 0)) {
    test=ringConditionsEquivalent(query->rKind1,query->rSize1,structure->rKind1,structure->rSize1);
    if (test) firstRingChecked=true; else {
      test=ringConditionsEquivalent(query->rKind1,query->rSize1,structure->rKind2,structure->rSize2);
	  if (test) secondRingChecked = true; else {
		  test = ringConditionsEquivalent(query->rKind1, query->rSize1, structure->rKind3, structure->rSize3);
		  if (test) thirdRingChecked = true; else {
			  test = ringConditionsEquivalent(query->rKind1, query->rSize1, structure->rKind4, structure->rSize4);
			  if (test) fourthRingChecked = true; else {
				  test = ringConditionsEquivalent(query->rKind1, query->rSize1, structure->rKind5, structure->rSize5);
				  if (test) fifthRingChecked = true; else {
					  test = ringConditionsEquivalent(query->rKind1, query->rSize1, structure->rKind6, structure->rSize6);
					  if (test) sixthRingChecked = true; else result = false;
				  };
			  };
		  };
	   };
    };
  };
  
  if (result) if ((query->rKind2 > 0) || (query->rSize2 > 0)) {
	  test = ringConditionsEquivalent(query->rKind2, query->rSize2, structure->rKind1, structure->rSize1);
	  if (test) firstRingChecked = true; else {
		  test = ringConditionsEquivalent(query->rKind2, query->rSize2, structure->rKind2, structure->rSize2);
		  if (test) secondRingChecked = true; else {
			  test = ringConditionsEquivalent(query->rKind2, query->rSize2, structure->rKind3, structure->rSize3);
			  if (test) thirdRingChecked = true; else {
				  test = ringConditionsEquivalent(query->rKind2, query->rSize2, structure->rKind4, structure->rSize4);
				  if (test) fourthRingChecked = true; else {
					  test = ringConditionsEquivalent(query->rKind2, query->rSize2, structure->rKind5, structure->rSize5);
					  if (test) fifthRingChecked = true; else {
						  test = ringConditionsEquivalent(query->rKind2, query->rSize2, structure->rKind6, structure->rSize6);
						  if (test) sixthRingChecked = true; else result = false;
					  };
				  };
			  };
		  };
	  };
  };
  
  if (result) if ((query->rKind3 > 0) || (query->rSize3 > 0)) {
	  test = ringConditionsEquivalent(query->rKind3, query->rSize3, structure->rKind1, structure->rSize1);
	  if (test) firstRingChecked = true; else {
		  test = ringConditionsEquivalent(query->rKind3, query->rSize3, structure->rKind2, structure->rSize2);
		  if (test) secondRingChecked = true; else {
			  test = ringConditionsEquivalent(query->rKind3, query->rSize3, structure->rKind3, structure->rSize3);
			  if (test) thirdRingChecked = true; else {
				  test = ringConditionsEquivalent(query->rKind3, query->rSize3, structure->rKind4, structure->rSize4);
				  if (test) fourthRingChecked = true; else {
					  test = ringConditionsEquivalent(query->rKind3, query->rSize3, structure->rKind5, structure->rSize5);
					  if (test) fifthRingChecked = true; else {
						  test = ringConditionsEquivalent(query->rKind3, query->rSize3, structure->rKind6, structure->rSize6);
						  if (test) sixthRingChecked = true; else result = false;
					  };
				  };
			  };
		  };
	  };
  };
  
  return	result;
}

bool TSingleAtom::atomEquivalent(TSingleAtom * structure, TSingleAtom * query,
	int nHStr, int nHQuery, bool chargeSensitivity, bool isotopeSensitivity,
	bool compareHard, bool rejectAnyInStructure, bool isUnsaturated, bool useAddInfo)
{
	bool result=false;
	int i,j,k;
	int naQuery;
	//bool test;
	TSingleAtom * queryCopy=NULL;

	if ((query->atomList) && (query->atomList->size() > 0)) {
	  queryCopy=query->clone();
	  delete queryCopy->atomList;
      queryCopy->atomList=NULL;
	  for (i=0; i<query->atomList->size(); i++) {
		queryCopy->na=(*(query->atomList))[i];
        result=atomEquivalent(structure,queryCopy,nHStr,nHQuery,chargeSensitivity,isotopeSensitivity,compareHard,rejectAnyInStructure,isUnsaturated,useAddInfo);
        if (query->rejectInList) result=(! result);
		if (result) break;
	  };
	  delete queryCopy;
	  if (result && useAddInfo) result = atomMolDynEquivalent(structure, query);
	  return result;
	};

	if (rejectAnyInStructure && (structure->na == ANY_ATOM)) {
		if (result && useAddInfo) result = atomMolDynEquivalent(structure, query);
		return result;
	};
	if ((structure == NULL) || (query == NULL)) {
		if (result && useAddInfo) result = atomMolDynEquivalent(structure, query);
		return result;
	};
	i=nHQuery;
	if (i>0) {   //if explicitly defined H is accociated with query atom
		j=structure->nv;
		k=abs(structure->nc);
		if (k>9) k=k-9;
		j=j-structure->currvalence-k;
		if (structure->rl != 0) j--;
		if (j<0) j=0;
		j=j+nHStr;
		if (i > j) {
			if (result && useAddInfo) result = atomMolDynEquivalent(structure, query);
			return result;
		};
		//on corresponding structure atom must be at least the same number of H}
	} else if ((query->na < NELEMMCDL) && (query->nv < hVal[query->na]) && ((query->special & NOOTHER_MASK) != 0))  {
		j=structure->nv;
		k=abs(structure->nc);
		if (k>9) k=k-9;
		j=j-structure->currvalence-k;
		if (structure->rl != 0) j--;
		if (j<0) j=0;
		if (j != 0) {
			if (result && useAddInfo) result = atomMolDynEquivalent(structure, query);
			return result;
		};
	}
	if (((structure->nc != query->nc) && chargeSensitivity)
		|| ((structure->rl != query->rl) && chargeSensitivity)
		|| ((structure->iz != query->iz) && isotopeSensitivity)) {
		if (result && useAddInfo) result = atomMolDynEquivalent(structure, query);
		return result;
	};
	//non-equivalency of Charge and Isotop is case of exact matching of these attributes

    if (structure->na == 255) return result;  //255 means NOT equal with any atom....
    if (compareHard) {
	  result = (structure->na == query->na) && (structure->special == query->special);
	  if (result && useAddInfo) result = atomMolDynEquivalent(structure, query);
      return result;
	};

  //Addition from January 2007 - if Structure atom is non-physical-then ONLY exact matching....
	if (structure->na > NPHYSICAL_ATOMS) if (structure->na != query->na) {
		if (result && useAddInfo) result = atomMolDynEquivalent(structure, query);
		return result;
	};

	i=query->special; //'NO OTHER' attribute on query
	if ((i & NOOTHER_MASK) != 0)  if ((structure->nb-nHStr) != query->nb) return result;
	//'NO OTHER' means, that number of neighbour in query must be the same as in structure
	if ((i & NONSATURATED_MASK) != 0) if (!isUnsaturated) {
		if (result && useAddInfo) result = atomMolDynEquivalent(structure, query);
		return result;
	};


	//Checking for equivalent position in the Periodic Table
	if (structure->na == query->na) result=true; else {

	  //test=true;
	  //nn=0;
	  //while (test  && (! result)) {
		//if ((query->atomList) && (query->atomList->size() > 0)) naQuery=(*(query->atomList))[nn]; else naQuery=query->na;
		naQuery=query->na;

		if (naQuery == METALL_ATOM) { //Metall checking
			for (i=0; i<NMETALS; i++) if (structure->na == metals[i]) {
				result=true;
				break;
			}
		}
		if (naQuery == HALOGEN_ATOM) { //Halogen checking
			for (i=0; i<NHALOGENS; i++) if (structure->na == halogens[i]) {
				result=true;
				break;
			}
		}
		if (naQuery == HETERO_ATOM) { //Hetero checking
			for (i=0; i<NHETERO; i++) if (structure->na == hetero[i]) {
				result=true;
				break;
			}
		}
		if ((naQuery == ANY_ATOM) || (naQuery == R_ATOM) || (naQuery == R_NOT_H_ATOM)) result=true; //Any atom
		if (naQuery > R_NOT_H_ATOM) for (i=0; i<RMG_QUERY_SIZE; i++) if (naQuery == structure->rmgAtom[i]) {
		  result=true;
		  break;
		};
	}

	if (result && useAddInfo) result=atomMolDynEquivalent(structure,query);
	return result;
}

void TSingleAtom::setRingProperty(int ringSize, bool isPlanar) {
  if (rSize1 == 0) {
    rSize1=ringSize;
	if (isPlanar) rKind1=1; else rKind1=2;
  } else if (rSize2 == 0) {
    rSize2=ringSize;
	if (isPlanar) rKind2=1; else rKind2=2;
  } else if (rSize3 == 0) {
	  rSize3 = ringSize;
	  if (isPlanar) rKind3 = 1; else rKind3 = 2;
  } else if (rSize4 == 0) {
	  rSize4 = ringSize;
	  if (isPlanar) rKind4 = 1; else rKind4 = 2;
  } else if (rSize5 == 0) {
	  rSize5 = ringSize;
	  if (isPlanar) rKind5 = 1; else rKind5 = 2;
  } else if (rSize6 == 0) {
    rSize6=ringSize;
	if (isPlanar) rKind6=1; else rKind6=2;
  };
};

TSingleAtom::TSingleAtom() {
	na = 6;
	nv = hVal[na];
	nc = 0;
	iz = 0;
	nb = 0;
	rl = 0;
	currvalence = 0;
	special = 0;
	astereo = 0;
	enumerator = 0;
	fragIndex = 0;
	rx=0;
	ry=0;
	anum="";
	for (int i=0; i<CONNMAX; i++) ac[i]=-1;
	for (int i=0; i<RMG_QUERY_SIZE; i++) rmgAtom[i]=0;
	atomList=NULL;
	
	//Molecular dynamic
	rejectInList=false;
  // processAddInfo=false;
	hybridization=0; 
	coordination = 0;
	aromatic=0; 
	chirality=0; 
	rKind1=0; 
	rSize1=0; 
	rKind2=0; 
	rSize2=0; 
	rKind3=0; 
	rSize3=0; 	
	rKind4=0;
	rSize4=0;
	rKind5=0;
	rSize5=0;
	rKind6=0;
	rSize6=0;
}

TSingleAtom::~TSingleAtom() {
  delete atomList;
};

ostream& operator<<(ostream& out, const TSingleAtom& sa)
{
	out << "Single atom:\n";
	out << sa.na << "  ";
	out << sa.nv << "  ";
	out << sa.nc << "  ";
	out << sa.iz << "\n";
	out << sa.rx << "  ";
	out << sa.ry << "  ";
	out << sa.rl << "  ";
	out << sa.nb << "  ";
	out << sa.currvalence << "  ";
	out << sa.special << " \n";

	for(int j = 0; j < CONNMAX; ++j)
		out << sa.ac[j] << " ";
	out << "\n";

	out << sa.astereo << "  ";
	out << sa.enumerator << "  ";
	out << sa.fragIndex << " \n";

	return out;
}


} // namespace MolStruct
