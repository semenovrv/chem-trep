
#include <vector>
#include <string>
#include <cmath>
#include <stdio.h>

#include "single_bond.h"

using namespace std;

namespace MolStruct {

const int bondValence[NBONDTYPES] = {1,2,3,1,1,0,0,0,1,1,1};

TSingleBond * TSingleBond::clone() const {
	TSingleBond * result;
	result=new TSingleBond();

	result->at[0]=this->at[0];
	result->at[1]=this->at[1];
	result->bstereo=this->bstereo;
	result->db=this->db;
	result->special=this->special;
	result->tb=this->tb;
	result->enumerator=this->enumerator;
	return result;
};

int TSingleBond::getValence() {
	int result=0;

	if ((this->tb <=NBONDTYPES) && (this->tb>0)) result=bondValence[this->tb-1];
	return result;
};

bool TSingleBond::bondEquivalent(TSingleBond * sBond, TSingleBond * qBond, bool compareHard){
	//for bond's number BONDNUMBER in structure (BOND) and for query bond's number
	//QBONDNUMBER (QBOND) determines, if they may be accociated one with other.
	//QBC-array, for each query bond contains explicitly defined cyclic conditions:
	// 0-no cyclic condition is implemented
	// 1-query bond and corresponding structure bond must be a Chain
	// 2-query bond and corresponding structure bond must be a Ring
	//Function has 'TRUE' value if bonds can be matched, FALSE otherwise

	bool test1;
	int bT,qBT,bD,qBD;
	bool result=false;

    if (compareHard) {
		result=(sBond->tb == qBond->tb) && (sBond->db == qBond->db) && (sBond->special == qBond->special);
        return result;
	};

    //Addition from January 2007 - if Structure bond is non-physical-then ONLY exact matching....
	test1=(sBond->tb == 1) || (sBond->tb == 2) || (sBond->tb == 3) || (sBond->tb == 4) || (sBond->tb == 9) || (sBond->tb == 10) || (sBond->tb == 11);
	if (! test1) if (sBond->tb != qBond->tb) return result;


	bT=sBond->tb; qBT=qBond->tb;
	bD=sBond->db; qBD=qBond->db;
	if ((bT >= 9) && (bT <= 11)) bT=1;
	if ((qBT >= 9) && (qBT <= 11)) qBT=1;

	if (sBond->special>=3) sBond->special=sBond->special-3;
	if (qBond->special>=3) qBond->special=qBond->special-3;
	if (((qBond->special==1) && (bD>1)) || ((qBond->special==2) && (bD<2))) return result;
	//Explicitly defined cyclic conditions are violated
	if ((qBD>1) && (bD<2)) return result;
	//Implicitly defined cyclic conditions are violated}
	if (qBT==ANY_BOND) {  //Any bond in query
		result=true;
		return result;
	};
	if ((qBT==AROMATIC_BOND) || (qBD==2) || (qBD==3)) { //Aromatic bond
		result=(bD==2) || (bD==3) || (bT==AROMATIC_BOND);
		return result;
	};
	if (checkAromatic) if ((bD==2) || (bD==3)) {
		//non-Aromatic in query, but Aromatic in structure
		result=false;
		return result;
	};
	if ((qBT==SINGLE_OR_DOUBLE_BOND) && ((bT==1) || (bT==2) || (bT==SINGLE_OR_DOUBLE_BOND))) return true;
	//SINGLE/DOUBLE query bond
	if (bT==qBT) return true;
	//Type of bonds in structure and in query are equivalent}

	//New bond types - for RMG
	if ((qBT==SINGLE_OR_TRIPLE_BOND) && ((bT==1) || (bT==3) || (bT==SINGLE_OR_TRIPLE_BOND))) return true;
	if ((qBT==UNSATURATED_BOND) && ((bT==2) || (bT==3) || (bT==UNSATURATED_BOND))) return true;
	return false;
};


int TSingleBond::bondConversion() {
	//generate a code for bond's nuber BNB in array BOND. It is took into consi-
	//deration bond type and cycle size
	int b1,b2;

	if (tb > 4) b1=0; else b1=tb;
	b2=7; //for other Switch value
	switch (db) {
 case 0 :
	 b2=5;
	 break;
 case 1 :
	 b2=0;
	 break;
 case 2 :
	 b1=4;
	 b2=1;
	 break;
 case 3 :
	 b1=4;
	 b2=2;
	 break;
 case 4 :
	 b2=3;
	 break;
 case 5 :
	 b2=4;
	 break;
 case 6 :
	 b2=6;
	 break;
	};
	b2=b2 << 2;
	b2=b2+b1;
	return b2;
}

void TSingleBond::bondCopy(TSingleBond * source) 
{
	this->at[0]   = source->at[0];
	this->at[1]   = source->at[1];
	this->bstereo = source->bstereo;
	this->db      = source->db;
	this->special = source->special;
	this->tb      = source->tb;
}

TSingleBond::TSingleBond(int atom0, int atom1, int type)
{
	at[0] = atom0;
	at[1] = atom1;
	tb = type;
	db = 0;    
	bstereo = 0;
	special = 0;
	enumerator = 0;
}

TSingleBond::TSingleBond(void)
{
	tb = 1;
	db = 0;    
	bstereo = 0;
	special = 0;
	at[0] = -1;
	at[1] = -1;
	enumerator = 0;
}

ostream& operator<<(ostream& out, const TSingleBond& sb)
{
	out << "Single bond:\n";

	out << sb.tb << "   ";
	out << sb.at[0] << "   " << sb.at[1] << "\n";
	out << sb.db << "   ";
	out << sb.bstereo << "   ";
	out << sb.special << "   ";
	out << sb.enumerator << "\n";

	return out;
}


} // namespace MolStruct


