

#define _USE_MATH_DEFINES

#include <vector>
#include <string>
#include <cstdio>
#include <cstring>
#include <iomanip>
#include <sstream>
#include <cmath> 
//#include "inchi.h"

#include <fstream>
#include <iostream>


#include "simple_molecule.h"

using namespace std;

#ifndef WIN32
#define _snprintf			snprintf
#endif

namespace MolStruct {

#define blDenominator 4.0   //Controls bond legth in bondEnlarge
#define nRotBondsMax 20     //Determines no. rotating bonds in correctOverlapped


void deleteIntElement(std::vector<int> * source, int index) 
{
	int i,n;
	std::vector<int> temp(source->size()-1);

	n=0;
	for (i=0; i<source->size(); i++) 
		if (i != index) {
			temp[n]=(*source)[i];
			n++;
		};
	source->resize(source->size()-1);
	for (i=0; i<source->size(); i++) (*source)[i]=temp[i];
};

std::string removeSpacesTab(const std::string data) {
  std::string result=data;
  int n;

  while ((result.size() > 0) && ((result.at(0) == ' ') || (result.at(0) == '\t'))) result=result.substr(1);
  while ((result.size() > 0) && ((result.at(result.size()-1) == ' ')  || (result.at(result.size()-1) == '\t'))) result=result.substr(0,result.size()-1);
  n=result.find('\t');
  while (n != std::string::npos) {
    result=result.substr(0,n)+' '+result.substr(n+1);
    n=result.find('\t');
  };
  return result;
};

bool compareAtoms(int a1, int a2, const std::vector<std::vector<int> *> aeqList) 
{
	std::vector<int> * l1;
	std::vector<int> * l2;
	int i;
	bool result=false;

	if ((a1 < 0) || (a2 < 0) || (a1 >= aeqList.size()) || (a2 >= aeqList.size())) return result;
	l1 = (std::vector<int> *)(aeqList.at(a1));
	l2 = (std::vector<int> *)(aeqList.at(a2));
	if ((l1 == NULL) || (l2 == NULL)) return result;
	if (l1->size() != l2->size()) return result;
	result=true;
	for (i=0; i<l1->size(); i++) if ((*l1)[i] != (*l2)[i]) {
		result=false;
		break;
	};
	return result;
};

bool incrementValues(std::vector<int>& currentValues, const std::vector<int> maxValues)
{
	int i,l;
	bool result=false;

	for (i=0; i<currentValues.size(); i++) {
		l=currentValues[i];
		l++;
		if (l <= maxValues[i]) {
			currentValues[i]=l;
			result=true;
			return result;
		} else currentValues[i]=0;
	};
	return result;
};

void writeSDF(const std::string & fieldName, const std::string & data, std::vector<std::string> & fileData) {
	if (data.size() == 0) return;
	fileData.push_back(">  <" + fieldName + ">");
	fileData.push_back(data);
	fileData.push_back("");
}


//sMol.pas
void packRingData(std::vector<std::vector<int> *> & ringData) {
	std::vector<std::vector<int> *> tempData;
	tempData.reserve(ringData.size());
	for (int i = 0; i < ringData.size(); i++) if (ringData[i] != NULL) tempData.push_back(ringData[i]);
	ringData.clear();
	for (int i = 0; i < tempData.size(); i++) ringData.push_back(tempData[i]);
};


bool listPresent(const std::vector<std::vector<int> *> & ringData, const std::vector<int> & tempList) {

	std::vector<int> *lR;
	int i, j;
	bool test;
	bool result = false;
	if (ringData.size() == 0) return result;
	for (i = 0; i < ringData.size(); i++) {
		lR = ringData[i];
		if ((lR->size() == tempList.size()) && (lR->size() > 0)) {
			test = true;
			for (j = 0; j < lR->size(); j++) if ((*lR)[j] != tempList[j]) {
				test = false;
				break;
			};
			result = test;
		};
		if (result) break;
	};
	return result;
};

bool bondPresent(std::vector<int> & bondList, int bn) {
	int i;
	int result = false;
	for (i = 0; i < bondList.size(); i++)if (bondList[i] == bn) {
		result = true;
		break;
	};
	return result;
};


double xDistPoint(double x1, double y1, double x2, double y2, double x0, double y0)
{
	//The function finds distance between FAtom AN and bond BN (including sign).
	// The bond is treated as segment, not as straight line!}
	double d, r, r1, yMin, yMax, xMin, xMax, xx, yy, result;

	if (y1 < y2) {
		yMin=y1;
		yMax=y2;
	} else {
		yMin=y2;
		yMax=y1;
	};
	xx=x1-x2;
	yy=y1-y2;
	r1=sqrt(xx*xx+yy*yy);
	yMin=yMin-0.1*r1; yMax=yMax+0.1*r1;
	d=y2-y1;
	if (abs(d) < 1E-8) {
		result=1E9;
		return result;
	};
	if ((y0 > yMin) && (y0 < yMax)) {
		r=x1+(y0-y1)*(x2-x1)/d;
		if (x1 > x2) {
			xMin=x2;
			xMax=x1;
		} else {
			xMin=x1;
			xMax=x2;
		};
		xMin=xMin-0.1*r1;
		xMax=xMax+0.1*r1;
		if (r < xMin) r=xMin;
		if (r > xMax) r=xMax;
		result=r-x0;
	} else result=1E9;
	return result;
};

bool overlapped(double x1A, double y1A, double x2A, double y2A,
					 double x1B, double y1B, double x2B, double y2B, double delta)
{

	double a1, b1, c1, a2, b2, c2, r, cX, cY, x, y, r1, r2;
	double xMin, xMax, yMin, yMax;
	bool result=false;

	r=y2A-y1A;
	if (abs(r) > 1E-9) {
		a1=1/r;
		cY=-y1A/r;
	} else {
		a1=1E9;
		cY=-y1A*1E9;
		if (r < 0) {
			a1=-a1;
			cY=-cY;
		};
	};
	r=x2A-x1A;
	if (abs(r) > 1E-9) {
		b1=1/r;
		cX=x1A/r;
	} else {
		b1=1E9;
		cX=x1A*1E9;
		if (r < 0) {
			b1=-b1;
			cX=-cX;
		};
	};
	b1=-b1;
	c1=cX+cY;
	r=y2B-y1B;
	if (abs(r) > 1E-9) {
		a2=1/r;
		cY=-y1B/r;
	} else {
		a2=1E9;
		cY=-y1B*1E9;
		if (r < 0) {
			a2=-a2;
			cY=-cY;
		};
	};
	r=x2B-x1B;
	if (abs(r) > 1E-9) {
		b2=1/r;
		cX=x1B/r;
	} else {
		b2=1E9;
		cX=x1B*1E9;
		if (r < 0) {
			b2=-b2;
			cX=-cX;
		};
	};
	b2=-b2;
	c2=cX+cY;
	r1=b1*c2-b2*c1;
	r2=a1*b2-a2*b1;
	if (abs(r2) > 1E-9) y=r1/r2; else {
		y=1E9;
		if (r1 < 0) y=-y;
	};
	r1=c1*a2-c2*a1;
	r2=a1*b2-a2*b1;
	if (abs(r2) > 1E-9) x=r1/r2; else {
		x=1E9;
		if (r1 < 0) x=-x;
	};
	if (x1A < x2A) {
		xMin=x1A;
		xMax=x2A;
	} else {
		xMin=x2A;
		xMax=x1A;
	};
	if (y1A < y2A) {
		yMin=y1A;
		yMax=y2A;
	} else {
		yMin=y2A;
		yMax=y1A;
	};
	xMin=xMin-delta;
	xMax=xMax+delta;
	yMin=yMin-delta;
	yMax=yMax+delta;
	result=((x >= xMin) && (x <= xMax) && (y >= yMin) && (y <= yMax));
	if (result) {
		if (x1B < x2B) {
			xMin=x1B;
			xMax=x2B;
		} else {
			xMin=x2B;
			xMax=x1B;
		};
		if (y1B < y2B) {
			yMin=y1B;
			yMax=y2B;
		} else {
			yMin=y2B;
			yMax=y1B;
		};
		xMin=xMin-delta;
		xMax=xMax+delta;
		yMin=yMin-delta;
		yMax=yMax+delta;
		result=((x >= xMin) && (x <= xMax) && (y >= yMin) && (y <= yMax));
	};
	if (! result) result=(abs(xDistPoint(x1A,y1A,x2A,y2A,x1B,y1B)) < delta);
	if (! result) result=(abs(xDistPoint(x1A,y1A,x2A,y2A,x2B,y2B)) < delta);
	if (! result) result=(abs(xDistPoint(x1B,y1B,x2B,y2B,x1A,y1A)) < delta);
	if (! result) result=(abs(xDistPoint(x1B,y1B,x2B,y2B,x2A,y2A)) < delta);
	return result;
};


//******************************************************************
//                      TSimpleMolecule
//------------------------------------------------------------------

void TSimpleMolecule::addAllHydrogens(bool forStereo) {
//It is assumed, that DefineConn were called early...

  int i,j,n,k;
  bool test;

  
  n=this->nAtoms();
  if (n > 0) for (j=0; j<nAtoms(); j++) {
    test=false;
	while (! test) {
	  k=getAtom(j)->nv;
	  k=k-getAtom(j)->currvalence-abs(getAtom(j)->nc)-getAtom(j)->rl;
      test=(k <= 0);
      if (! test) addDefaultHydrogen(j,forStereo);
	};
  };
};

bool TSimpleMolecule::isCoordinatesBad(int atn) {
  int i;
  double r,r1,x,y;
  bool result=false;

  r=averageBondLength();
  if (r == 0) return result;
  for (i=0; i<nAtoms(); i++) if (i != atn) {
	x=getAtom(i)->rx-getAtom(atn)->rx;
	y=getAtom(i)->ry-getAtom(atn)->ry;
	r1=sqrt(x*x+y*y);
    if (r1/r < 0.1) result=true;
    if (result) break;
  };
};


bool TSimpleMolecule::addDefaultHydrogen(int atn, bool forStereo){
  
  bool result=false;
  TSingleAtom * sa;
  TSingleBond * sb;
  double x,y,r,xa,ya;
  double xc,yc,xn,yn,xo,yo,r1;
  int i,n;
  bool testCoincided;
  int cHA,cHB;

  if (/*(nAtoms() >= NATOMSMAXSMBIG) || */ (atn < 0) || (atn >= nAtoms()) || (getAtom(atn)->nb >= CONNMAX)) return result;
  result=true;
  
  sa=new TSingleAtom();
  sa->na=1;
  sa->nv=1;
  unitVector(atn,x,y);
  if (forStereo) for (i=0; i<nBonds(); i++) {
    sb=getBond(i);
	if (((sb->tb == 9) || (sb->tb == 10)) && (sb->at[0] == atn)) {
	  xa=getAtom(sb->at[1])->rx-getAtom(atn)->rx;
	  ya=getAtom(sb->at[1])->ry-getAtom(atn)->ry;
      r=sqrt(xa*xa+ya*ya);
      if (r > 0) {
        x=(getAtom(sb->at[1])->rx-getAtom(atn)->rx)/r;
		y=(getAtom(sb->at[1])->ry-getAtom(atn)->ry)/r;
	  };
	};
  };
  if (getAtom(atn)->nb > 0) {
	n=getAtom(atn)->ac[0];
	xa=getAtom(atn)->rx-getAtom(n)->rx;
	ya=getAtom(atn)->ry-getAtom(n)->ry;
    r=sqrt(xa*xa+ya*ya);
  } else r=averageBondLength();
  if (r <= 0) r=30;
  sa->rx=getAtom(atn)->rx+x*r;
  sa->ry=getAtom(atn)->ry+y*r;
  sa->nb=1;
  sa->ac[0]=atn;
  sa->currvalence=1;
  addAtom(sa);
  getAtom(atn)->nb++;
  getAtom(atn)->ac[getAtom(atn)->nb-1]=nAtoms()-1;
  getAtom(atn)->currvalence++;
  sb=new TSingleBond();
  sb->at[0]=atn;
  sb->at[1]=nAtoms()-1;
  addBond(sb);
  //Testing if atoms are coincided
  testCoincided=isCoordinatesBad(nAtoms()-1);
  if (testCoincided && (getAtom(atn)->nb == 2)) {
    //Try to flip
    cHA=nAtoms()-1;
    cHB=atn;
    xc=getAtom(cHA)->rx-getAtom(cHB)->rx; //Axes direction calculation
    yc=getAtom(cHA)->ry-getAtom(cHB)->ry;
    r1=sqrt(xc*xc+yc*yc);
    xc=xc/r1;
    yc=yc/r1;
    xo=xc*xc-yc*yc; //Sqr(Xc)-Sqr(Yc);
    yo=2*xc*yc;
    //rotation at angle Pi around axes for all atoms in the LIST fragment
    xc=getAtom(cHA)->rx-getAtom(cHB)->rx;
    yc=getAtom(cHA)->ry-getAtom(cHB)->ry;
    xn=xc*xo+yc*yo;
    yn=xc*yo-yc*xo;
    getAtom(cHA)->rx=getAtom(cHB)->rx+xn;
    getAtom(cHA)->ry=getAtom(cHB)->ry+yn;
    testCoincided=isCoordinatesBad(nAtoms()-1);
  };
  if (testCoincided) {
    //Rescale bond length
    r=r/2;
    getAtom(nAtoms()-1)->rx=getAtom(atn)->rx+x*r;
    getAtom(nAtoms()-1)->ry=getAtom(atn)->ry+y*r;
  };
  return result;
};


int TSimpleMolecule::getExplicitH(int atomNo) const {
  int result=0;
  int i,n;

  for (i=0; i<getAtom(atomNo)->nb; i++) {
	n=getAtom(atomNo)->ac[i];
	if (getAtom(n)->na == 1) result++;
  };
  return result;
};


int TSimpleMolecule::getNH(int atomNo) const {
	int result=0;
	const TSingleAtom * atom;
	int i,n;

	if (atomNo >= nAtoms()) return result;
	atom=getAtom(atomNo);
	result=atom->nv;
	result=result-(atom->currvalence)+(atom->nc*TSingleAtom::chargeDeltaValency(atom->na))-(atom->rl);
	if (result < 0) result=0;
	if (atom->nb > 0) for (i=0; i<atom->nb; i++) {
		n=atom->ac[i];
		if (getAtom(n)->na == 1) result=result+1;
	};
	return result;
};

bool TSimpleMolecule::bondConnectsCH2(int bondNo){
  bool result=false;
  int an;

  if (getBond(bondNo)->tb != 1) return result;
  an=getBond(bondNo)->at[0];
  if ((getAtom(an)->na != 6) || (getAtom(an)->nv != 4) || (getAtom(an)->rl != 0)) return result;
  an=getBond(bondNo)->at[1];
  if ((getAtom(an)->na != 6) || (getAtom(an)->nv != 4) || (getAtom(an)->rl != 0)) return result;
  if ((getNH(getBond(bondNo)->at[0]) != 2) || (getNH(getBond(bondNo)->at[1]) != 2)) return result;
  return true;
};

void TSimpleMolecule::getHomolog(TSimpleMolecule * sMol){
  int i,j,k,n;

  sMol->clear();
  sMol->moleculeCopy(*this);
  sMol->defineAtomConn();
  sMol->allAboutCycles();

  //Acyclic
  for (i=sMol->nBonds()-1; i>=0; i--) if ((sMol->getBond(i)->db == 0) && (sMol->bondConnectsCH2(i))) {
	n=sMol->getBond(i)->at[0];
	k=-1;
	for (j=0; j<sMol->nBonds(); j++) if ((i != j) && ((sMol->getBond(j)->at[0] == n) || (sMol->getBond(j)->at[1] == n))) {
	  if ((sMol->getAtom(sMol->getBond(j)->at[0])->na == 6) && (sMol->getAtom(sMol->getBond(j)->at[0])->na == 6)) {
        k=j;  
	    break;
	  };
	};
	if (k >= 0) {
      if (sMol->getBond(k)->at[0] == n) sMol->getBond(k)->at[0]=sMol->getBond(i)->at[1]; else
      if (sMol->getBond(k)->at[1] == n) sMol->getBond(k)->at[1]=sMol->getBond(i)->at[1];
      sMol->deleteBond(i);	
      sMol->defineAtomConn();
	};
  };
  //Cyclic
  for (i=sMol->nBonds()-1; i>=0; i--) if ((sMol->getBond(i)->db > 4) && (sMol->bondConnectsCH2(i))) {
	n=sMol->getBond(i)->at[0];
	k=-1;
	for (j=0; j<sMol->nBonds(); j++) if ((i != j) && ((sMol->getBond(j)->at[0] == n) || (sMol->getBond(j)->at[1] == n))) {
	  if ((sMol->getAtom(sMol->getBond(j)->at[0])->na == 6) && (sMol->getAtom(sMol->getBond(j)->at[0])->na == 6)) {
        k=j;  
	    break;
	  };
	};
	if (k >= 0) {
      if (sMol->getBond(k)->at[0] == n) sMol->getBond(k)->at[0]=sMol->getBond(i)->at[1]; else
      if (sMol->getBond(k)->at[1] == n) sMol->getBond(k)->at[1]=sMol->getBond(i)->at[1];
      sMol->deleteBond(i);	
      sMol->defineAtomConn();
      sMol->allAboutCycles();
	};
  };
  //hydrogens in CH4
  for (i=sMol->nBonds()-1; i>=0; i--) if ((sMol->getAtom(sMol->getBond(i)->at[0])->na == 1) || (sMol->getAtom(sMol->getBond(i)->at[0])->na == 1)) 
  if ((sMol->getAtom(sMol->getBond(i)->at[0])->na == 6) || (sMol->getAtom(sMol->getBond(i)->at[0])->na == 6)) {
    sMol->deleteBond(i);	
    sMol->defineAtomConn();    
  };
  //Delete unconnected atoms
  sMol->defineAtomConn();    
  for (i=sMol->nAtoms()-1; i>=0; i--) if (sMol->getAtom(i)->nb == 0) sMol->deleteAtom(i);
  sMol->defineAtomConn();
  sMol->allAboutCycles();
};

bool TSimpleMolecule::isHomolog(TSimpleMolecule * other){
  bool result=false;
  TSimpleMolecule sm;
  int c1,c2,c3;
  std::string inChIKey;

  this->getHomolog(&sm);
  sm.defineAtomConn();
  sm.allAboutCycles();
  sm.calculateAllIndeces();

  c1=sm.unC1;
  c2=sm.unC2;
  c3=sm.unC3;
  inChIKey=sm.inChIKey;

  other->getHomolog(&sm);
  sm.defineAtomConn();
  sm.allAboutCycles();
  sm.calculateAllIndeces();

  result=((c1 == sm.unC1) && (c2 == sm.unC2) && (c3 == sm.unC3) && (inChIKey.compare(sm.inChIKey) == 0));
  return result;
};


int TSimpleMolecule::singleAtomicDescriptor(int aNumber,int bNumber, bool useEnumerator) {
	//Bond BNumber has to have type=9 or 10
	int result=0;
	int an [4];
	int j,m,k;
	double x[3];
	double y[3];
	double r,r1,r2;
	bool testBad;
	double tn [2];
	bool rsn;
	bool isInvert;

	if ((getAtom(aNumber)->nb < 3) || (getAtom(aNumber)->nb > 4)) return result;
	for (j=0; j<getAtom(aNumber)->nb; j++) an[j]=getAtom(aNumber)->ac[j];
	if (useEnumerator) {
		//fragIndex as main base
		for (j=0; j<(getAtom(aNumber)->nb-1); j++) for (m=j+1; m<getAtom(aNumber)->nb; m++) 
			if (getAtom(an[j])->fragIndex > getAtom(an[m])->fragIndex) {
				k=an[j];
				an[j]=an[m];
				an[m]=k;
			};
	} else {
		for (j=0; j<(getAtom(aNumber)->nb-1); j++) for (m=j+1; m<getAtom(aNumber)->nb; m++) 
			if (an[j] > an[m]) {
				k=an[j];
				an[j]=an[m];
				an[m]=k;
			};
	};

	//  return 1;

	for (j=0; j<3; j++) {
		k=an[j];
		x[j]=getAtom(k)->rx-getAtom(aNumber)->rx;
		y[j]=getAtom(k)->ry-getAtom(aNumber)->ry;
	};

	isInvert=true;
	if (getAtom(aNumber)->nb == 4) if (getBond(bNumber)->at[1] != an[3]) {
		isInvert=false;
		for (j=0; j<3; j++) if (an[j] == getBond(bNumber)->at[1]) {
			k=an[3];
			x[j]=getAtom(k)->rx-getAtom(aNumber)->rx;
			y[j]=getAtom(k)->ry-getAtom(aNumber)->ry;
		};
	};
	testBad=false;
	for (j=0; j<3; j++) {
		r=sqrt(x[j]*x[j]+y[j]*y[j]);
		if (r == 0) testBad=true; else {
			x[j]=x[j]/r;
			y[j]=y[j]/r;
		};
	};

	//  return 1;

	if (! testBad) {
		for (j=0; j<2; j++) { //determination of the sign of rotation
			r1=x[0]*y[j+1]-y[0]*x[j+1];     //Sin
			r2=x[0]*x[j+1]+y[0]*y[j+1];     //Cos
			if (r1 < 0) r2=-2-r2;           // 0 degress = +1; Pi/2 degrees = 0; Pi degrees = -1; 3*Pi/2 degrees = -2; 2*Pi degrees = -3
			tn[j]=r2;
		};
		rsn=( tn[0] > tn[1]);
		//return 1;

		//Analizing cases....
		if (getBond(bNumber)->tb == 10) rsn=! rsn;
		if (getAtom(aNumber)->nb == 4) if (isInvert) rsn=! rsn;
		if (tn[0] != tn[1]) {
			if (rsn) result=1; else result=2;
		};
	};
	return result;
};


void TSimpleMolecule::bondUnitVector(int bn, double& xv, double& yv) const {
	//for bon's number BN in the structure, described by ATOM, BOND, CONN calculates
	//the unit vector (Xv,Yv on output). The vector is best direction to add new fra-
	//gment to current structure. The procedure is used by TEMPLATE and MAKEPOLI}

	double si,r1,r2,r3,s1,s2,s3,s4,x1,y1,y2;
	int na1,na2,i;

	na1=getBond(bn)->at[0];
	na2=getBond(bn)->at[1];
	s1=getAtom(na1)->rx;
	s2=getAtom(na1)->ry;
	s3=getAtom(na2)->rx;
	s4=getAtom(na2)->ry;
	r1=s1-s3;
	r2=s2-s4;
	r3=sqrt(r1*r1+r2*r2);
	r1=r1/r3; //COS
	r2=r2/r3; //SIN
	si=0;
	for (i=0; i<getAtom(na1)->nb; i++) if (getAtom(na1)->ac[i] != na2) {
		x1=getAtom(getAtom(na1)->ac[i])->rx-s1;
		y1=getAtom(getAtom(na1)->ac[i])->ry-s2;
		y2=x1*r2-y1*r1;
		if (y2 != 0) si=si+y2/abs(y2);
	}
	for (i=0; i<getAtom(na2)->nb; i++) if (getAtom(na2)->ac[i] != na1) {
		x1=getAtom(getAtom(na2)->ac[i])->rx-s3;
		y1=getAtom(getAtom(na2)->ac[i])->ry-s4;
		y2=x1*r2-y1*r1;
		if (y2 != 0) si=si+y2/abs(y2);
	}
	if (si != 0) si=si/abs(si); else si=1;
	xv=-r2*si; yv=r1*si;
}


bool TSimpleMolecule::threeBondResolve(int an, int bondExcluded, double& xv, double& yv, std::vector<adjustedlist> * bkExt) const {
	//Addition from 16 April 2006
	int bondNoList [3];
	bool result=false;
	double centerX[3];
	double centerY[3];
	int nBondNo;
	std::vector<adjustedlist> * bk=NULL;
	std::vector<int>bondList(listarSize());
	std::vector<int>* blStore=NULL;
	bool testBad,testOK;
	int m,j,n1,n2,n,k;
	int rs;
	double dist,x,y,minDist;
	int i;
	bool test;
	int at;
	double r1,r2,s1,s3;

	if (bkExt != NULL) bk=bkExt; else {
		//bk=(neigbourlist *)malloc((CONNMAX+1)*/*NBONDSMAXSMBIG*/this->nAtoms()*4);
		bk = new std::vector<adjustedlist>;
		bk->resize(listarSize());
		defineBondConn(*bk);
	};
	nBondNo=0;
	testBad=false;
	testOK=true;

	for (i=0; i<(*bk)[an].nb; i++) {
		n=(*bk)[an].adjusted[i];
		if (n != bondExcluded) {
			vaweBond(n,*bk,rs,bondList);
			if (rs > 0) {
				//I have to analize bondList to determine second cycle to exclude adamanthane
				if (blStore == NULL) {
					//Save
					blStore= new std::vector<int>(rs);
					for (j=0; j<blStore->size(); j++) (*blStore)[j]=bondList[j];
				} else {
					m=0;
					for (j=0; j<blStore->size(); j++) {
						test=false;
						for (k=0; k<rs; k++) if (bondList[k] == (*blStore)[j]) {
							test=true;
							break;
						};
						if (test) m++;
					};
					if ((m > 1) && (m < blStore->size())) {
						testOK=false;
					};
				};
				//center determination
				centerX[nBondNo]=0;
				centerY[nBondNo]=0;
				for (j=0; j<rs; j++) {
					m=bondList[j];
					n1=getBond(m)->at[0];
					n2=getBond(m)->at[1];
					centerX[nBondNo]=centerX[nBondNo]+getAtom(n1)->rx+getAtom(n2)->rx;
					centerY[nBondNo]=centerY[nBondNo]+getAtom(n1)->ry+getAtom(n2)->ry;
				};
				centerX[nBondNo]=centerX[nBondNo]/(2*rs);
				centerY[nBondNo]=centerY[nBondNo]/(2*rs);
				bondNoList[nBondNo]=n;
				nBondNo++;
			};
		};// else testBad:=true;
		if ((nBondNo == 3) || (testBad)) break;
	};
	if (nBondNo < 2) testBad=true;

	if ((! testBad) && testOK) {
		dist=0;
		if (nBondNo == 2) {
			//single bond, attached to ring. Re-definition center....
			x=0; y=0;
			for (i=0; i<nAtoms(); i++) {
				x=x+getAtom(i)->rx;
				y=y+getAtom(i)->ry;
			};
			x=x/nAtoms();
			y=y/nAtoms();
			for (i=0; i<nBondNo; i++) {
				centerX[i]=x;
				centerY[i]=y;
			};
		};
		for (i=0; i<nBondNo; i++) {
			at=getBond(bondNoList[i])->at[0];
			if (at == an) at=getBond(bondNoList[i])->at[1];
			r1=getAtom(at)->rx-getAtom(an)->rx;
			r2=getAtom(at)->ry-getAtom(an)->ry;
			s1=sqrt(r1*r1+r2*r2);
			if (s1 == 0) testBad=true; else {
				r1=r1/s1; r2=r2/s1;
				x=-r1;
				y=-r2;
				minDist=1000000000;
				for (j=0; j<nBondNo; j++) {
					r1=(x*s1+getAtom(an)->rx-centerX[j]);
					r2=(y*s1+getAtom(an)->ry-centerY[j]);
					s3=sqrt(r1*r1+r2*r2);
					if (s3 < minDist) minDist=s3;
				};
				if (minDist > dist) {
					dist=minDist;
					xv=x;
					yv=y;
				};
			};
		};
		if (! testBad) result=true;
	};
	//end addition from 16 April 2006
	if (bkExt == NULL) delete(bk);
	if (blStore != NULL) delete(blStore);
	return result;
};

bool TSimpleMolecule::unitVectorCoincident(int aN, double xV, double yV) const {
	bool result=false;
	int i, aT;
	double r1,r2,s1;

	for (i=0; i<getAtom(aN)->nb; i++) {
		aT=getAtom(aN)->ac[i];
		r1=getAtom(aT)->rx-getAtom(aN)->rx;
		r2=getAtom(aT)->ry-getAtom(aN)->ry;
		s1=sqrt(r1*r1+r2*r2);
		if (s1 > 0.00001) {
			r1=r1/s1;
			r2=r2/s1;
			if ((abs(r1-xV) < 0.1) && (abs(r2-yV) < 0.1)) result=true;
		};
		if (result) break;
	};
	return result;
};

void TSimpleMolecule::unitVector(int aN, double& xV, double& yV) const {
	//For atom's number AN in the array ATOM (the same structure is described with
	//bond's array BOND and bond-connection matrix invariants CONN) unit vector
	//is calculated (XV, YV) on output. Unit vector shows the best direction to
	//make new bond}
	double sc[4] = {0,1/2,1.7320508/2,1}; 
	double cc[4] = {1,1.7320508/2,1/2,0};  
	double sQ3=sqrt(3.0)/2.0;
	double si,r1,r2,r3,r4,s1,s2,s3,s4,fi;
	double xm[3];  //Initial dimensions 1..3
	double ym[3];  //Initial dimensions 1..3
	int i,i1,i2;
	int nB1,aT,aT1;
	bool mK[CONNMAX];
	bool test;

	nB1=1;
	if (getAtom(aN)->nb == 0) {
		xV=1;
		yV=0;               //if no atoms connected-horizontal unit vector
		return;
	} else if (getAtom(aN)->nb == 1) {  //single atom connected-start calculations}
		aT=getAtom(aN)->ac[0];
		r1=getAtom(aT)->rx-getAtom(aN)->rx;
		r2=getAtom(aT)->ry-getAtom(aN)->ry;
		s1=sqrt(r1*r1+r2*r2);
		if (s1 == 0) {
			xV=0.3826834;
			yV=-0.9238795;
			return;
		}
		r1=r1/s1; r2=r2/s1;
		nB1=getAtom(aT)->nb;
		si=1;
		if (nB1 > 1) {
			aT1=getAtom(aT)->ac[0];
			if (aT1 == aN) aT1=getAtom(aT)->ac[1];
			r3=getAtom(aT1)->rx-getAtom(aN)->rx;
			r4=getAtom(aT1)->ry-getAtom(aN)->ry;
			s2=sqrt(r3*r3+r4*r4);
			if (s2 == 0) s2=1;
			r3=r3/s2; r4=r4/s2;
			si=r2*r3-r4*r1;
			if (si != 0) si=si/abs(si); else si=1;
		}
		s3=sqrt(3.0)/2.0;
		s4=-0.5;
		xV=r1*s4-r2*si*s3;  //single atom connected-unit vector at 2Pi/3 angle
		yV=r2*s4+r1*si*s3;
	} else {            //multiply connection case
		xV=0;
		yV=0;
		for (i=0; i<getAtom(aN)->nb; i++) mK[i]=true;
		if (getAtom(aN)->nb == 3) {
			//three-bond connection-search for projec. of tetrahedron

			if (threeBondResolve(aN,-1,xV,yV,NULL)) {
				return;
			};

			for (i=0; i<3; i++) {
				aT=getAtom(aN)->ac[i];
				xm[i]=getAtom(aT)->rx-getAtom(aN)->rx;
				ym[i]=getAtom(aT)->ry-getAtom(aN)->ry;
				r1=sqrt(xm[i]*xm[i]+ym[i]*ym[i]);
				xm[i]=xm[i]/r1; ym[i]=ym[i]/r1;
			}
			i=0;
			do {      //Search for bond, which make 30 degree angle with another
				i1=1; i2=2;
				switch (i) {
		 case 0: {
			 i1=1; i2=2; break;
					}
		 case 1: {
			 i1=0; i2=2; break;
					}
		 case 2: {  //Else value - only three
			 i1=0; i2=1; break;
					}
				}
				s1=xm[i1]*xm[i2]+ym[i1]*ym[i2];
				s2=xm[i]*xm[i1]+ym[i]*ym[i1];
				s3=xm[i]*xm[i2]+ym[i]*ym[i2];
				test=false;
				if (abs(s1+0.5) < 0.05) {  //two bonds form 120 degrees angle
					test=(abs(s2-sQ3)<0.05) || (abs(s3-sQ3)<0.05);
					if (! test) test=((abs(s2) < 0.05) && (abs(s3+sQ3)<0.05))
						|| ((abs(s3)<0.05) && (abs(s2+sQ3)<0.05));
				}
				i++;
			} while (! (test || (i == 3)));
			if (test) mK[i-1]=false;  //In the original Java code without -1
		}
		for (i=0; i<getAtom(aN)->nb; i++) if (mK[i]) {
			aT=getAtom(aN)->ac[i];
			r1=getAtom(aT)->rx-getAtom(aN)->rx;
			r2=getAtom(aT)->ry-getAtom(aN)->ry;
			s1=sqrt(r1*r1+r2*r2);
			if (s1 < 0.05) s1=1;
			r1=r1/s1; r2=r2/s1;
			xV=xV-r1;
			yV=yV-r2;
		}
		r1=sqrt(xV*xV+yV*yV);
		if (r1 > 0.05) {
			//unit vector may be calculated from current connection
			xV=xV/r1;
			yV=yV/r1;
		} else {
			//impossible to calculate unit vector-more definitions
			if (getAtom(aN)->nb == 2) {    //linear existing connection-Pi/2 angle
				aT=getAtom(aN)->ac[0];
				r1=getAtom(aT)->rx-getAtom(aN)->rx;
				r2=getAtom(aT)->ry-getAtom(aN)->ry;
				s1=sqrt(r1*r1+r2*r2);
				if (s1 < 0.05) {
					xV=0.3826834;
					yV=-0.9238795;
					return;
				}
				r1=r1/s1; r2=r2/s1;
				xV=-r2; yV=-r1;
			} else if (getAtom(aN)->nb == 3) { //3-bonds connection - projection of tetrahedrone is calculated
				//Which of bonds has the horizontal direction
				r3=100000;
				for (i=0; i<3; i++) {
					aT=getAtom(aN)->ac[i];
					r1=getAtom(aT)->rx-getAtom(aN)->rx;
					r2=getAtom(aT)->ry-getAtom(aN)->ry;
					s1=sqrt(r1*r1+r2*r2);
					if (s1 < 0.05) {
						xV=0.3826834;
						yV=-0.9238795;
						return;
					}
					r1=abs(r1/s1); r2=abs(r2/s1);
					if (r1 < r3) {r3=r1; nB1=i;}
					if (r2 < r3) {r3=r2; nB1=i;}
				}
				aT=getAtom(aN)->ac[nB1];       //Horizontal or vertical bond found}
				r1=getAtom(aT)->rx-getAtom(aN)->rx;
				r2=getAtom(aT)->ry-getAtom(aN)->ry;
				s1=sqrt(r1*r1+r2*r2);
				r1=r1/s1; r2=r2/s1;
				xV= r1*cos(150*M_PI/180)+r2*sin(150*M_PI/180);
				yV=-r1*sin(150*M_PI/180)+r2*cos(150*M_PI/180);
			} else {
				xV=sqrt(2.0)/2.0; //other case - Pi/4 to horizont angle
				yV=-xV;
			}
		}

		//checking if coincided xV and uV with existing bonds
		if (unitVectorCoincident(aN,xV,yV)) {
			fi=15*M_PI/180.0;
			r3=cos(fi);
			r4=sin(fi);
			r1= xV*r3+yV*r4;
			r2=-xV*r4+yV*r3;
			if (unitVectorCoincident(aN,r1,r2)) {
				r1=xV*r3-yV*r4;
				r2=xV*r4+yV*r3;
			};
			if (unitVectorCoincident(aN,r1,r2)) {
				fi=7.5*M_PI/180.0;
				r3=cos(fi);
				r4=sin(fi);
				r1= xV*r3+yV*r4;
				r2=-xV*r4+yV*r3;
				if (unitVectorCoincident(aN,r1,r2)) {
					r1=xV*r3-yV*r4;
					r2=xV*r4+yV*r3;
				};
			};
			xV=r1;
			yV=r2;
		};
	};
	//Correction of angles, closed to 0, 30, 60, 90 degrees to exact values
	if (getAtom(aN)->nb < 6) for (i=0; i<4; i++) if ((abs(abs(xV)-sc[i])<0.04)
		&& (abs(abs(yV)-cc[i]) < 0.04)) {
			if (xV < 0) xV=-sc[i]; else xV=sc[i];
			if (yV < 0) yV=-cc[i]; else yV=cc[i];
	}
}

void TSimpleMolecule::clear() {
  TSingleAtom * sa;
  TSingleBond * sb;
  for (int i=0; i<fAtom.size(); i++) {
    sa=fAtom[i];
	fAtom[i]=NULL;
    delete sa;
  };
  fAtom.clear();

  for (int i=0; i<fBond.size(); i++) {
	sb=fBond[i];
	fBond[i]=NULL;
	delete sb;
  };
  fBond.clear();
};

void TSimpleMolecule::assignBondFromStereodescriptor(int atomNo, int atomChirality, bool ignoreH) {
  int i,j,bnDouble,n1,n2,n,k;
  bool test;
  double r,max,x1,y1,x2,y2;
  int bc[4];

  if (atomChirality == 0) return;
  if ((getAtom(atomNo)->nb < 3) || (getAtom(atomNo)->nb > 4)) return;  //Two connection - cannot contain atomic descriptos...
  bnDouble=0;
  n=0;
  for (i=0; i<nBonds(); i++) if ((getBond(i)->at[0] == atomNo) || (getBond(i)->at[1] == atomNo)) {  //Is a stereobond, pointed to Atom No exists-immediate return
    if (((getBond(i)->tb == 9) || (getBond(i)->tb == 10)) && (getBond(i)->at[0] == atomNo)) return;
    if (getBond(i)->tb == 2) bnDouble=i;
    bc[n]=i;
    n++;
  };
  if (bnDouble != 0) {        //Double bond-can be on S or P only!!!
	test=((getAtom(atomNo)->na == 15) || (getAtom(atomNo)->na == 16));
    test=(test && (getAtom(atomNo)->nb == 4));
    if (! test) return;
  };
  //Search for pair of bonds with smallest angle...
  if (getAtom(atomNo)->nb == 3) {
    //Search fo acyclic bond if possible...
    n=-1;
	for (i=0; i<getAtom(atomNo)->nb; i++) if (getBond(bc[i])->tb == 1) if (getBond(bc[i])->db <= 2) {
      n=i;
      break;
	};
    //Search for any single bond if cannot found acyclie
	if (n == -1) for (i=0; i<getAtom(atomNo)->nb; i++) if (getBond(bc[i])->tb == 1) {
      n=i;
      break;
	};
  } else {
    //Search for bond pairs with smallest angles....
    n1=-1; n2=-1; max=-1; n=-1;
    for (i=0; i<(getAtom(atomNo)->nb-1); i++) for (j=i+1; j<getAtom(atomNo)->nb; j++) {
      k=bc[i];
	  if (getBond(k)->at[0] == atomNo) k=getBond(k)->at[1]; else k=getBond(k)->at[0];
	  x1=getAtom(k)->rx-getAtom(atomNo)->rx;
      y1=getAtom(k)->ry-getAtom(atomNo)->ry;
      r=sqrt(x1*x1+y1*y1);
      if (r != 0) {
        x1=x1/r; y1=y1/r;
	  };
      k=bc[j];
      if (getBond(k)->at[0] == atomNo) k=getBond(k)->at[1]; else k=getBond(k)->at[0];
      x2=getAtom(k)->rx-getAtom(atomNo)->rx;
      y2=getAtom(k)->ry-getAtom(atomNo)->ry;
      r=sqrt(x2*x2+y2*y2);
      if (r != 0) {
        x2=x2/r; y2=y2/r;
	  };
      r=x1*x2+y1*y2;
      if (r > max) {
        n1=i; n2=j; max=r;
	  };
	};
    if ((n1 >= 0) && (n2 >= 0)) {
      if ((getBond(bc[n1])->tb == 1) && (getBond(bc[n1])->db < 2)) n=n1; else
      if ((getBond(bc[n2])->tb == 1) && (getBond(bc[n2])->db < 2)) n=n2; else
      if (getBond(bc[n1])->tb == 1) n=n1; else
      if (getBond(bc[n2])->tb == 1) n=n2; else
      if (getBond(bc[n1])->tb == 2) n=n1; else
      if (getBond(bc[n2])->tb == 2) n=n2;
	};
  };
  if (n >= 0) {
    n=bc[n];
    if (getBond(n)->at[0] == atomNo) k=getBond(n)->at[1]; else k=getBond(n)->at[0];
    getBond(n)->at[0]=atomNo;
    getBond(n)->at[1]=k;
    if (getBond(n)->tb == 2) {  //Make semipolar bond from double
      getBond(n)->tb=1;
	  getAtom(atomNo)->nc=getAtom(atomNo)->nc+1;
	  getAtom(k)->nc=getAtom(k)->nc-1;
	};
    getBond(n)->tb=9;
    k=singleAtomicDescriptor(atomNo,n,false);
    if (k != atomChirality) getBond(n)->tb=10;
    k=singleAtomicDescriptor(atomNo,n,false);
    if (k != atomChirality) getBond(n)->tb=1;
  };
};


void TSimpleMolecule::setCoordinatesString(const std::string value){
	//set back above prcedure
	TSingleAtom * sa;
	TSingleBond * sb;
	int i,n,k,kk;
	std::string s;

	clear();
	n=0;
	s=value.substr(n,2); //2 symbols - number of atoms
	n=n+2;
	kk=strtol(s.c_str(),NULL,10);
	for (i=0; i<kk; i++) {
		sa=new TSingleAtom();
		s=value.substr(n,1);
		if (s == "1") sa->na=101; else
			if (s == "2") sa->na=102; else
				if (s == "3") sa->na=103; else
					if (s == "A") sa->na=ANY_ATOM; else
						if (s == "E") {
							sa->na=ANY_ATOM;
							sa->special=1;
						} else sa->na=6;//Atom.positionofAtom(s);
		n=n+1;
		s=value.substr(n,4);
		n=n+4;
		k=strtol(s.c_str(),NULL,10);
		sa->rx=(double)k/10000.0;
		s=value.substr(n,4);
		n=n+4;
		k=strtol(s.c_str(),NULL,10);
		sa->ry=(double)k/10000.0;
		addAtom(sa);
	};
	s=value.substr(n,2);  //2 symbols - number of bonds
	n=n+2;
	kk=strtol(s.c_str(),NULL,10);
	for (i=0; i<kk; i++) {
		sb=new TSingleBond();
		sb->tb=ANY_BOND;  //any bond
		s=value.substr(n,2);
		n=n+2;
		k=strtol(s.c_str(),NULL,10);
		sb->at[0]=k-1;
		s=value.substr(n,2);
		n=n+2;
		k=strtol(s.c_str(),NULL,10);
		sb->at[1]=k-1;
		addBond(sb);
	};
	defineAtomConn();
	allAboutCycles();
};

void TSimpleMolecule::setCoordinatesStringEx(const std::string value){
	//set back above prcedure
	TSingleAtom * sa;
	TSingleBond sb;
	int i,n,k,kk,an;
	std::string s;

	clear();
	n=0;
	s=value.substr(n,3); //2 symbols - number of atoms
	n=n+3;
	kk=strtol(s.c_str(),NULL,10);
	for (i=0; i<kk; i++) {
		sa=new TSingleAtom();
		//memset(sa,0,sizeof(TSingleAtom));

		s=value.substr(n,2);  //s=value.substr(n,1);
		an=strtol(s.c_str(),NULL,10);
		sa->na=an;//Atom.positionofAtom(s);
		n=n+2;  //n+1;
		s=value.substr(n,1);
		if (s=="1") sa->nc=1; else
		if (s=="2") sa->nc=2; else
		if (s=="3") sa->nc=3; else
		if (s=="4") sa->nc=4; else
		if (s=="5") sa->nc=-1; else
		if (s=="6") sa->nc=-2; else
		if (s=="7") sa->nc=-3; else
   		if (s=="8") sa->nc=-4; else
		if (s=="9") sa->rl=1;
        n=n+1;  //!!! New 
		s=value.substr(n,4);
		n=n+4;
		k=strtol(s.c_str(),NULL,10);
		sa->rx=(double)k/10000.0;
		s=value.substr(n,4);
		n=n+4;
		k=strtol(s.c_str(),NULL,10);
		sa->ry=(double)k/10000.0;
		sa->nv=hVal[sa->na];
		addAtom(sa);
	};
	s=value.substr(n,3);  //2 symbols - number of bonds
	n=n+3;
	kk=strtol(s.c_str(),NULL,10);
	for (i=0; i<kk; i++) {
		//sb=new TSingleBond();
		//memset(sb,0,sizeof(TSingleBond));
		s=value.substr(n,1);
		n=n+1;  //New
		k=strtol(s.c_str(),NULL,10)+1;
		if (k == 8) k=11;
		sb.tb=k;  //any bond
		s=value.substr(n,3);
		n=n+3;
		k=strtol(s.c_str(),NULL,10);
		sb.at[0]=k-1;
		s=value.substr(n,3);
		n=n+3;
		k=strtol(s.c_str(),NULL,10);
		sb.at[1]=k-1;
		addBond(sb.tb,sb.at[0],sb.at[1]);
	};
	defineAtomConn();
	allAboutCycles();
};

int TSimpleMolecule::setResourceString(const std::string value){
	//set back above prcedure
	TSingleAtom * sa;
	TSingleBond sb;
	int i,n,k,kk,an,nc,rl;
	std::string s;

	clear();
	n=0;
	s=value.substr(n,2); //3 symbols - number of atoms
	n=n+2;
	kk=strtol(s.c_str(),NULL,16);
	for (i=0; i<kk; i++) {
	  sa=new TSingleAtom();
	  //memset(sa,0,sizeof(TSingleAtom));

	  s=value.substr(n,2);  //s=value.substr(n,1);
	  an=strtol(s.c_str(),NULL,16);
	  if (an > 200) an=an-100; //Delphi has other definition of Halogen, Metal, Any, etc.
	  sa->na=an;
	  n=n+2; 

	  s=value.substr(n,2);
	  nc=strtol(s.c_str(),NULL,16);
	  sa->nc=nc-127;
	  n=n+2;

	  s=value.substr(n,1);
	  rl=strtol(s.c_str(),NULL,16);
	  sa->rl=rl;
	  n=n+1;

  	  s=value.substr(n,4);
	  k=strtol(s.c_str(),NULL,16);
	  sa->rx=(double)k/65500.0;
	  n=n+4;

	  s=value.substr(n,4);
	  k=strtol(s.c_str(),NULL,16);
	  sa->ry=(double)k/10000.0;
	  n=n+4;
		
	  sa->nv=hVal[sa->na];
	  addAtom(sa);
	};
	s=value.substr(n,2);  //3 symbols - number of bonds
	n=n+2;
	kk=strtol(s.c_str(),NULL,16);
	for (i=0; i<kk; i++) {
  	  s=value.substr(n,1);
	  k=strtol(s.c_str(),NULL,16);
	  sb.tb=k;  //any bond
	  n=n+1;  //New

	  s=value.substr(n,2);
	  k=strtol(s.c_str(),NULL,16);
	  sb.at[0]=k-1;  //Numerations of atoms is started from 1
	  n=n+2;

	  s=value.substr(n,2);
	  k=strtol(s.c_str(),NULL,16);
	  sb.at[1]=k-1;
	  n=n+2;
	  addBond(sb.tb,sb.at[0],sb.at[1]);
	};
	defineAtomConn();
	allAboutCycles();
	return n;
};

int TSimpleMolecule::setJMEString(const std::string value){
	TSingleAtom * sa;
	TSingleBond sb;
	int n,i,na,nb,nc,an;
	std::string s,data;
	double rx,ry;

	clear();
	data=value;

	n=data.find(' ');
	s=data.substr(0,n);
	data=data.substr(n+1);
	na=strtol(s.c_str(),NULL,10);

	n=data.find(' ');
	s=data.substr(0,n);
	data=data.substr(n+1);
	nb=strtol(s.c_str(),NULL,10);

	for (i=0; i<na; i++) {
	  n=data.find(' ');
	  s=data.substr(0,n);
	  data=data.substr(n+1);
	  
	  nc=0;
	  if (s.at(n-1) == '+') {
		nc=1; 
		s=s.substr(0,n-1);
	  } else if (s.at(n-1) == '-') {
        nc=-1;
		s=s.substr(0,n-1);
	  };
	  an=getAtomPosition(s);

	  n=data.find(' ');
	  s=data.substr(0,n);
	  data=data.substr(n+1);
	  rx=strtod(s.c_str(),NULL);

	  n=data.find(' ');
	  if (n == std::string::npos) {
        s=data;
		data="";
	  } else {
	    s=data.substr(0,n);
	    data=data.substr(n+1);
	  };
	  ry=strtod(s.c_str(),NULL);

	  sa=new TSingleAtom();
	  //memset(sa,0,sizeof(TSingleAtom));

	  sa->na=an;
	  sa->nc=nc;
	  sa->rx=rx;
	  sa->ry=ry;
	  sa->nv=hVal[sa->na];
	  addAtom(sa);
	};


	for (i=0; i<nb; i++) {

	  n=data.find(' ');
	  s=data.substr(0,n);
	  data=data.substr(n+1);
      n=strtol(s.c_str(),NULL,10);
	  sb.at[0]=n-1;

	  n=data.find(' ');
	  s=data.substr(0,n);
	  data=data.substr(n+1);
      n=strtol(s.c_str(),NULL,10);
	  sb.at[1]=n-1;

	  n=data.find(' ');
	  if (n == std::string::npos) {
        s=data;
		data="";
	  } else {
	    s=data.substr(0,n);
	    data=data.substr(n+1);
	  };
      n=strtol(s.c_str(),NULL,10);
	  sb.tb=n;
      addBond(sb.tb,sb.at[0],sb.at[1]);
	};


	defineAtomConn();
	allAboutCycles();

	return 0;

};




void TSimpleMolecule::addAtom(TSingleAtom * sa){
	fAtom.push_back(sa);
};

void TSimpleMolecule::addAtom(int na, int charge, double rx, double ry) {
	TSingleAtom * sa;
	sa=new TSingleAtom();
	sa->na=na;
	sa->nv=hVal[na];
	sa->nc=charge;
	sa->rx=rx;
	sa->ry=ry;
	fAtom.push_back(sa);
};

void TSimpleMolecule::addBond(TSingleBond * sb){
	fBond.push_back(sb);
};

void TSimpleMolecule::addBond(int tb, int at1, int at2){
	TSingleBond * sb;
	sb=new TSingleBond();
	sb->tb=tb;
	sb->at[0]=at1;
	sb->at[1]=at2;
	fBond.push_back(sb);
};

void TSimpleMolecule::addFragmentAtPos(TSimpleMolecule & otherMol, int an, double posX, double posY) {
  double thisLength,otherLength,x,y;
  int i,naStore;
  TSingleAtom * sa;
  TSingleBond * sb;

  thisLength=this->averageAtomDistance();
  if (thisLength == 0) thisLength=1;
  otherLength=otherMol.averageAtomDistance();
  if (otherLength == 0) otherLength=1;

  naStore=this->nAtoms();

  x=0; y=0;
  if (an>=0) {
    x=otherMol.getAtom(an)->rx;
    y=otherMol.getAtom(an)->ry;
  };
  for (i=0; i<otherMol.nAtoms(); i++) {
    sa=otherMol.getAtom(i);
	this->addAtom(sa->clone());
	sa=this->getAtom(this->nAtoms()-1);
	sa->rx=posX+(sa->rx-x)*thisLength/otherLength;
	sa->ry=posY+(sa->ry-y)*thisLength/otherLength;
  };
  for (i=0; i<otherMol.nBonds(); i++) {
	sb=otherMol.getBond(i);
	this->addBond(sb->tb,sb->at[0]+naStore,sb->at[1]+naStore);
  };
};

void TSimpleMolecule::moleculeCopy(TSimpleMolecule & source) {
	TSingleAtom * sa;
	TSingleBond * sb;

	clear();
	for (int i=0; i<source.nAtoms(); i++) {
		sa=source.getAtom(i);
		fAtom.push_back(sa->clone());
	};
	for (int i=0; i<source.nBonds(); i++) {
		sb=source.getBond(i);
		fBond.push_back(sb->clone());
	};
	defineAtomConn();
	allAboutCycles();
	if (refofs == NULL) refofs=source.refofs;
};

void TSimpleMolecule::moleculeCopy(const TSimpleMolecule & source) {
	//TSingleAtom * sa;
	//TSingleBond * sb;

	clear();
	for (int i=0; i<source.nAtoms(); i++) {
		const TSingleAtom * sa = source.getAtom(i);
		fAtom.push_back(sa->clone());
	};
	for (int i=0; i<source.nBonds(); i++) {
		const TSingleBond * sb = source.getBond(i);
		fBond.push_back(sb->clone());
	};
	defineAtomConn();
	allAboutCycles();
	if (refofs == NULL) refofs=source.refofs;
};

void TSimpleMolecule::deleteBond(int index) {
	std::vector<TSingleBond*> tempBond(nBonds() + nAtoms());
	int i,n;

	n=0;
	//bond deletion
	for (i=0; i<nBonds(); i++) if (i != index) {
		tempBond[n]=getBond(i);
		n++;
	} else {
		delete(getBond(i));
		fBond[i]=NULL;
	};
	fBond.resize(n);
	for (i=0; i<n; i++) fBond[i]=tempBond[i];
	tempBond.clear();
};


void TSimpleMolecule::deleteAtom(int index) {
	std::vector<TSingleAtom*> tempAtom(nAtoms()-1);
	std::vector<TSingleBond*> tempBond(nBonds());
	int i,n;

	n=0;
	//atom deletion
	for (i=0; i<nAtoms(); i++) if (i != index) {
		tempAtom[n]=getAtom(i);
		n++;
	} else {
		delete(getAtom(i));
		fAtom[i]=NULL;
	};
	fAtom.resize(n);
	for (i=0; i<n; i++) fAtom[i]=tempAtom[i];
	tempAtom.clear();
	//connected bonds deletion
	n=0;
	for (i=0; i<tempBond.size(); i++) {
		if ((getBond(i)->at[0] == index) || (getBond(i)->at[1] == index)) {
			delete(getBond(i));
			fBond[i]=NULL;
		} else {
			if (getBond(i)->at[0] > index) getBond(i)->at[0]=getBond(i)->at[0]-1;
			if (getBond(i)->at[1] > index) getBond(i)->at[1]=getBond(i)->at[1]-1;
			tempBond[n]=getBond(i);
			n++;
		};
	};
	fBond.resize(n);
	for (i=0; i<n; i++) fBond[i]=tempBond[i];
	tempBond.clear();
};


double TSimpleMolecule::bondLength(int index) const {
	double result=1;
	double x,y;
	int k1,k2;
	k1=getBond(index)->at[0];
	k2=getBond(index)->at[1];
	x=getAtom(k1)->rx-getAtom(k2)->rx;
	y=getAtom(k1)->ry-getAtom(k2)->ry;
	result=sqrt(x*x+y*y);
	return result;
};

int TSimpleMolecule::listarSize() const {
	int result=10;  //Minimal vector size 10
	if (nAtoms() > result) result=nAtoms();
	if (nBonds() > result) result=nBonds();
	return result;
};

void TSimpleMolecule::getFragment(int & an, int & bn, TSimpleMolecule & fragment) {
  int i,j,k,n,at1,at2,atn;
  //std::vector<int> atomList(nAtoms());
  bool test=true;
  TSingleAtom * sa;

  fragment.clear();
  if((an < 0) && (bn < 0)) return;

  this->defineAtomConn();
  for (i=0; i<nAtoms(); i++) this->getAtom(i)->enumerator=-1;
  atn=an;
  if (atn < 0) atn=this->getBond(bn)->at[0];
  this->getAtom(atn)->enumerator=1;
  n=1;
  while (test) {
    test=false;
	for (i=0; i<nAtoms(); i++) if (this->getAtom(i)->enumerator == 1) {
	  sa=this->getAtom(i);
	  sa->enumerator=2;
	  for (j=0; j<sa->nb; j++) {  //All neighbours
		k=sa->ac[j];  
		if (this->getAtom(k)->enumerator == -1) {
		  this->getAtom(k)->enumerator=1;
		  test=true;
		  n++;
		};
	  };
	};
  };
  if (n == this->nAtoms()) {  //Molecule is then single fragment. No atom and bond enumeration is required
	fragment.moleculeCopy(* this);
	return;
  };
  for (i=0; i<nAtoms(); i++) if (this->getAtom(i)->enumerator != -1) {
    sa=this->getAtom(i);
	fragment.addAtom(sa->na,sa->nc,sa->rx,sa->ry);
	sa->enumerator=fragment.nAtoms()-1;
	if (i == an) an=fragment.nAtoms()-1;  //Enumeration atom number
  };
  for (i=0; i<nBonds(); i++) if (this->getAtom(this->getBond(i)->at[0])->enumerator != -1) {
    at1=this->getAtom(this->getBond(i)->at[0])->enumerator;
    at2=this->getAtom(this->getBond(i)->at[1])->enumerator;
	fragment.addBond(this->getBond(i)->tb,at1,at2);
	if (bn == i) bn=fragment.nBonds()-1;   //Enumeration bond number
  };
  fragment.defineAtomConn();
  fragment.allAboutCycles();
};

bool TSimpleMolecule::makeFragment(int& n, std::vector<int>& list, int aT, int aTEx) const {
	//Starting from atom's number AT in bond-connection matrix invariants array CONN
	//the all connected atom's (through one or more bonds) numbers are inserted into
	//the LIST array. Total number of such atoms is stored in N. If ATEX non-equal
	//zero, than atom with this number and all other, which is connected through
	//atom ATEX, don't include in LIST array. Boolean variable TEST has TRUE value
	//on output if between AT and ATEX exist cyclic bond, false otherwise}
	int i,j,k,l,an;
	bool test, test1;

	test=false;
	list.resize(listarSize());
	if ((nAtoms() ==0 ) || (aT < 0) || (aT >= nAtoms())) {
		n=0;
		list[0]=-1;
		list.resize(0);
		return test;
	}
	n=1;
	for (i = 0; i < list.size(); i++) list[i] = -1;
	list[0]=aT;
	for (i=0; i<getAtom(aT)->nb; i++) if (getAtom(aT)->ac[i] != aTEx) {
		list[n]=getAtom(aT)->ac[i];
		n++;
	}
	if (n == 1) {
		list.resize(1);
		return test;                      //exit in AT non-connected atom
	};
	k=0;
	do {                                 //iterational algorithm to define all}
		an = list[k];
		for (i=0; i<getAtom(an)->nb; i++) {
			test1=false;
			l=getAtom(an)->ac[i];
			if (l != aTEx) {
				for (j = 0; j < n; j++) if (list[j] == l) {
					test1 = true;
					break;
				};
			} else test=true;
			if ((! test1) && (l != aTEx) && (l >= 0) && (l < nAtoms())) {
				list[n]=l;
				n++;
			}
		}
		k=k+1;                             //atoms which are connected with AT}
	} while (k<n);
	list.resize(n);
	return test;
}

void TSimpleMolecule::readConnectionMatrix(const std::vector<int>iA1, const std::vector<int>iA2, int nAtoms, int nBonds) {
	TSingleAtom * sa;
	TSingleBond * sb;
	int i;

	clear();
	srand(30000);
	for (i=1; i<=nAtoms; i++) {
		sa=new TSingleAtom();
		sa->na=6;
		sa->nv=hVal[sa->na];
		sa->rx=(double)rand()/1000.0;
		sa->ry=(double)rand()/1000.0;
		fAtom.push_back(sa);
	};

	for (i=0; i<nBonds; i++) {
		sb=new TSingleBond();
		sb->at[0]=iA1[i];
		sb->at[1]=iA2[i];
		sb->tb=1;
		fBond.push_back(sb);
	};
	defineAtomConn();
	allAboutCycles();
};

void TSimpleMolecule::readConnectionMatrix(const std::vector<int>iA1, const std::vector<int>iA2, const std::vector<double>rx, const std::vector<double>ry, int nAtoms, int nBonds) {
	readConnectionMatrix(iA1,iA2,nAtoms,nBonds);
	for (int i=0; i<nAtoms; i++) {
		getAtom(i)->rx=rx[i];
		getAtom(i)->ry=ry[i];
	};
};


void TSimpleMolecule::normalizeCoordinates(double aveBL) {
	int i;
	double d,xMin,yMin;

	if (nAtoms() == 0) return;

	d=averageBondLength();
	if ((d > 0) && (aveBL>0)) for (i=0; i<nAtoms(); i++) {
		getAtom(i)->rx=getAtom(i)->rx*aveBL/d;
		getAtom(i)->ry=getAtom(i)->ry*aveBL/d;
	};
	//Normalizing X and Y coordinates
	xMin=getAtom(0)->rx;
	yMin=getAtom(0)->ry;
	for (i=0; i<nAtoms(); i++) {
		if (getAtom(i)->rx < xMin) xMin=getAtom(i)->rx;
		if (getAtom(i)->ry < yMin) yMin=getAtom(i)->ry;
	};
	for (i=0; i< nAtoms(); i++) {
		getAtom(i)->rx=getAtom(i)->rx-xMin+aveBL;
		getAtom(i)->ry=getAtom(i)->ry-yMin+aveBL;
	};
};


int TSimpleMolecule::nAtoms() const {
	return fAtom.size();
};

int TSimpleMolecule::nBonds() const {
	return fBond.size();
};


TSingleAtom * TSimpleMolecule::getAtom(int index) {
	if ((index<0) || (index >= fAtom.size())) {
      throw;
	} else return (TSingleAtom *)fAtom.at(index);
};

TSingleBond * TSimpleMolecule::getBond(int index) {
	return (TSingleBond *)fBond.at(index);
};

const TSingleAtom * TSimpleMolecule::getAtom(int index) const{
	return (const TSingleAtom *)fAtom.at(index);
};

const TSingleBond * TSimpleMolecule::getBond(int index) const{
	return (const TSingleBond *)fBond.at(index);
};


void TSimpleMolecule::defineAtomConn() {
	//for each atom returns list of adjusted atoms in atomConnection
	int i,n1, n2; 
	TSingleAtom * sa;

	for (i=0; i<nAtoms(); i++) {
		getAtom(i)->nb=0;
		getAtom(i)->currvalence=0;
	};
	for (i=0; i<nBonds(); i++) {
		n1=getBond(i)->at[0]; n2=getBond(i)->at[1];
		sa=getAtom(n1);
		sa->ac[sa->nb]=n2;
		sa->nb++;
		sa->currvalence=sa->currvalence+getBond(i)->getValence();
		sa=getAtom(n2);
		sa->ac[sa->nb]=n1;
		sa->nb++;
		sa->currvalence=sa->currvalence+getBond(i)->getValence();
	};
};

void TSimpleMolecule::defineBondConn(std::vector<adjustedlist> & bondConnection) const {
	//for each atom returns list of adjusted bonds in bondConnection
	int i,n1,n2; 

	for (i=0; i<nAtoms(); i++) bondConnection[i].nb=0;
	for (i=0; i<nBonds(); i++) {
		n1=getBond(i)->at[0]; n2=getBond(i)->at[1];
		bondConnection[n1].adjusted[bondConnection[n1].nb]=i;
		bondConnection[n1].nb++;
		bondConnection[n2].adjusted[bondConnection[n2].nb]=i;
		bondConnection[n2].nb++;
	};
};

void TSimpleMolecule::calculateRMGTypes(std::vector<adjustedlist> * bC) {
	//for each atom returns list of adjusted bonds in bondConnection
	int i,j,n1, n2, nn, nArom, nDouble, nTriple, nSingle, aCount, bCount; 
	TSingleAtom * atom;
	bool test,testCreate;
	std::vector<adjustedlist> * bondConnection=NULL;


	testCreate=bC == NULL;
	if (testCreate) {
		bondConnection = new std::vector<adjustedlist>;   //(neigbourlist *)malloc((1+CONNMAX)*NBONDSMAXSMBIG*4);
		bondConnection->resize(listarSize());

	  /*
	  for (i=0; i<=nAtoms(); i++) (*bondConnection)[i].nb=0;
	  for (i=0; i<nBonds(); i++) {
		n1=getBond(i)->at[0]; n2=getBond(i)->at[1];
		(*bondConnection)[n1].adjusted[(*bondConnection)[n1].nb]=i;
		(*bondConnection)[n1].nb++;
		(*bondConnection)[n2].adjusted[(*bondConnection)[n2].nb]=i;
		(*bondConnection)[n2].nb++;
	  };
	  */
	  defineBondConn(*bondConnection);
	} else bondConnection=bC;
	//RMG types determination
	aCount=nAtoms();
	for (i=0; i < aCount; i++) {
	  atom=getAtom(i);
  	  nn=0;
	  atom->rmgAtom[nn]=R_ATOM;
	  nn++;
	  switch (atom->na) {
	    case 6:
		  //CBF_ATOM and  queryin
		  nArom=0; nDouble=0; nTriple=0;
		  bCount=(*bondConnection)[i].nb;
          for (j=0; j<bCount; j++) {
			n2=getBond((*bondConnection)[i].adjusted[j])->db;
			if ((n2 == 2) || (n2 == 3)) nArom++; else {
  			  n2=getBond((*bondConnection)[i].adjusted[j])->tb;
			  if (n2 == 2) nDouble++; else if (n2 == 3) nTriple++;
			};
		  };
		  if (nArom == 2) {
			atom->rmgAtom[nn]=CB_ATOM;
			nn++;
		  } else if (nArom == 3) {
			atom->rmgAtom[nn]=CBF_ATOM;
			nn++;
		  } else if (nDouble == 2) {
			atom->rmgAtom[nn]=CDD_ATOM;
			nn++;
			//CO checking  ???????? 
		  } else if (nTriple == 1) {
			atom->rmgAtom[nn]=CT_ATOM;
			nn++;
		  } else if (nDouble == 0) {
			atom->rmgAtom[nn]=CS_ATOM;
			nn++;
		  } else if (nDouble == 1) {
			//CO checking
			test=false;
			for (j=0; j<(*bondConnection)[i].nb; j++) {
			  n1=(*bondConnection)[i].adjusted[j];
			  if (getBond(n1)->tb == 2) {          
				n2=getBond(n1)->at[0];
			    if (getAtom(n2)->na == 8) {
				  test=true;
				  break;
				};
				n2=getBond(n1)->at[1];
			    if (getAtom(n2)->na == 8) {
				  test=true;
				  break;
				};
			  };
			};
			if (test) {
  			  atom->rmgAtom[nn]=CO_ATOM;
			  nn++;
			} else {
			//?????  the below 
			  atom->rmgAtom[nn]=CD_ATOM;
			  nn++;
			};
		  };          
          break;
		case 8:
		  nDouble=0; nSingle=0;
          for (j=0; j<(*bondConnection)[i].nb; j++) {
  		    n2=getBond((*bondConnection)[i].adjusted[j])->tb;
			if (n2 == 2) nDouble++; else if (n2 == 1) nSingle++;
		  };
		  if (nDouble == 1) {
		    atom->rmgAtom[nn]=OD_ATOM;
			nn++;
		  } else {
		    atom->rmgAtom[nn]=OS_ATOM;
			nn++;
		  };
		  break;
		case 16:
		  nDouble=0; nSingle=0;
          for (j=0; j<(*bondConnection)[i].nb; j++) {
  		    n2=getBond((*bondConnection)[i].adjusted[j])->tb;
			if (n2 == 2) nDouble++; else nSingle++;
		  };
		  if (nDouble == 0) {
		    atom->rmgAtom[nn]=SS_ATOM;
			nn++;
		  };
		  break;
	  };
	};
	if (testCreate) delete bondConnection;
};


void TSimpleMolecule::newB(const std::vector<adjustedlist> & bk, int bnum, int anum, int & total, int * e, int * e1) const {
	/*for atom's number ANUM calculates the bond's list, which are connected
	with the ANUM, excluding bond BNUM. The bond's list is inserted into array
	L, TOTAL-number of components in the array. Array L1 contains atom's number
	which is connected with atom ANUM through corresponding bond in the array L.*/
	int n,i;

	total=0;
	for (i=0; i<bk[anum].nb; i++) if (bk[anum].adjusted[i] != bnum) {
		n=bk[anum].adjusted[i];  //bond number zero-based
		e[total]=n;
		if (getBond(n)->at[0] == anum) e1[total]=getBond(n)->at[1]; else e1[total]=getBond(n)->at[0]; 
		total++;
	};
};

void TSimpleMolecule::singleVawe(const std::vector<adjustedlist> & bk, std::vector<int> & alreadyDefined,
	std::vector<int> & prevBond, std::vector<int> & prevAtom,
	int & nPrev, std::vector<int> & curBond, std::vector<int> & curAtom) const {
		/*
		a vawe algorithm is used to calculate the neighbour sphere. PREVBOND and
		PREVATOM contain on input the bond's and atom's list correspondingly, the
		atom has been took into consideration at previous step. NPREV-number of
		components in PREVBOND and PREVATOM on input.  ALREADYDEFINED is the global
		array for procedure, which calls the SINGLEVAWE. It contains for each bond's
		number either zero, if the bond has not been took into consideration yet or
		bond's number from previous neighbour sphere. The subroutine makes next:
		1) calculates the new neighbour sphere for each bond from PREVBOND.
		2) if some bonds were not defined, they are labelled in ALREADYDEFINED,
		they numbers being stored into PREVBOND and PREVATOM. So, on output
		the arrays contain information for new neighbour sphere.
		3) On output NPrev is changed so new number of component in PREVBOND
		is lighted}
		*/
		int e[CONNMAX];
		int e1[CONNMAX];
		int i,j,ncur;
		int n;

		ncur=0;
		n=0;
		for (i=0; i<nPrev; i++) {
			newB(bk,prevBond[i],prevAtom[i],n,e,e1);
			if (n != 0) for (j=0; j<n; j++) if (alreadyDefined[e[j]]<0) {
				curBond[ncur]=e[j];
				curAtom[ncur]=e1[j];
				ncur++;
				alreadyDefined[e[j]]=prevBond[i];
			};
		};
		nPrev=ncur;
		for (i=0; i<ncur; i++) prevBond[i]=curBond[i]; 
		for (i=0; i<ncur; i++) prevAtom[i]=curAtom[i]; 
};

void TSimpleMolecule::atomDistanceList(int atomNo, int prohibitedAtomNo, std::vector<int> & atomDistance) const {
	std::vector<int> currentList, accumulatedList;
	int sphereN, i, j, n, nAddedToCalculate;

	atomDistance.resize(nAtoms());
	for (i = 0; i < atomDistance.size(); i++) atomDistance[i] = -1;
	atomDistance[atomNo] = 0;
	if (prohibitedAtomNo>=0) atomDistance[prohibitedAtomNo] = 1000000; 
	for (sphereN = 0; sphereN < nAtoms(); sphereN++) {
		nAddedToCalculate = 0;
		for (i = 0; i < nAtoms(); i++) if (atomDistance[i] == sphereN) {
			for (j = 0; j < getAtom(i)->nb; j++) {
				n = getAtom(i)->ac[j];
				if (atomDistance[n] == -1) {
					atomDistance[n] = sphereN + 1;
					nAddedToCalculate++;
				}
			}
		}

		if (nAddedToCalculate == 0) break;
	}
};


void TSimpleMolecule::vaweBond(int bondN, const std::vector<adjustedlist> & bk,
	int & ringSize, std::vector<int> & bondList, int maxSize) const {
		/*The procedure for bond's number BONDN determines, whether or not the bond
		belongs to ring (RINGSIZE=0 if not ring bond and RINGSIZE=number of atoms
		in the smallest ring, to which the bond is assigned). If the bond belongs
		to ring, array BONDLIST contains on output the bond's numbers, which are
		also assigned to the same ring */

		int i,j,k,aA,i1,i2;
		std::vector<int>pA(listarSize());
		std::vector<int>pB(listarSize());
		std::vector<int>oD(listarSize());
		std::vector<int>curBond(listarSize());
		std::vector<int>curAtom(listarSize());
		bool test;
		int nP;


		for (j=0; j<nBonds(); j++) oD[j]=-1;
		oD[bondN]=65500;
		ringSize=1;
		nP=1;  //OK for zer-based in C
		pB[0]=bondN; //pB.setValue(1,bondN);
		pA[0]=getBond(bondN)->at[0];//pA.setValue(1,fBond.getAT(bondN,1));
		aA=getBond(bondN)->at[1];//fBond.getAT(bondN,2);

		k=0;
		test=false;
		while (! ((nP==0) || test)) {                        //recursion
			ringSize++;
			singleVawe(bk,oD,pB,pA,nP,curBond,curAtom);      //new neighbour sphere generation}
			test=false;
			if (nP>0) for (j=0; j<nP; j++) if (pA[j] == aA) {
				test=true;
				k=pB[j];
			};
			if ((maxSize > 0) && (ringSize >= maxSize)) {
				break;
			};
		}; //recursion finishing
		if (test) {                    //bond belongs to ring
			bondList.resize(ringSize);
			bondList[ringSize-1]=bondN;         //bond BONDN must be last in the list
			for (j=1; j<ringSize; j++) {        //BONDLIST formation
				bondList[j-1]=k;
				k=oD[k];
			};
			//Sorting remaining bonds in ascending order - need for cycle comparison
			k=2;                   //all bonds are sorted 
			//remains last bond in cycle description=bondN (input parameter)
			for (i=0; i<(ringSize-k); i++) for (j=i+1; j<=(ringSize-k); j++) {
				i1=bondList[i];
				i2=bondList[j];
				if (i1 > i2) {
					bondList[i]=i2;	 
					bondList[j]=i1;	 
				}; 
			};
		} else ringSize=0;
};


//This method must be called ONLY from AllAboutCycles
bool TSimpleMolecule::aromatic(int cycleSize, const std::vector<int> bondList, 
	std::vector<int>& arom) {

		/*for a cycle with the dimension of CYCLESIZE determines, whether or not the
		cycle is aromatic. BONDLIST contains the bond's numbers, which are belongs
		to the cycle. AROM-global array for the structure under study, for each bond
		in the structure contains label, whether the bond is aromatic or not. On out-
		put function AROMATIC contains TRUE if the ring is aromatic, FALSE otherwise
		*/
		bool test;
		int i,j,k,n,n1,m;
		int doubleDetected[11];
		int atomDetected[11];
		bool result=false;

		if ((cycleSize<5) || (cycleSize>6)) return result;
		//procedure define aromatic only for 5 or 6 membered cycles}
		//Simple Test}
		n=0;
		for (i=0; i<cycleSize; i++) {
			k=bondList[i];

			if ((getBond(k)->tb == 1) && (arom[k]==0)) doubleDetected[i]=0; else
				if ((getBond(k)->tb == 2) || (getBond(k)->tb == 4) || (arom[k]>0)) {
					n=n+1;                 //number of double or aromatic bonds calculation}
					doubleDetected[i]=1;
				} else return result;
		};
		if (cycleSize==6) {     //6-membered ring}
			if (n>=3) {
				result=true;
				if (n<=4) for (i=0; i<cycleSize; i++) if (doubleDetected[i]==0) {
					//Checking if not connected
					for (j=0; j<cycleSize; j++) if ((j != i) && (doubleDetected[j]==0)) {
						test=(getBond(bondList[i])->at[0] == getBond(bondList[j])->at[0]) || (getBond(bondList[i])->at[0] == getBond(bondList[j])->at[1])
							|| (getBond(bondList[i])->at[1] == getBond(bondList[j])->at[0]) || (getBond(bondList[i])->at[1] == getBond(bondList[j])->at[1]);
						if (test) result=false;
						if (! result) break;
					};
					if (! result) break;
				};
			};
			return result;
		};

		if (n<2) return result;        //5-membered ring
		m=0;
		for (i=0; i<cycleSize; i++) if (doubleDetected[i]==0) {
			atomDetected[m]=getBond(bondList[i])->at[0];
			m++;
			atomDetected[m]=getBond(bondList[i])->at[1];
			m++;
		};
		n=-1;
		if (m>0) for (i=0; i<(m-1); i++) for (j=i+1; j<m; j++) if (atomDetected[i]==atomDetected[j]) n=i;
		//Addition from January 2001
		if (n>=0) {
			n1=-1;
			for (i=0; i<(m-1); i++) if (i != n) for (j=i+1; j<m; j++) if (atomDetected[i]==atomDetected[j]) n1=i;
			if (n1>=0) return result;
		};
		if (n < 0) {   //No common atoms...
			result=m==4;
			return result;
		};
		test=false;
		i=0;
		j=getAtom(atomDetected[n])->na;
		for (i=0; i<NAROMMAX; i++) if (j == possibleAromatic[i]) {
			test=true;
			break;
		};
		if (! test) test=(j==6) && (getAtom(atomDetected[n])->nc < 0);
		//cyclopentadienyl checking}
		result=test;
		return result;
};

void TSimpleMolecule::allAboutCycles(int maxSize) {
	/*
	for each bond in the structure cyclic conditions are calculated and stored
	into the BOND[I].DB variables. Possible values:
	0 - chain;
	1 - chain, appended to cycle;
	2 - aromatic, 5-membered ring;
	3 - aromatic, 6-membered ring;
	4...more = non-aromatic ring size+1}
	*/
	std::vector<int>bondTested(listarSize());
	std::vector<int>bondList(listarSize());
	std::vector<int>ar(listarSize());
	std::vector<int>rSize(listarSize());
	std::vector<int>ar1(listarSize());
	std::vector<int>cycleDescription(3* listarSize());
	std::vector<int>cycleAddress(listarSize());

	int i,n,j,k,at1,at2,an;
	int cycleSize;
	std::vector<adjustedlist> bk;
	bool test,isPlanar;



	if (nBonds() == 0) return;
	bk.resize(listarSize());//  (neigbourlist *)malloc((1 + CONNMAX)*NBONDSMAXSMBIG * 4);

	defineBondConn(bk);
	//initial values for MainList-number of bonds, connected to each atom and their numbers in array BOND
	for (i=0; i<nBonds(); i++) {
		getBond(i)->db=0;
		ar[i]=0;;
		bondTested[i]=0;
	};
	n=0;
	cycleAddress[0]=0;
	for (i = 0; i < nBonds(); i++) if (bondTested[i] == 0) {
		at1 = getBond(i)->at[0];
		at2 = getBond(i)->at[1];
		test = ((getAtom(at1)->nb == 1) || (getAtom(at1)->nb == 1));
		if (!test) {
	  	  vaweBond(i, bk, cycleSize, bondList, maxSize); //Is I-th bond cyclic?
		  if (cycleSize > 0) {       //Yes, cyclic
			  canonizeCycle(cycleSize, bondList);

			  if (false) { //(savePlanar) {
				  isPlanar = isPlanarRing(cycleSize, bondList, bk);
				  an = getBond(i)->at[0];
				  getAtom(an)->setRingProperty(cycleSize, isPlanar);
				  an = getBond(i)->at[1];
				  getAtom(an)->setRingProperty(cycleSize, isPlanar);

			  }

			  n++;
			  cycleAddress[n] = cycleAddress[n - 1] + cycleSize;  //store cycle definitions
			  rSize[n - 1] = cycleSize;
			  for (j = 0; j < cycleSize; j++) cycleDescription[cycleAddress[n - 1] + j] = bondList[j];
			  ar1[n - 1] = 0;
			  if (aromatic(cycleSize, bondList, ar)) {       //store, if aromatic cycle}
				  ar1[n - 1] = 1;
				  for (j = 0; j < cycleSize; j++) ar[bondList[j]] = 1;
			  };
			  for (j = 0; j < cycleSize; j++) bondTested[bondList[j]] = 1;
			  //marking all bonds in the cycle}
		  };
	    };
		bondTested[i]=1;                                //marking I-th bond
	};
	bondList.resize(listarSize());
	if (n>0) {
		j=1;
		while (j != 0) {
			j=0;
			for (i=0; i<n; i++) if (ar1[i]==0) {   //search for condenced aromatics
				cycleSize=rSize[i];
				if (cycleSize<=6) for (k=0; k<cycleSize; k++) bondList[k]=cycleDescription[cycleAddress[i]+k];
				if (aromatic(cycleSize,bondList,ar)) {
					j++;
					ar1[i]=1;
					for (k=0; k<cycleSize; k++) ar[bondList[k]]=1;
				};
			};
		}; //until J=0;          {until new aromatics are not detected under iteration}
		for (i=0; i<n; i++) if (ar1[i]==0) {       //Label of non-aromatic cycles in BOND
			cycleSize=rSize[i];
			for (k=0; k<cycleSize; k++) bondList[k]=cycleDescription[cycleAddress[i]+k];  //!!!
			for (j=0; j<cycleSize; j++) {
				k=cycleSize+1;
				if ((getBond(bondList[j])->db==0) || (k<getBond(bondList[j])->db)) getBond(bondList[j])->db=k;
			};
		};
		for (i=0; i<n; i++) if (ar1[i]==1)  {     //Label of aromatic cycle in BOND
			cycleSize=rSize[i];
			if (cycleSize > bondList.size()) bondList.resize(cycleSize);
			for (k=0; k<cycleSize; k++) bondList[k]=cycleDescription[cycleAddress[i]+k];  //!!!
			for (j=0; j<cycleSize; j++) {
				if (cycleSize==5) getBond(bondList[j])->db=2;
				if ((cycleSize==6) && (getBond(bondList[j])->db != 2)) getBond(bondList[j])->db=3;
			};
		};
	};
	//free(bk);
};


void TSimpleMolecule::redrawMolecule() {
	int i;
	int atomClean;
	int bondClean;
	std::vector<int>listAtomClean(nAtoms());
	std::vector<int>listBondClean(nBonds());

	if (nAtoms()==0) return;
	for (i=0; i<nAtoms(); i++) listAtomClean[i]=i;
	for (i=0; i<nBonds(); i++) listBondClean[i]=i;
	atomClean=nAtoms();
	bondClean=nBonds();
	redraw(listAtomClean,listBondClean,atomClean,bondClean,1,1,0,false);
};


void TSimpleMolecule::twoAtomUnitVector(int na1, int na2, double & xv, double & yv, const std::vector<int>atomDefine) {
	/*Using atom's numbers NA1 and NA2 in the array ATOMS calculates the direction
	for bond between NA1 and NA2 (XV,YV on output). At this direction new
	substitutor may be added with optimal place.*/
	double si,r1,r2,r3,s1,s2,s3,s4,x1,y1,y2;
	int i,k;

	if ((getAtom(na1)->rx==getAtom(na2)->rx) && (getAtom(na1)->ry==getAtom(na2)->ry)) {
		xv=1;
		yv=1;
		return;
	};
	s1=getAtom(na1)->rx;
	s2=getAtom(na1)->ry;
	s3=getAtom(na2)->rx;
	s4=getAtom(na2)->ry;
	r1=s1-s3;
	r2=s2-s4;
	r3=sqrt(r1*r1+r2*r2);
	r1=r1/r3;  r2=r2/r3;
	si=0;
	for (i=0; i<getAtom(na1)->nb; i++) {
		k=getAtom(na1)->ac[i];
		if (k != na2) if (atomDefine[k] > 0) {
			x1=getAtom(k)->rx-s1;
			y1=getAtom(k)->ry-s2;
			y2=x1*r2-y1*r1;
			if (y2 != 0) si=si+y2/abs(y2);
		};
	};
	for (i=0; i<getAtom(na2)->nb; i++) {
		k=getAtom(na2)->ac[i]; 
		if (k != na1) if (atomDefine[k] > 0) {
			x1=getAtom(k)->rx-s3;
			y1=getAtom(k)->ry-s4;
			y2=x1*r2-y1*r1;
			if (y2 != 0) si=si+y2/abs(y2);
		};
	};
	if (si != 0) si=si/abs(si); else si=1;
	xv=-r2*si; yv=r1*si;
};

void TSimpleMolecule::defC(int& currNumDef, int baseCycle, int atomClean, std::vector<int>& cycleDefine,
	std::vector<int>& atomDefine, std::vector<int>& atomCycle, std::vector<int>& cycleAddress, 
	std::vector<int>& cycleSize, std::vector<int>& dsATN, std::vector<int>& dsTP, 
	std::vector<int>& dsSC, std::vector<int>& dsNA1, std::vector<int>& dsNA2) {
		//The procedure create priority list formation for cleaning of molecule. Recur-
		// sive calls to the procedure are required together with DefA. After recursion
		// have been finished, the DEFINESEQUENCE array will be created. The procedure
		// analyze cyclic fragments
		int aA,cN,maxAtDef,i,j,atDef;
		bool test;

		test=true;
		while (test) {
			cN=-1;
			if ((baseCycle==0) || (currNumDef==atomClean)) return;
			maxAtDef=0;
			for (i=0; i<baseCycle; i++) if (cycleDefine[i]==0) {
				//        {for all cycles, which were not included into priority list}
				atDef=0;
				//        {search for maximal priority level}
				for (j=0; j<cycleSize[i]; j++) if (atomDefine[atomCycle[cycleAddress[i]+j]]>0) atDef=atDef+1;
				if (atDef > maxAtDef) {
					//          {  search for cycles, for which maximal number of atoms were inserted in the priority list}
					maxAtDef=atDef;
					cN=i;
				} else if ((maxAtDef > 0) && (atDef==maxAtDef)) {
					//         {search for minimal cycle size}
					if (cycleSize[i]<cycleSize[cN]) cN=i;
				};
			};
			if (cN>=0) if (maxAtDef==cycleSize[cN]) {
				cycleDefine[cN]=1;
				cN=-1;
			};

			//          {CN-number of selected cycle with max.number of atoms, already defined, or with minimal size}
			if (cN>=0) {
				cycleDefine[cN]=1;  //{Label, that CN-th cycle is inserted into priority list}

				test=(atomDefine[atomCycle[cycleAddress[cN]]]>0) &&
					(atomDefine[atomCycle[cycleAddress[cN]+cycleSize[cN]-1]]==0);
				//        {if first atom from cyclic fragment was inserted into list previously and last-no}
				while (! test) {
					//          {enumeration of atom's sequence in cyclic fragment}
					aA=atomCycle[cycleAddress[cN]];
					for (j=1; j<cycleSize[cN]; j++) atomCycle[cycleAddress[cN]+j-1]=atomCycle[cycleAddress[cN]+j];
					atomCycle[cycleAddress[cN]+cycleSize[cN]-1]=aA;
					test=(atomDefine[atomCycle[cycleAddress[cN]]]>0) &&
						(atomDefine[atomCycle[cycleAddress[cN]+cycleSize[cN]-1]]==0);
				};
				//        {now first atom in the cycle definition list must be inserted into the priority list early, and the last atom-no}
				for (j=0;j<cycleSize[cN]-maxAtDef; j++) {
					//          for each undefined atom from selected cycle}
					//atom is added to priority list}
					dsATN[currNumDef]=atomCycle[cycleAddress[cN]+maxAtDef+j];

					atomDefine[dsATN[currNumDef]]=1;
					if (maxAtDef<3) dsTP[currNumDef]=maxAtDef+1; else dsTP[currNumDef]=4;
					//type of clean procedure definition}
					dsNA1[currNumDef]=atomCycle[cycleAddress[cN]+maxAtDef-1];
					dsNA2[currNumDef]=atomCycle[cycleAddress[cN]];
					dsSC[currNumDef]=cycleSize[cN]-maxAtDef;
					currNumDef++;

					/*       Initial Java codes
					//atom is added to priority list}
					dsATN.setValue(currNumDef.value,atomCycle.getValue(cycleAddress.getValue(cN)+maxAtDef+j));
					atomDefine.setValue(dsATN.getValue(currNumDef.value),1);
					if (maxAtDef<3) dsTP.setValue(currNumDef.value,maxAtDef+1); else dsTP.setValue(currNumDef.value,4);
					//type of clean procedure definition}
					dsNA1.setValue(currNumDef.value,atomCycle.getValue(cycleAddress.getValue(cN)+maxAtDef));
					dsNA2.setValue(currNumDef.value,atomCycle.getValue(cycleAddress.getValue(cN)+1));
					dsSC.setValue(currNumDef.value,cycleSize.getValue(cN)-maxAtDef);
					currNumDef.value=currNumDef.value+1;
					*/
				};
			};
			test=(cN<0) || (currNumDef==atomClean);
			test=! test;
		};
};

void TSimpleMolecule::defA(int& currNumDef, int atomClean, int sPN, int baseCycle, std::vector<int>& atomDefine, const std::vector<int> listAtomClean,
	std::vector<int>& cycleDefine, std::vector<int>& cycleSize, std::vector<int>& cycleAddress, std::vector<int>& atomCycle,
	std::vector<int>& dsATN, std::vector<int>& dsTP, std::vector<int>& dsNA1, std::vector<int>& dsNA2) {

		//The procedure create priority list formation for cleaning of molecule. Recur-
		// sive calls to the procedure are required together with DefC. After recursion
		//have been finished, the DEFINESEQUENCE array will be created. The procedure
		//analyze chain fragments
		int i,j,k,rC;



		if (currNumDef==atomClean) return;

		for (i=0; i<atomClean; i++) if (atomDefine[listAtomClean[i]]==0)
			if (getAtom(listAtomClean[i])->nb > 0) for (j=0; j<getAtom(listAtomClean[i])->nb; j++) {
				k=getAtom(listAtomClean[i])->ac[j]; 
				if (atomDefine[k]>0) {
					//if atom has a neighbour, which has already been inserted into priority
					// list, then it is inserted into the priority list}
					dsATN[currNumDef]=listAtomClean[i];
					atomDefine[dsATN[currNumDef]]=1;
					dsTP[currNumDef]=1;                  //type of clean}
					dsNA1[currNumDef]=k;
					dsNA2[currNumDef]=-1;
					currNumDef++;
					//only one atom is added into priority list, then it is necessary to
					// analyze cycles by DEFC procedure!}
					return;
				};
			};
		if ((sPN<3) || (sPN==4)) {   //{next may not be implemented to CleanGroup command}
			//{initializing priority list for clean of a new fragment}
			rC=0;
			j=100000;
			//minimal cycle size searching}
			if (baseCycle>0) for (i=0; i<baseCycle; i++) if (cycleDefine[i]==0) if (cycleSize[i]<j) {
				rC=i;
				j=cycleSize[i];
			};
			if (rC>0) i=atomCycle[cycleAddress[rC]]; else {
				//        {not found-first undefined atom is selected}
				i=0;
				while (atomDefine[listAtomClean[i]] != 0) i++;
				i=listAtomClean[i];
			};
		} else i=listAtomClean[atomClean-1];
		//for group-first atom, from which group started, is selected}
		dsATN[currNumDef]=i;
		atomDefine[dsATN[currNumDef]]=1;
		dsTP[currNumDef]=0;                  //type of clean}
		dsNA1[currNumDef]=-1;
		dsNA2[currNumDef]=-1;
		currNumDef++; //Addition the atom selected to the list}
};

void TSimpleMolecule::canonizeCycle(int ringSize, std::vector<int> & bondList) {
	//Order of nonds in cycle description in so way, that bond with minimal number 
	//is the first in this description. Then is coming bond, which is appended to 
	//atom with maximal number, then next and so on
	std::vector<int> bondUsed(ringSize);
	std::vector<int> newBondList(ringSize);
	int i,j,n,m,currentAtom;

	n=bondList[0];
	m=0;
	for (i=0; i<ringSize; i++) {
		bondUsed[i]=0;
		if (bondList[i]<n) {  //minimal bond number
			n=bondList[i];
			m=i;
		};
	};
	currentAtom=getBond(n)->at[0];
	if (getBond(n)->at[1] > currentAtom) currentAtom=getBond(n)->at[1];
	newBondList[0]=n;
	bondUsed[m]=1;
	n=1;
	for (i=1; i<ringSize; i++) {  //cycle size-1 because of 1-st bond already selected
		for (j=0; j<ringSize; j++) if (bondUsed[j] == 0) {
			m=bondList[j];
			if ((getBond(m)->at[0] == currentAtom) || (getBond(m)->at[1] == currentAtom)) {
				bondUsed[j]=1;
				newBondList[n]=m;
				n++;
				if (getBond(m)->at[0] == currentAtom) currentAtom=getBond(m)->at[1]; else currentAtom=getBond(m)->at[0];
				break;
			};
		};
	};
	//copy new set
	for (i=0; i<ringSize; i++)bondList[i]=newBondList[i];
};

/*
bool TSimpleMolecule::threeBondResolve(int an, int bondExcluded, double& xv, double& yv, neigbourlist & bk) {
//Addition from 16 April 2006
int bondNoList[3];
bool result=false;
double centerX[3];
double centerY[3];
int nBondNo;
std::vector<int> bondList(nBonds());
std::vector<int> blStore(nBonds());
bool testBad,testOK;
int m,j,n1,n2,n,k;
int rs;
double dist,x,y,minDist;
int i;
bool test,testStore;
int at;
double r1,r2,s1,s3;

xv=0; yv=1;
nBondNo=0;
testBad=false;
testOK=true;
testStore=false;
for (i=0; i<bk[an].nb; i++) {
n=bk[an].adjusted[i];
if (n != bondExcluded) {
vaweBond(n,bk,rs,bondList);
if (rs > 0) {
//I have to analize bondList to determine second cycle to exclude adamanthane
if (! testStore) {
//Save
blStore.clear();
for (j=0; j<rs; j++) blStore.push_back(bondList[j]);
testStore=true;
} else {
m=0;
for (j=0; j<blStore.size(); j++) {
test=false;
for (k=0; k<rs; k++) if (bondList[k] == blStore[j]) {
test=true;
break;
};
if (test) m++;
};
if ((m > 1) && (m < blStore.size())) {
testOK=false;
};
};
//center determination
centerX[nBondNo]=0;
centerY[nBondNo]=0;
for (j=0; j<rs; j++) {
m=bondList[j];
n1=getBond(m)->at[0];
n2=getBond(m)->at[1];
centerX[nBondNo]=centerX[nBondNo]+getAtom(n1)->rx+getAtom(n2)->rx;
centerY[nBondNo]=centerY[nBondNo]+getAtom(n1)->ry+getAtom(n2)->ry;
};
centerX[nBondNo]=centerX[nBondNo]/(2*rs);
centerY[nBondNo]=centerY[nBondNo]/(2*rs);
bondNoList[nBondNo]=n;
nBondNo++;
};
};// else testBad:=true;
if ((nBondNo == 3) || (testBad)) break;
};
if (nBondNo < 2) testBad=true;
if ((! testBad) && testOK) {
dist=0;
if (nBondNo == 2) {
//single bond, attached to ring. Re-definition center....
x=0; y=0;
for (i=0; i<nAtoms(); i++) {
x=x+getAtom(i)->rx;
y=y+getAtom(i)->ry;
};
x=x/(double)nAtoms();
y=y/(double)nAtoms();
for (i=0; i<nBondNo; i++) {
centerX[i]=x;
centerY[i]=y;
};
};
for (i=0; i<nBondNo; i++) {
at=getBond(bondNoList[i])->at[0];
if (at== an) getBond(bondNoList[i])->at[1];
r1=getAtom(at)->rx-getAtom(an)->rx;
r2=getAtom(at)->ry-getAtom(an)->ry;
s1=sqrt(r1*r1+r2*r2);
if (s1 == 0) testBad=true; else {
r1=r1/s1; r2=r2/s1;
x=-r1;
y=-r2;
minDist=1000000000;
for (j=0; j<nBondNo; j++) {
r1=(x*s1+getAtom(an)->rx-centerX[j]);
r2=(y*s1+getAtom(an)->ry-centerY[j]);
s3=sqrt(r1*r1+r2*r2);
if (s3 < minDist) minDist=s3;
};
if (minDist > dist) {
dist=minDist;
xv=x;
yv=y;
};
};
};
if (! testBad) result=true;
};
//end addition from 16 April 2006
return result;
};
*/

void TSimpleMolecule::redraw(const std::vector<int>listAtomClean, const std::vector<int>listBondClean, 
	int & atomClean, int & bondClean, int spn, int sCHA1, int sCHB1, bool iOPT7) {

	std::vector<int> dsATN(nBonds()+nAtoms());  //(NBONDSMAXSM);
	std::vector<int> dsTP(nBonds() + nAtoms());//(NBONDSMAXSM);
	std::vector<int> dsSC(nBonds() + nAtoms());// (NBONDSMAXSM);
	std::vector<int> dsNA1(nBonds() + nAtoms());// (NBONDSMAXSM);
	std::vector<int> dsNA2(nBonds() + nAtoms());// (NBONDSMAXSM);
		/*
		the arrays are the priority list for clean. Minimal element of record has
		the maximal priority. Each element of record contains:
		ATN-atom's number in the ATOM array, coordinate of which must be calculated
		at current step.
		TP -type of calculation of coordinates for the atom ATN:
		=0 - simple put at any position
		=1 - calculation from coordinates of previously defined atom NA1. NA2
		value hasn't effect in the case.
		=2 - for cyclic fragment, coordinates of all cyclic atoms must be cal-
		culated from the ones for atom NA1-base atom in cycle. This type
		means fragment: chain bond, attached to the cycle. NA2 value hasn't
		effect in the case.
		=3 - for cyclic fragment, coordinates of all cyclic atoms must be cal-
		culated from the ones for pair of atoms NA1 and NA2. The pair is
		connected with bond. This type means fragment in condenced-cyclic
		systems.
		=4 - for cyclic fragment, coordinates of all cyclic atoms must be cal-
		culated from the ones for pair of atoms NA1 and NA2. The pair is
		connected through more, then one bond. This type means fragment
		in polycyclic systems.
		AN1,AN2-auxiliary atom's numbers (see TP description).
		SC-cycle size minus number of already defined atoms. The variable is taken
		into consideration for TP=2-4.
		*/
		bool test;
		int cs;
		int add,nb,k,k1,k2,l,i,j,baseCycle,lx,ly;
		int currNumDef;
		int atomSecond;
		std::vector<int> atomCycle(3*(nBonds()+nAtoms()));
		std::vector<int> cycleDescription(3*(nBonds()+nAtoms()));
		std::vector<int> cycleAddress(nBonds() + nAtoms());
		std::vector<int> atomDefine(listarSize());
		std::vector<int> cycleSize(listarSize());
		std::vector<int> cycleDefine(listarSize());
		std::vector<int> tempArray(listarSize());
		double r,cf,fi,ux,uy,ux1,uy1,ux2,uy2,uvX,uvY,c,s;
		double xCenterOld,yCenterOld,xCenterNew,yCenterNew,bondLengthOld,bondLengthNew;
		std::vector<adjustedlist> bk;
		int n,n1,n2;
		double r1;
		bool isCycle,isChainFour;
		int mm1,mm2,mAny,bnEx;

		//Lx,Ly,Button:integer;

		if ((atomClean<1) || (bondClean==0)) return;
		defineAtomConn();
		allAboutCycles();

		test=true;
		bk.resize(listarSize());
		defineBondConn(bk);
		//{Start clean for LISTATOMCLEAN and LISTBONDCLEAN atoms and bonds}
		cycleAddress[0]=0;
		baseCycle=0;
		cs=0;
		ux=0;
		uy=0;
		uvX=0;
		uvY=0;
		atomSecond=0;

		for (i=0; i<bondClean; i++) {
			vaweBond(listBondClean[i],bk,cs,atomDefine); //If I-th bond belongs to cycle}

			if (cs>0) {
				canonizeCycle(cs,atomDefine);
				test=false;
				for (j=0; j<baseCycle; j++) if(cycleSize[j] == cs) {
					test=true;
					for (k=0; k<cycleSize[j]; k++) if (atomDefine[k] != cycleDescription[cycleAddress[j]+k]) {
						test=false;
						break;
					};
					if (test) break;
				};
				if (! test) {     //{not found-the cycle is added to cycle's list}
					baseCycle++;
					cycleAddress[baseCycle]=cycleAddress[baseCycle-1]+cs;
					cycleSize[baseCycle-1]=cs;
					for (j=0; j<cs; j++) cycleDescription[cycleAddress[baseCycle-1]+j]=atomDefine[j];
				};
			};

		};

		for (i=0; i<baseCycle; i++) { //making atom list in those order, as they will be defined at cycles}
			cycleDefine[i]=0;
			k=cycleDescription[cycleAddress[i]+0];
			j=cycleDescription[cycleAddress[i]+1];

			if ((getBond(k)->at[0]==getBond(j)->at[0]) || (getBond(k)->at[0]==getBond(j)->at[1])) atomCycle[cycleAddress[i]+0]=getBond(k)->at[0];
			else atomCycle[cycleAddress[i]+0]=getBond(k)->at[1];
			n=atomCycle[cycleAddress[i]+0]; //previoulsly putted atom
			for (j=1; j<cycleSize[i]; j++) {
				k=cycleDescription[cycleAddress[i]+j];
				if (getBond(k)->at[0]==n) atomCycle[cycleAddress[i]+j]=getBond(k)->at[1]; else atomCycle[cycleAddress[i]+j]=getBond(k)->at[0];
				n=atomCycle[cycleAddress[i]+j];
			};
		};

		//calculation of the coordinates of fragment's center and scaling factor}

		if ((spn<3) || (spn==4)) {
			xCenterOld=0;
			yCenterOld=0;
			bondLengthOld=0;
			if (spn == 4) {
				for (i=0; i<nAtoms(); i++) tempArray[i]=0;
				for (i=0; i<atomClean; i++) {
					k=listAtomClean[i];
					tempArray[k]=1;
				};
				n=0;
				for (i=0; i<nAtoms(); i++) if (tempArray[i]==0) {
					xCenterOld=xCenterOld+getAtom(i)->rx;
					yCenterOld=yCenterOld+getAtom(i)->ry;
					n++;
				};
				xCenterOld=xCenterOld/(double)n; yCenterOld=yCenterOld/(double)n;
				for (i=0; i<nBonds(); i++) tempArray[i]=0;
				for (i=0; i<bondClean; i++) {
					k=listBondClean[i];
					tempArray[k]=1;
				};
				n=0;
				bondLengthOld=1E10;
				r1=0;
				for (i=0; i<nBonds(); i++) if (tempArray[i] == 0) {
					r=this->bondLength(i);    
					if (r < bondLengthOld) bondLengthOld=n;
					r1=r1+r;
					n++;
				};
				r1=r1/n;
				if (5*bondLengthOld < r1) bondLengthOld=r1;
			} else {
				for (i=0; i<atomClean; i++) {
					k=listAtomClean[i];
					xCenterOld=xCenterOld+getAtom(k)->rx;
					yCenterOld=yCenterOld+getAtom(k)->ry;
				};
				xCenterOld=xCenterOld/(double)atomClean; yCenterOld=yCenterOld/(double)atomClean;
				for (i=0; i<bondClean; i++) {
					k=listBondClean[i];
					bondLengthOld=bondLengthOld+this->bondLength(k); 
				};
				bondLengthOld=bondLengthOld/(double)bondClean;
			};
			if (bondLengthOld<0.01) bondLengthOld=1;
		} else {

			xCenterOld=getAtom(sCHA1)->rx; //{CHA1 - To atom connected, CHB1-splitting bond}
			yCenterOld=getAtom(sCHA1)->ry;
			lx=getBond(sCHB1)->at[0];
			ly=getBond(sCHB1)->at[1];
			if (ly == sCHA1) {
				lx=getBond(sCHB1)->at[1];
				ly=getBond(sCHB1)->at[0];
			};
			atomSecond=ly;
			//for group it is necessary to calculate unit vector direction
			uvX=getAtom(ly)->rx-getAtom(lx)->rx;
			uvY=getAtom(ly)->ry-getAtom(lx)->ry;
			bondLengthOld=sqrt(uvX*uvX+uvY*uvY);
			if (bondLengthOld<0.01) bondLengthOld=1;
			uvX=uvX/bondLengthOld;
			uvY=uvY/bondLengthOld;

		};

		if (atomDefine.size() < nAtoms()) atomDefine.resize(nAtoms());
		for (i=0; i<nAtoms(); i++) atomDefine[i]=0; //flags-zero value OK
		if ((spn<3) || (spn==4)) {     //checking already cleaned atoms...
			for (i=0; i<nAtoms(); i++) atomDefine[i]=1;
			for (i=0; i<atomClean; i++) {
				k=listAtomClean[i];
				atomDefine[k]=0;
			};
		};
		currNumDef=0;
		i=0;
		while (currNumDef<atomClean) {
			//start recursion-priority list formation}
			defC(currNumDef,baseCycle,atomClean,cycleDefine,atomDefine,atomCycle,cycleAddress,cycleSize,dsATN,dsTP,dsSC,dsNA1,dsNA2);
			defA(currNumDef,atomClean,spn,baseCycle,atomDefine,listAtomClean,cycleDefine,cycleSize,cycleAddress,atomCycle,dsATN,dsTP,dsNA1,dsNA2);
			i++;
		}; //end recursion}

		for (i=0; i<nAtoms(); i++) atomDefine[i]=0;

		if ((spn<3) || (spn==4)) {     //checking already cleaned atoms...
			for (i=0; i<nAtoms(); i++) atomDefine[i]=1;
			for (i=0; i<atomClean; i++) {
				k=listAtomClean[i];
				atomDefine[k]=0;
			};
			for (i=0; i<nAtoms(); i++) if (atomDefine[i]==1) {
				getAtom(i)->rx=(xCenterOld+(getAtom(i)->rx-xCenterOld)/bondLengthOld);
				getAtom(i)->ry=(yCenterOld+(getAtom(i)->ry-yCenterOld)/bondLengthOld);
			};
		};

		i=0;

		if (atomClean>=1) while (i<atomClean) {

			switch (dsTP[i]) {
	case 0:{
		atomDefine[dsATN[i]]=1; //old coordinates of the atom are used}
		getAtom(dsATN[i])->rx=xCenterOld+(getAtom(dsATN[i])->rx-xCenterOld)/bondLengthOld;
		getAtom(dsATN[i])->ry=yCenterOld+(getAtom(dsATN[i])->ry-yCenterOld)/bondLengthOld;
		break;
			 }
	case 1:
	case 2: if (atomDefine[dsATN[i]] == 0) {
		k1=dsNA1[i];
		k=getAtom(k1)->nb;
		nb=0;
		ux=0;
		uy=0;
		ux1=0;
		uy1=0;
		//bond direction calculation}
		for (j=0; j<k; j++) if (atomDefine[getAtom(k1)->ac[j]]>0) {
			k2=getAtom(k1)->ac[j];
			nb=nb+1;           //number of already cleaned atoms}
			ux=ux+(getAtom(k1)->rx-getAtom(k2)->rx);
			uy=uy+(getAtom(k1)->ry-getAtom(k2)->ry);
			if (dsTP[i]==1) for (l=0; l<getAtom(k2)->nb; l++) if (atomDefine[getAtom(k2)->ac[l]]>0) {
				ux1=ux1+(getAtom(k2)->rx-getAtom(getAtom(k2)->ac[l])->rx);
				uy1=uy1+(getAtom(k2)->ry-getAtom(getAtom(k2)->ac[l])->rx);
			};
		};
		if (dsTP[i]==1) {     //chain fragment
			if ((abs(ux) <= 0.00001) && (abs(uy) <= 0.00001)) {
				bnEx=-1;
				for (j=0; j<bk[k1].nb; j++) {
					n=bk[k1].adjusted[j];
					n1=getBond(n)->at[0];
					if (n1 == k1) n1=getBond(n)->at[1];
					if (n1 == dsATN[i]) bnEx=n;
					if (bnEx >= 0) break;
				};
				test= (bnEx >= 0);
				if (test) test=threeBondResolve(k1,bnEx,ux,uy,&bk);
				if (! test) {
					ux=ux1; uy=uy1;
				};
			};

			if ((ux==0) && (uy==0)) ux=1;
			fi=(k-nb-1)*M_PI/(double)k;
			//Addition for fluorine chain fragments...
			isChainFour=false;
			n=0; mm1=-1; mm2=-1; mAny=-1;
			if ((k == 4) && (getAtom(k1)->na == 6))  for (j=0; j<getAtom(k1)->nb; j++) {
				n1=getAtom(k1)->ac[j];
				if (atomDefine[n1] == 0) {
					mAny=n1;
					if (getAtom(n1)->nb == 1) n++; else {
						if (mm1 == -1) mm1=n1; else mm2=n1;
					};
				};
			};
			if ((mm1 == -1) && (n == 3)) {
				mm1=mAny;
				n--;
			};
			isChainFour=(n == 2) && (mm1 >= 0);
			//End addition
			isCycle=false;
			if ((nb == 2) && ((k == 4) || (k == 5))) {
				n=0;
				for (j=0; j<bk[dsNA1[i]].nb; j++) {
					k2=bk[dsNA1[i]].adjusted[j];
					if (getBond(k2)->db > 1) {
						n1=getBond(k2)->at[0];
						if (n1 == k1) n1=getBond(k2)->at[1];
						if (atomDefine[n1] != 0) n++;
					};
				};
				isCycle=(n == 2);
			};
			if (isChainFour) {
				n=0;
				fi=M_PI/3;
				ux1= ux*cos(fi)+uy*sin(fi);
				uy1=-ux*sin(fi)+uy*cos(fi);
				cf=sqrt(ux1*ux1+uy1*uy1);
				if (cf != 0) {
					ux1=ux1/cf; uy1=uy1/cf;
				};
				getAtom(mm1)->rx=getAtom(k1)->rx+ux1; //coordinates
				getAtom(mm1)->ry=getAtom(k1)->ry+uy1;
				atomDefine[mm1]=1;
				if (mm2 >= 0) {
					fi=-M_PI/3;
					ux1= ux*cos(fi)+uy*sin(fi);
					uy1=-ux*sin(fi)+uy*cos(fi);
					cf=sqrt(ux1*ux1+uy1*uy1);
					if (cf != 0) {
						ux1=ux1/cf; uy1=uy1/cf;
					};
					getAtom(mm2)->rx=getAtom(k1)->rx+ux1; //coordinates
					getAtom(mm2)->ry=getAtom(k1)->ry+uy1; //coordinates
					atomDefine[mm2]=1;
				};
				ux=0; uy=0;  //New cosines....
				for (j=0; j<getAtom(k1)->nb; j++) if (atomDefine[getAtom(k1)->ac[j]] > 0) {
					k2=getAtom(k1)->ac[j];
					ux=ux+(getAtom(k1)->rx-getAtom(k2)->rx);
					uy=uy+(getAtom(k1)->ry-getAtom(k2)->ry);
				};
				fi=M_PI/6;
				for (j=0; j<getAtom(k1)->nb; j++) {
					n1=getAtom(k1)->ac[j];
					if (atomDefine[n1] == 0) {
						ux1= ux*cos(fi)+uy*sin(fi);
						uy1=-ux*sin(fi)+uy*cos(fi);
						cf=sqrt(ux1*ux1+uy1*uy1);
						if (cf != 0) {
							ux1=ux1/cf; uy1=uy1/cf;
						};
						getAtom(n1)->rx=getAtom(k1)->rx+ux1; //coordinates}
						getAtom(n1)->ry=getAtom(k1)->ry+uy1; //coordinates}
						fi=-fi;
						atomDefine[n1]=1;
					};
				};
			} else if (isCycle) {
				if (k == 4) fi=M_PI/6; else fi=M_PI/3;
				ux1= ux*cos(fi)+uy*sin(fi);
				uy1=-ux*sin(fi)+uy*cos(fi);
				cf=sqrt(ux1*ux1+uy1*uy1);
				if (cf != 0) {
					ux1=ux1/cf; uy1=uy1/cf;
				};
				getAtom(dsATN[i])->rx=getAtom(k1)->rx+ux1; //coordinates
				getAtom(dsATN[i])->ry=getAtom(k1)->ry+uy1; //coordinates
				atomDefine[dsATN[i]]=1;
				for (j=0; j<getAtom(dsNA1[i])->nb; j++) {
					k2=getAtom(dsNA1[i])->ac[j];
					if (atomDefine[k2] == 0) {
						ux1=ux*cos(fi)-uy*sin(fi);
						uy1=ux*sin(fi)+uy*cos(fi);
						if (cf != 0) {
							ux1=ux1/cf; uy1=uy1/cf;
						};
						getAtom(k2)->rx=getAtom(k1)->rx+ux1;
						getAtom(k2)->ry=getAtom(k1)->ry+uy1;
						fi=0;
						atomDefine[k2]=1;
					};
				};
			} else {

				if (bk[dsNA1[i]].nb == 2) {  //two-connected fragment
					//Search for triple bond or for two double-bonds
					n1=bk[dsNA1[i]].adjusted[0];
					n2=bk[dsNA1[i]].adjusted[1];
					if ((getBond(n1)->tb == 3) || (getBond(n2)->tb == 3)) k=3; else
						if ((getBond(n1)->tb == 2) && (getBond(n2)->tb == 2)) {
							if (getAtom(dsNA1[i])->na == 6) k=3;
						};
				};
				//Two-conncted fragment like pyrophosphate
				if ((k== 2) && (getAtom(dsNA1[i])->na != 6)) {
					n1=getAtom(dsNA1[i])->ac[0];
					n2=getAtom(dsNA1[i])->ac[1];
					if ((getAtom(n1)->nb >= 4) && (getAtom(n2)->nb >= 4) && (getAtom(n1)->na != 6) && (getAtom(n2)->na != 6)) k=3;
				};
				//End pyrophosphate
				if (k==2) { //120 degrees fragment}
					if ((ux1==0) && (uy1==0)) ux1=1;
					fi=M_PI/3;
					cf=uy*ux1-ux*uy1;
					if (cf != 0) fi=fi*cf/abs(cf);
				};
				ux1=ux*cos(fi)+uy*sin(fi);
				uy1=-ux*sin(fi)+uy*cos(fi);
				cf=sqrt(ux1*ux1+uy1*uy1);
				if (cf != 0) {
					ux1=ux1/cf; uy1=uy1/cf;
				};
				getAtom(dsATN[i])->rx=getAtom(k1)->rx+ux1; //coordinates
				getAtom(dsATN[i])->ry=getAtom(k1)->ry+uy1; //coordinates
				atomDefine[dsATN[i]]=1;
			};
		} else {                   //cyclic fragment}

			if ((ux==0) && (uy==0)) ux=1;
			fi=(k-nb-2)*M_PI/(double)k;
			ux1=ux*cos(fi)+uy*sin(fi);
			uy1=-ux*sin(fi)+uy*cos(fi);
			cf=sqrt(ux1*ux1+uy1*uy1);
			if (cf != 0) {
				ux1=ux1/cf; uy1=uy1/cf;
			};
			nb=dsSC[i];
			cs=nb+1;
			cf=1/(2*sin(M_PI/(double)cs));
			ux=getAtom(dsNA1[i])->rx+cf*ux1;
			uy=getAtom(dsNA1[i])->ry+cf*uy1;
			fi=2*M_PI/(double)cs;
			ux1=getAtom(dsNA1[i])->rx-ux;
			uy1=getAtom(dsNA1[i])->ry-uy;
			//coordinates of all atoms for the cycle under study are calculated}

			for (j=0; j<nb; j++) {
				getAtom(dsATN[i+j-0])->rx=ux+ux1*cos((j+1)*fi)+uy1*sin((j+1)*fi);
				getAtom(dsATN[i+j-0])->ry=uy-ux1*sin((j+1)*fi)+uy1*cos((j+1)*fi);
				atomDefine[dsATN[i+j-0]]=1;
			};

			i=i+nb-1;
		};
		break;
			  } else break;
	case 3:
	case 4: {

		twoAtomUnitVector(dsNA1[i],dsNA2[i],ux,uy,atomDefine);
		//calculation of an optimal side to add new fragment
		ux1=getAtom(dsNA1[i])->rx-getAtom(dsNA2[i])->rx;
		uy1=getAtom(dsNA1[i])->ry-getAtom(dsNA2[i])->ry;
		if (dsTP[i]==3)   {  //angle calc. for condenced cycle}
			nb=dsSC[i];
			cs=nb+2;
			fi=sqrt(ux1*ux1+uy1*uy1);
			cf=fi/(2*sin(M_PI/(double)cs)/cos(M_PI/(double)cs));
			add=0;
		} else {     //angle calc. for polycycle
			r=sqrt(ux1*ux1+uy1*uy1);
			add=(int)ceil(r);
			nb=dsSC[i];
			cs=nb+2+add;
			cf=r/(2*sin(M_PI*(nb+1)/(double)cs));
			r=r/2.0;
			cf=sqrt(cf*cf-r*r);
		};
		ux=(getAtom(dsNA1[i])->rx+getAtom(dsNA2[i])->rx)/2.0+ux*cf;
		uy=(getAtom(dsNA1[i])->ry+getAtom(dsNA2[i])->ry)/2.0+uy*cf;
		ux1=getAtom(dsNA1[i])->rx-ux;
		uy1=getAtom(dsNA1[i])->ry-uy;
		fi=2*M_PI/(double)cs;
		ux2=ux+ux1*cos((1.0+add)*fi)+uy1*sin((1.0+add)*fi);
		uy2=uy-ux1*sin((1.0+add)*fi)+uy1*cos((1.0+add)*fi);
		if ((abs(getAtom(dsNA2[i])->rx-ux2)<0.01) && (abs(getAtom(dsNA2[i])->ry-uy2)<0.01)) fi=-fi;
		//FAtom's coordinates calculation
		for (j=0; j<nb; j++) {
			getAtom(dsATN[i+j-0])->rx=ux+ux1*cos((j+1)*fi)+uy1*sin((j+1)*fi);
			getAtom(dsATN[i+j-0])->ry=uy-ux1*sin((j+1)*fi)+uy1*cos((j+1)*fi);
			atomDefine[dsATN[i+j-0]]=1;
		}
		i=i+nb-1;
		break;
			  }
			}
			i++;
		}
		bondLengthNew=0;

		//Rescaling and shift of structure
		for (i=0; i<bondClean; i++)  bondLengthNew=bondLengthNew+bondLength(listBondClean[i]);
		bondLengthNew=bondLengthNew/(double)bondClean;
		xCenterNew=0; yCenterNew=0;
		for (i=0; i<atomClean; i++) {
			xCenterNew=xCenterNew+getAtom(listAtomClean[i])->rx;
			yCenterNew=yCenterNew+getAtom(listAtomClean[i])->ry;
		};
		xCenterNew=xCenterNew/(double)atomClean; yCenterNew=yCenterNew/(double)atomClean;
		bondLengthNew=0.15*bondLengthOld;
		//IOPT[7]-controlles, whether or not the atom's shifts should be created, if coordinates of pair of atoms are identical
		if (spn==3) {                     //Rescaling and shift coordinates for group (CODE=3)
			xCenterNew=getAtom(sCHA1)->rx;
			yCenterNew=getAtom(sCHA1)->ry;
			ux=getAtom(atomSecond)->rx-getAtom(sCHA1)->rx;
			uy=getAtom(atomSecond)->ry-getAtom(sCHA1)->ry;
			r=sqrt(ux*ux+uy*uy);
			ux=ux/r; uy=uy/r;
			c=ux*uvX+uy*uvY;
			s=ux*uvY-uy*uvX;
			for (i=0; i<atomClean; i++) {
				ux1=getAtom(listAtomClean[i])->rx-xCenterNew;
				uy1=getAtom(listAtomClean[i])->ry-yCenterNew;
				getAtom(listAtomClean[i])->rx=c*ux1-s*uy1;
				getAtom(listAtomClean[i])->ry=s*ux1+c*uy1;
			};
			xCenterNew=0;
			yCenterNew=0;

		};
		if (spn != 4) for (i=0; i<atomClean; i++) {   //New screen coordinates
			getAtom(listAtomClean[i])->rx=getAtom(listAtomClean[i])->rx-xCenterNew+xCenterOld;
			getAtom(listAtomClean[i])->ry=getAtom(listAtomClean[i])->ry-yCenterNew+yCenterOld;
		};
};


int TSimpleMolecule::loadAllFromStream(std::istream & data) {
  int result=0;

  return result;
};

/*
function TSimpleMolecule.LoadAllFromStream(S:TStream):integer;
var
  N,NW:word;
  I:integer;
  SA:TSingleAtom;
  SB:TSingleBond;
{$IFNDEF NIST}
  R1,R2,R3:single;
{$ENDIF}
begin
  MoleculeZero;
  result:=0;
  S.Read(N,sizeof(N));
  FAtom.ReDim(N);
  result:=result+sizeof(N);
  S.Read(N,sizeof(N));
  FBond.ReDim(N);
  result:=result+sizeof(N);
  S.Read(NW,sizeof(NW));
  result:=result+sizeof(N);
  for I:=1 to FAtom.MaxIndex do begin
    S.Read(SA,sizeof(SA));
    FAtom[I]:=SA;
    result:=result+sizeof(SA);
  end;
  for I:=1 to FBond.MaxIndex do begin
    S.Read(SB,sizeof(SB));
    FBond[I]:=SB;
    result:=result+sizeof(SB);
  end;
  S.Read(FUnC1,sizeof(FUnC1));
  result:=result+sizeof(FUnC1);
  S.Read(FUnC2,sizeof(FUnC2));
  result:=result+sizeof(FUnC1);
  S.Read(FMolweight,sizeof(FMolweight));
  result:=result+sizeof(FMolweight);
{$IFDEF NIST}
  S.read(FXCoorShift,sizeof(FXCoorShift));
  result:=result+sizeof(FXCoorShift);
  S.Read(FYCoorShift,sizeof(FYCoorShift));
  result:=result+sizeof(FYCoorShift);
  S.Read(FCoorScale,sizeof(FCoorScale));
  result:=result+sizeof(FCoorScale);
{$ELSE}
  //rescaling if necessary
  S.read(R1,sizeof(R1));
  result:=result+sizeof(R1);
  S.read(R2,sizeof(R2));
  result:=result+sizeof(R2);
  S.read(R3,sizeof(R3));
  result:=result+sizeof(R3);
  if (R1<>0) or (R2<>0) or (R3<>1) then for I:=1 to FAtom.NAtoms do begin
    FAtom.RX[I]:=FAtom.RX[I]*R3+R1;
    FAtom.RY[I]:=FAtom.RY[I]*R3+R2;
  end;
{$ENDIF}
end;
*/


bool TSimpleMolecule::readMolfile(const std::string fileName) {
	bool result = false;
	std::ifstream fileMol(fileName);
	if (fileMol.is_open()) result=readMolfile(fileMol);
	return result;

}

void TSimpleMolecule::saveMolfile(const std::string fileName) {
	std::ofstream fileMol(fileName);
	if (fileMol.is_open()) {
		getMolfile(fileMol);
		fileMol.close();
	}
	
}

bool TSimpleMolecule::readMolfile(std::istream & data) {
  TSingleAtom * sa;
  TSingleBond * sb;
  std::string s,s1;
  int i,j,na,nb,nv,nc,iz,anum;
  double rx,ry,rz;
  bool result=false;
  int at1,at2,tb,n,n1,n2;
  bool test;  


  this->clear();
  getline(data,s);
  if (! data.good()) return result;
  getline(data,s);
  if (! data.good()) return result;
  getline(data,s);
  if (! data.good()) return result;
  getline(data,s);
  if (! data.good()) return result;
  //Number of atoms
  s1=s.substr(0,3);
  s=s.substr(3);
  na=atoi(s1.c_str());
  //Number of bonds
  s1=s.substr(0,3);
  s=s.substr(3);
  nb=atoi(s1.c_str());

  for (i=0; i<na; i++) {
    getline(data,s);
	s1=s.substr(0,10);
	rx=atof(s1.c_str());
	s1=s.substr(10,10);
	ry=atof(s1.c_str());
	s1=s.substr(20,10);
	rz=atof(s1.c_str());
	s1=s.substr(30,3);
	//removing spaces
	if (s1.at(0) == ' ') s1=s1.substr(1);
	if (s1.at(0) == ' ') s1=s1.substr(1);
	if (s1.at(s1.length()-1) == ' ') s1=s1.substr(0,s1.length()-1);
	if (s1.at(s1.length()-1) == ' ') s1=s1.substr(0,s1.length()-1);
    anum=0;
	for (j=1; j<NELEMMCDL; j++) if (s1 == aSymb[j]) {
	  anum=j;
	  break;
	};
	if (anum == 0) return result;
	s1=s.substr(33,3);
	iz=atoi(s1.c_str());
	s1=s.substr(36,3);
	nc=atoi(s1.c_str());
	s1=s.substr(48,3);
	nv=atoi(s1.c_str());
	
    sa=new TSingleAtom();
	sa->na=anum;
	if (nc != 0) {
	  if (nc == 4) sa->rl=1; else if (nc == 9) sa->rl=2; else sa->nc=4-nc;
	}
	sa->iz=iz;
	if (nv == 0) sa->nv=hVal[anum]; else {
      if (nv != 15) sa->nv=nv; else sa->nv=0;  //=15 - zero
	}
	sa->rx=rx;
	sa->ry=ry;

	fAtom.push_back(sa);
  };
  for (i=0; i<nb; i++) {
    getline(data,s);

	s1=s.substr(0,3);
	at1=atoi(s1.c_str());

	s1=s.substr(3,3);
	at2=atoi(s1.c_str());

	s1=s.substr(6,3);
	tb=atoi(s1.c_str());

	s1=s.substr(9,3);
	n=atoi(s1.c_str());
	if (n == 1) tb=9; else
	if (n == 6) tb=10; else
	if (n == 4) tb=11;

	if (s.size() >= 15) {
	  s1=s.substr(15,3);
	  n=atoi(s1.c_str());
	} else n=0;

    sb=new TSingleBond();
	if (n == 1) sb->special=2; else if (n == 2) sb->special=1;
	sb->at[0]=at1-1;
	sb->at[1]=at2-1;
	sb->tb=tb;



	fBond.push_back(sb);
  }
  test=true;
  while (test && data.good()) {
	getline(data,s);
	test=! ((s.find("M") != std::string::npos) && (s.find("END") != std::string::npos));
	if (test) {
	  if (s.find("CHG") != std::string::npos) {
        n=s.find("CHG");
		s=s.substr(n+3);
		s=removeSpacesTab(s);
		n=s.find(" ");
		s1=s.substr(0,n);
		s=s.substr(n+1);
		s=removeSpacesTab(s);  
		while (s.size() > 0) {
		  n=s.find(" ");
		  if (n > 0) {
			s1=s.substr(0,n);
			s=s.substr(n+1);
			n1=atoi(s1.c_str());
			s=removeSpacesTab(s);
		    n=s.find(" ");
		    if (n > 0) {
			  s1=s.substr(0,n);
			  s=s.substr(n+1);
			} else {
			  s1=s;
			  s="";
			};
			n2=atoi(s1.c_str());
			if ((n1 > 0) && (n1 <= fAtom.size())) {
			  getAtom(n1-1)->nc=n2;
			  if (getAtom(n1-1)->nv < hVal[getAtom(n1-1)->na]) {
                getAtom(n1-1)->nv=getAtom(n1-1)->nv+abs(getAtom(n1-1)->nc);
				if (getAtom(n1-1)->nv > hVal[getAtom(n1-1)->na]) getAtom(n1-1)->nv=hVal[getAtom(n1-1)->na];
			  };
			};
		  } else s=""; 
		};
	  } else if (s.find("RAD") != std::string::npos) {
        n=s.find("RAD");
		s=s.substr(n+3);
		s=removeSpacesTab(s);
		n=s.find(" ");
		s1=s.substr(0,n);
		s=s.substr(n+1);
		s=removeSpacesTab(s);  
		while (s.size() > 0) {
		  n=s.find(" ");
		  if (n > 0) {
			s1=s.substr(0,n);
			s=s.substr(n+1);
			n1=atoi(s1.c_str());
			s=removeSpacesTab(s);
		    n=s.find(" ");
		    if (n > 0) {
			  s1=s.substr(0,n);
			  s=s.substr(n+1);
			} else {
			  s1=s;
			  s="";
			};
			n2=atoi(s1.c_str());
			if ((n1 > 0) && (n1 <= fAtom.size()) && (n2 >= 1) && (n2 <= 3)) {
			  if (n2 == 2) getAtom(n1-1)->rl=1; else getAtom(n1-1)->rl=2;
			  if (getAtom(n1-1)->nv < hVal[getAtom(n1-1)->na]) {
                getAtom(n1-1)->nv=getAtom(n1-1)->nv+getAtom(n1-1)->rl;
				if (getAtom(n1-1)->nv > hVal[getAtom(n1-1)->na]) getAtom(n1-1)->nv=hVal[getAtom(n1-1)->na];
			  };
			};
		  } else s=""; 
		};
	  };
	};
  };

  result=true;
  defineAtomConn();
  allAboutCycles();
  return result;
};

bool TSimpleMolecule::readMolfile(const std::vector<std::string> & data) {
	TSingleAtom * sa;
	TSingleBond * sb;
	std::string s, s1;
	int i, j, na, nb, nv, nc, iz, anum, lineStart;
	double rx, ry, rz;
	bool result = false;
	int at1, at2, tb, n, n1, n2;
	bool test;


	this->clear();
	if (data.size() < 2) return result;
	//search for line with data

	lineStart = -1;
	n = 6;
	if (n > data.size()) n = data.size();

	for (i = 0; i < n; i++) {
		s = data[i];
		s = s.substr(0, 3);
		na = atoi(s.c_str());
		if (na > 0) {
			lineStart = i;
			break;
		}
	}
	if (lineStart < 0) return result;

	s = data[lineStart];
	//Number of atoms
	s1 = s.substr(0, 3);
	s = s.substr(3);
	na = atoi(s1.c_str());
	//Number of bonds
	s1 = s.substr(0, 3);
	s = s.substr(3);
	nb = atoi(s1.c_str());

	for (i = 0; i<na; i++) {
		s = data[lineStart+1 + i];
		s1 = s.substr(0, 10);
		rx = atof(s1.c_str());
		s1 = s.substr(10, 10);
		ry = atof(s1.c_str());
		s1 = s.substr(20, 10);
		rz = atof(s1.c_str());
		s1 = s.substr(30, 3);
		//removing spaces
		if (s1.at(0) == ' ') s1 = s1.substr(1);
		if (s1.at(0) == ' ') s1 = s1.substr(1);
		if (s1.at(s1.length() - 1) == ' ') s1 = s1.substr(0, s1.length() - 1);
		if (s1.at(s1.length() - 1) == ' ') s1 = s1.substr(0, s1.length() - 1);
		anum = 0;
		for (j = 1; j<NELEMMCDL; j++) if (s1 == aSymb[j]) {
			anum = j;
			break;
		};
		if (anum == 0) return result;
		s1 = s.substr(33, 3);
		iz = atoi(s1.c_str());
		s1 = s.substr(36, 3);
		nc = atoi(s1.c_str());
		s1 = s.substr(48, 3);
		nv = atoi(s1.c_str());

		sa = new TSingleAtom();
		sa->na = anum;
		if (nc != 0) {
			if (nc == 4) sa->rl = 1; else if (nc == 9) sa->rl = 2; else sa->nc = 4 - nc;
		}
		sa->iz = iz;
		if (nv == 0) sa->nv = hVal[anum]; else {
			if (nv != 15) sa->nv = nv; else sa->nv = 0;  //=15 - zero
		}
		sa->rx = rx;
		sa->ry = ry;

		fAtom.push_back(sa);
	};
	for (i = 0; i<nb; i++) {
		s = data[na + lineStart+1 + i];
//		getline(data, s);

		s1 = s.substr(0, 3);
		at1 = atoi(s1.c_str());

		s1 = s.substr(3, 3);
		at2 = atoi(s1.c_str());

		s1 = s.substr(6, 3);
		tb = atoi(s1.c_str());

		s1 = s.substr(9, 3);
		n = atoi(s1.c_str());
		if (n == 1) tb = 9; else
			if (n == 6) tb = 10; else
				if (n == 4) tb = 11;

		if (s.size() >= 15) {
			s1 = s.substr(15, 3);
			n = atoi(s1.c_str());
		}
		else n = 0;

		sb = new TSingleBond();
		if (n == 1) sb->special = 2; else if (n == 2) sb->special = 1;
		sb->at[0] = at1 - 1;
		sb->at[1] = at2 - 1;
		sb->tb = tb;



		fBond.push_back(sb);
	}
	/*
	test = true;
	while (test && data.good()) {
		getline(data, s);
		test = !((s.find("M") != std::string::npos) && (s.find("END") != std::string::npos));
		if (test) {
			if (s.find("CHG") != std::string::npos) {
				n = s.find("CHG");
				s = s.substr(n + 3);
				s = removeSpacesTab(s);
				n = s.find(" ");
				s1 = s.substr(0, n);
				s = s.substr(n + 1);
				s = removeSpacesTab(s);
				while (s.size() > 0) {
					n = s.find(" ");
					if (n > 0) {
						s1 = s.substr(0, n);
						s = s.substr(n + 1);
						n1 = atoi(s1.c_str());
						s = removeSpacesTab(s);
						n = s.find(" ");
						if (n > 0) {
							s1 = s.substr(0, n);
							s = s.substr(n + 1);
						}
						else {
							s1 = s;
							s = "";
						};
						n2 = atoi(s1.c_str());
						if ((n1 > 0) && (n1 <= fAtom.size())) {
							getAtom(n1 - 1)->nc = n2;
							if (getAtom(n1 - 1)->nv < hVal[getAtom(n1 - 1)->na]) {
								getAtom(n1 - 1)->nv = getAtom(n1 - 1)->nv + abs(getAtom(n1 - 1)->nc);
								if (getAtom(n1 - 1)->nv > hVal[getAtom(n1 - 1)->na]) getAtom(n1 - 1)->nv = hVal[getAtom(n1 - 1)->na];
							};
						};
					}
					else s = "";
				};
			}
			else if (s.find("RAD") != std::string::npos) {
				n = s.find("RAD");
				s = s.substr(n + 3);
				s = removeSpacesTab(s);
				n = s.find(" ");
				s1 = s.substr(0, n);
				s = s.substr(n + 1);
				s = removeSpacesTab(s);
				while (s.size() > 0) {
					n = s.find(" ");
					if (n > 0) {
						s1 = s.substr(0, n);
						s = s.substr(n + 1);
						n1 = atoi(s1.c_str());
						s = removeSpacesTab(s);
						n = s.find(" ");
						if (n > 0) {
							s1 = s.substr(0, n);
							s = s.substr(n + 1);
						}
						else {
							s1 = s;
							s = "";
						};
						n2 = atoi(s1.c_str());
						if ((n1 > 0) && (n1 <= fAtom.size()) && (n2 >= 1) && (n2 <= 3)) {
							if (n2 == 2) getAtom(n1 - 1)->rl = 1; else getAtom(n1 - 1)->rl = 2;
							if (getAtom(n1 - 1)->nv < hVal[getAtom(n1 - 1)->na]) {
								getAtom(n1 - 1)->nv = getAtom(n1 - 1)->nv + getAtom(n1 - 1)->rl;
								if (getAtom(n1 - 1)->nv > hVal[getAtom(n1 - 1)->na]) getAtom(n1 - 1)->nv = hVal[getAtom(n1 - 1)->na];
							};
						};
					}
					else s = "";
				};
			};
		};
	};
	*/
	result = true;
	defineAtomConn();
	allAboutCycles();
	return result;
};


/*
std::string intToStr(int value, int length){
  std::string result="";
  char buffer [16];

  itoa(value,buffer,10);
  result=buffer;
  while (result.size() < length) result=" "+result;
  return result;
};


std::string floatToStr(double value, int length, int afterPoint){
  std::string result="";
  char buffer [16];

  //ftoa(value,buffer,afterPoint);
  result=printf("f10.4",value);
  //result=buffer;
  //while (result.size() < length) result=" "+result;
  return result;
};

*/

std::string TSimpleMolecule::getSimplePolymerMolBlock() const {
	const std::string endln = "\n";
	const std::string SAL = "M  SAL   1"; //atoms inside PRU, wirhout *   'M  SAL   1  5   2   3   4   5   6'
	const std::string SBL = "M  SBL   1  2"; //Crossed bond numbers //'M  SBL   1  2   1   6'
	const std::string HEAD = "M  STY  1   1 SRU" + endln+ "M  SLB  1   1   1" + endln + "M  SCN  1   1 HT" + endln;
	const std::string SDI = "M  SDI   1  4";//    5.6200   -1.0800    4.7900   -1.0800';
	int i, n, an1, an2, bn1, bn2, nSal;
	std::string s, ss;
	std::string result = "";

	n = 0;
	an1 = -1; an2 = -1;
	for (i = 0; i < nAtoms(); i++) if ((getAtom(i)->na == ID_ZVEZDA) && (getAtom(i)->nb == 1)) {
		n++;
		if (an1 < 0)  an1 = i; else an2 = i;
	};
	if (n != 2) return result;
	bn1 = -1; bn2 = -1;
	for (i = 0; i < nBonds(); i++) if ((getBond(i)->at[0] == an1) || (getBond(i)->at[1] == an1) || (getBond(i)->at[0] == an2) || (getBond(i)->at[1] == an2)) {
		if (bn1 < 0) bn1 = i; else bn2 = i;
	};
	

	result = HEAD;
	nSal = 0;
	ss = "";
	for (i = 0; i < nAtoms(); i++) {
		if ((i != an1) && (i != an2)) {
			nSal++;
			s = std::to_string(i+1);
			while (s.length() < 4) s = " " + s;
			ss = ss + s;
		};
		
		if ((((nSal % 15) == 0) || (i == (nAtoms()-1))) && (ss.length() > 0)) {
			n = nSal % 15;
			if (n == 0) n = 15;
			s=std::to_string(n);
			while (s.length() < 3) s = " " + s;
			ss = SAL + s + ss;
			result = result + ss + endln;
			ss = "";
			s = "";
		};
		
	};
	s = std::to_string(bn1+1);
	while (s.length() < 4) s = " " + s;
	result = result + SBL + s;
	s = std::to_string(bn2+1);
	while (s.length() < 4) s = " " + s;
	result = result + s;

	return result;
};

std::string TSimpleMolecule::getMolfile() const {
	std::string result = "";
	ostringstream ssData;
	getMolfile(ssData, false, NULL);
	result = ssData.str();
	return result;
};



void TSimpleMolecule::getMolfile(std::ostream & data, bool writeend, std::vector<double> * zCoor) const {
	//char buff[BUFF_SIZE];
	const TSingleAtom * sa;
	const TSingleBond * sb;
	int charge,bondType,stereoType,nv,k,i;
	std::string asym;
    std::vector<int> sChg,sRad;



	data<<endl<<endl<<endl;  //three empty lines
	//data<< intToStr(nAtoms(),3) << intToStr(nBonds(),3) << intToStr(0,3) << intToStr(0,3) << intToStr(0,3) << intToStr(0,3) << intToStr(0,3) << intToStr(0,3) << intToStr(0,3) << intToStr(0,3) << intToStr(999,3) << " V2000";
	data << std::setw(3) << nAtoms() << std::setw(3) << nBonds() << std::setw(3) << 0 << std::setw(3) << 0 << std::setw(3) << 0;
	data << std::setw(3) << 0 << std::setw(3) << 0 << std::setw(3) << 0 << std::setw(3) << 0 << std::setw(3) << 0 << std::setw(3) << 999 << " V2000" << endl;
	//_sprintf(buff,BUFF_SIZE,"%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d V2000",nAtoms(),nBonds(),0,0,0,0,0,0,0,0,999);
	//data << buff << endl;
	for (int i=0; i<nAtoms(); i++) {
		sa=getAtom(i);
		switch (sa->nc) {
	      case 1: charge = 3; break;
	      case 2: charge = 2; break;
	      case 3: charge = 1; break;
	      case -1: charge = 5; break;
	      case -2: charge = 6; break;
	      case -3: charge = 7; break;
	      default: charge=0; break;
		};
		if (sa->rl == 1) charge=4; else if (sa->rl == 2) charge=9;
		if (sa->nc != 0) {
		  sChg.push_back(i+1);
		  sChg.push_back(sa->nc);
		};
		if (sa->rl != 0) {
		  sRad.push_back(i+1);
		  if (sa->rl == 1) k=2; else k=3;
		  sRad.push_back(k);
		};

		nv=0;
        if (getAtom(i)->nv == hVal[getAtom(i)->na]) nv=0; else {
          nv=getAtom(i)->nv;
          if (nv == 0) nv=15; //Addition from February 2006
		};
        if ((getAtom(i)->na == 7) && (getAtom(i)->nv > 3) && (getAtom(i)->nc == 1)) nv=0;  //N
        if ((getAtom(i)->na == 5) && (getAtom(i)->nv > 3) && (getAtom(i)->nc == -1)) nv=0; //B

		double rx=sa->rx;
		rx=(floor(rx*100000)/100000);
		double ry=sa->ry;
		ry=(floor(ry*100000)/100000);
		double rz = 0.0;
		if (zCoor) rz = (*zCoor)[i];

		if (abs(rx)<0.000001) rx=0.0;
		if (abs(ry)<0.000001) ry=0.0;

		data << std::setw(10) << std::setprecision(7) << rx << std::setw(10) << std::setprecision(7) << ry << std::setw(10) << std::setprecision(7) << rz;
		asym=aSymb[sa->na];
		if (asym.size() == 1) asym=asym+" ";
		data << std::setw(3) << asym << std::setw(3) << sa->iz << std::setw(3) << charge << std::setw(3) << 0 << std::setw(3) << 0 << std::setw(3) << 0 << std::setw(3) << nv;
		data << 0 << std::setw(3) << 0 << std::setw(3) << 0 << std::setw(3) << 0 << endl;

		//_sprintf(buff, BUFF_SIZE, "%10.4f%10.4f%10.4f %-3s%2d%3d%3d%3d%3d",	sa->rx, sa->ry, 0.0, (aSymb[sa->na]).c_str(), 0,charge,0,0,0);    
		//data << buff << endl;
		//data << endl;
	};
	for (int i=0; i<nBonds(); i++) {
		sb=getBond(i);
		bondType=sb->tb;
		stereoType=0;
		if (bondType == 9) {
			bondType=1;
			stereoType=1;
		} else if (bondType == 10) {
			bondType=1;
			stereoType=6;
		} else if (bondType == 11) {
			bondType=1;
			stereoType=4;
		};
		data << std::setw(3) << sb->at[0]+1 << std::setw(3) << sb->at[1]+1 << std::setw(3) << bondType << std::setw(3) << stereoType << std::setw(3) << 0 << std::setw(3) << 0 <<endl;
		//_sprintf(buff, BUFF_SIZE, "%3d%3d%3d%3d%3d%3d",(sb->at[0]+1), (sb->at[1]+1), bondType, stereoType, 0, 0);
		//data << buff << endl;
	};
	if (sChg.size() > 0) {
	  data << "M  CHG" << std::setw(3) << sChg.size();
	  for (i=0; i<sChg.size(); i++) data << std::setw(3) << sChg[i];
	  data << endl;
	};
	if (sRad.size() > 0) {
	  data << "M  RAD" << std::setw(3) << sRad.size();
	  for (i=0; i<sRad.size(); i++) data << std::setw(3) << sRad[i];
	  data << endl;
	};
	std::string s = getSimplePolymerMolBlock();
	if (s.length() > 0) data << s << endl;
	data << "M  END" << endl;
	if (writeend) {
	  data << "$$$$" << endl;
	};
};

//*******************************************************************************
//Below routines from ChainRotate-flipping acyclic bonds and enlarging bonds to get fine picture
//*******************************************************************************

double TSimpleMolecule::atomDistanceMetric(int an) {
	int i, n;
	double r, rr;
	double x1, y1, x2, y2, d, result;

	if (getAtom(an)->nb == 0) return 0;
	result=0.01;
	n=getAtom(an)->ac[0];
	x1=getAtom(an)->rx-getAtom(n)->rx;
	y1=getAtom(an)->ry-getAtom(n)->ry;
	rr=sqrt(x1*x1+y1*y1);
	for (i=0; i<nAtoms(); i++) if ((i != an) && (i != n)) {
		x2=getAtom(i)->rx-getAtom(n)->rx;
		y2=getAtom(i)->ry-getAtom(n)->ry;
		d=rr*sqrt(x2*x2+y2*y2);
		if (d == 0) r=0; else r=(x1*x2+y1*y2)/d;
		if (r > 0) result=result+r;
	};
	return result;
};

bool TSimpleMolecule::bondsOverlapped(int bN1, int bN2, double delta) {
	double x1A,y1A,x2A,y2A,x1B,y1B,x2B,y2B;
	bool result;

	x1A=getAtom(getBond(bN1)->at[0])->rx;
	y1A=getAtom(getBond(bN1)->at[0])->ry;
	x2A=getAtom(getBond(bN1)->at[1])->rx;
	y2A=getAtom(getBond(bN1)->at[1])->ry;

	x1B=getAtom(getBond(bN2)->at[0])->rx;
	y1B=getAtom(getBond(bN2)->at[0])->ry;
	x2B=getAtom(getBond(bN2)->at[1])->rx;
	y2B=getAtom(getBond(bN2)->at[1])->ry;

	result=false;
	if (((x1A>(x1B+delta)) && (x2A>(x1B+delta)) && (x1A>(x2B+delta)) && (x2A>(x2B+delta))) ||
		((x1A<(x1B-delta)) && (x2A<(x1B-delta)) && (x1A<(x2B-delta)) && (x2A<(x2B-delta))) ||
		((y1A>(y1B+delta)) && (y2A>(y1B+delta)) && (y1A>(y2B+delta)) && (y2A>(y2B+delta))) ||
		((y1A<(y1B-delta)) && (y2A<(y1B-delta)) && (y1A<(y2B-delta)) && (y2A<(y2B-delta)))) return result;
	result=overlapped(x1A,y1A,x2A,y2A,x1B,y1B,x2B,y2B,delta);
	return result;
};


int TSimpleMolecule::hasOverlapped(double delta, bool findFirst) {
	int i, j;
	bool test;
	double r,xx,yy;
	int result=0;

	for (i=1; i<(nBonds()-1); i++) for (j=i+1; j<nBonds(); j++) {
		test=(getBond(i)->at[0] == getBond(j)->at[0]) || (getBond(i)->at[0] == getBond(j)->at[1]) ||
			(getBond(i)->at[1] == getBond(j)->at[0]) || (getBond(i)->at[1] == getBond(j)->at[1]);
		if (! test) if (bondsOverlapped(i,j,delta)) {
			result++;
			if (findFirst) return result;
		};
	};
	for (i=0; i<(nAtoms()-1); i++) for (j=i+1; j<nAtoms(); j++) {
		xx=getAtom(i)->rx-getAtom(j)->rx;
		yy=getAtom(i)->ry-getAtom(j)->ry;
		r=sqrt(xx*xx+yy*yy);
		if (r < (2*delta)) {
			result++;
			if (findFirst) return result;
		};
	};
	return result;
};

bool TSimpleMolecule::checkOverlapped() {
	bool result;
	result=hasOverlapped(averageBondLength()/32,true)>0;
	return result;
};

double TSimpleMolecule::averageBondLength() const {
	double result=0;
	if (nBonds() == 0) return result;
	for (int i=0; i<nBonds(); i++) result=result+bondLength(i);
	result=result/(double)nBonds();
	return result;
};

double TSimpleMolecule::averageAtomDistance() {
  int i,j,n;
  double result=0;
  double d;
  TSingleAtom * sa1;
  TSingleAtom * sa2;

  if (nAtoms() <= 1) return result;
  n=0;
  for (i=0; i<(nAtoms()-1); i++) for (j=i+1; j<nAtoms(); j++)  {
	sa1=getAtom(i);
    sa2=getAtom(j);
	d=sqrt((sa1->rx-sa2->rx)*(sa1->rx-sa2->rx)+(sa1->ry-sa2->ry)*(sa1->ry-sa2->ry));
	result=result+d;
	n++;
  };
  result=result/n;
  return result;
};


void TSimpleMolecule::flipSmall(int cHB) {
	//the bond number CHB is used to divide the fragment into two parts. For lower
	//part flip is executed. Case CHB-ring bond-no effects}

	int cHA1, cHA2, n;
	bool test;
	std::vector<int>list1(listarSize());
	double r, xc, yc, xo, yo, xn, yn;
	int i;
	int nB1;

	if (cHB < 0) return;

	test=makeFragment(nB1,list1,getBond(cHB)->at[1],getBond(cHB)->at[0]);
	if (nB1 > 1) {      
		//One of the atoms haven't neighbours-flip unavalable
		cHA1=getBond(cHB)->at[0];
		cHA2=getBond(cHB)->at[1];
		xc=getAtom(cHA2)->rx-getAtom(cHA1)->rx; //Axes direction calculation
		yc=getAtom(cHA2)->ry-getAtom(cHA1)->ry;
		r=sqrt(xc*xc+yc*yc);
		xc=xc/r;
		yc=yc/r;
		xo=xc*xc-yc*yc;
		yo=2*xc*yc;

		for (i=0; i<nB1; i++) {
			//rotation at angle Pi around axes for all atoms in the LIST fragment
			n=list1[i];
			xc=getAtom(n)->rx-getAtom(cHA1)->rx;
			yc=getAtom(n)->ry-getAtom(cHA1)->ry;
			xn=xc*xo+yc*yo;
			yn=xc*yo-yc*xo;
			getAtom(n)->rx=getAtom(cHA1)->rx+xn;
			getAtom(n)->ry=getAtom(cHA1)->ry+yn;
		};
	};
};

void TSimpleMolecule::bondEnlarge(int bN) {
	int nB;
	double r;
	std::vector<int> list(listarSize());
	int cHA1, cH1, cH2, n;
	bool test;
	double xc, yc, xc1, yc1;
	int i;

	for (i=0; i<nAtoms(); i++) list[i]=i;
	cHA1=getBond(bN)->at[0];
	test=makeFragment(nB,list,cHA1,getBond(bN)->at[1]);
	if (list[0] == getBond(bN)->at[0]) { //center definition
		cH1=getBond(bN)->at[0];
		cH2=getBond(bN)->at[1];
	} else {
		cH1=getBond(bN)->at[1];
		cH2=getBond(bN)->at[2];
	};
	xc=getAtom(cH1)->rx-getAtom(cH2)->rx;
	yc=getAtom(cH1)->ry-getAtom(cH2)->ry;
	r=sqrt(xc*xc+yc*yc);
	xc=xc/r;
	yc=yc/r;
	xc1=getAtom(cH2)->rx-getAtom(cH1)->rx;
	yc1=getAtom(cH2)->ry-getAtom(cH1)->ry;
	r=r*2;
	for (i=0; i<nB; i++) { //coordinates change
		n=list[i];
		getAtom(n)->rx=getAtom(n)->rx+xc1+r*xc;
		getAtom(n)->ry=getAtom(n)->ry+yc1+r*yc;
	};
};

int TSimpleMolecule::fragmentSecond(int sphere, int att, int secAt, const std::vector<int> a, 
	const std::vector<int> b, const std::vector<adjustedlist> & bk, std::vector<int>& wSphere) {

		/*Generate a longinteger number, which is characterize fragment with atom
		CHA in center. Up to second neigbour's sphere are took into account. Array
		A contains the atom's codes, generated by ALLATATOM, array B-bond's codes,
		generated by BONDCONVERSION*/

		int i,j,k,n,i1,i2,i3,jj;
		std::vector<int> aDist(listarSize());
		std::vector<int> aList(listarSize());
		std::vector<int> atomoBond(listarSize());
		std::vector<int> * aPrev [3];
		bool test,testOK;
		int w2;
		int l2,l3,tl1,tl2;
		std::vector<int> ab(listarSize());
		int result=0;
		//    String s;

		if (nAtoms()==0) return result;


		for (j=0; j<3; j++) aPrev[j]=new std::vector<int>(listarSize());
		for (i=0; i<nAtoms(); i++) {
			aDist[i]=65700;
			atomoBond[i]=0;
			for (j=0; j<3; j++) (*aPrev[j])[i]=-1;
		};
		for (i=0; i<sphere; i++) wSphere[i]=0;
		aDist[att]=0;
		k=0;
		atomoBond[att]=a[att];
		testOK=false;
		test=true;

		while (test) {
			test=false;
			k++;
			for (i1=0; i1<nAtoms(); i1++) if (aDist[i1]==(k-1)) { //K-sphere number, K-1-previous sphere
				if (bk[i1].nb > 0) for (j=0; j<bk[i1].nb; j++) {
					//for each neighbour
					i3=bk[i1].adjusted[j];//[I1].BN[J];
					i2=getBond(i3)->at[0];
					if (i2 == i1) i2=getBond(i3)->at[1];

					if ((aDist[i2] >= k) && (aDist[i2] <= 65700)) {
						//if the neighbour has not been added yet to list}
						aDist[i2]=k;          //Mark neighbour's sphere number
						jj=0;
						test=true;
						for (jj=1; jj<=3; jj++) {
							if (((*aPrev[jj-1])[i2] == -1) || (jj == 3)) {
								(*aPrev[jj-1])[i2]=i1;
								break;
							};
						};
						w2=a[i2];
						w2=w2 << 6;
						w2=w2+b[i3];
						w2=w2 ^ atomoBond[i2];  //atomoBond[I2]:=AtomoBond[I2] xor W2; !!!!!!!
						atomoBond[i2]=w2;
					};
				};
			};
			if (secAt < 0) testOK=(k==sphere); else testOK=((k>=sphere) && (aDist[secAt] < 65700));
			if (testOK) test=false;
		}; //end WHILE
		// At this point I have sphere numbers in Adist, corresponding AtomoBond word, number of atom(s), from
		//which the desirable atom were generated

		for (i=0; i<nAtoms(); i++) ab[i]=atomoBond[i];
		for (i1=k; i1>=1; i1--) {
			n=0;
			for (i=0; i<nAtoms(); i++) if (aDist[i] == i1) {
				aList[n]=i;
				n++;
			};
			if (n > 0) {
				i2=0;
				for (i=0; i<n; i++) {
					for (j=1; j<=3; j++) if ((*aPrev[j-1])[aList[i]] >= 0) {
						tl1=ab[aList[i]];
						tl2=atomoBond[(*aPrev[j-1])[aList[i]]]; // [Aprev[J][AList[I]]];
						tl2=tl2 << 9;
						tl1=tl1 ^ tl2;
						ab[aList[i]]=tl1;
					};
				};
				//Sorting
				if (n > 1) for (i=0; i<(n-1); i++) for (j=i+1; j<n; j++) if (ab[aList[i]] < ab[aList[j]]) {
					i2=aList[i]; aList[i]=aList[j]; aList[j]=i2;
				};
				if (n > 1) i2=15/(n-1);
				l2=ab[aList[0]];
				if (n > 1) for (i=1; i<n; i++) {
					l3=ab[aList[i]];
					l3=l3 << ((i+1)*i2);
					l2=l2 ^ l3;
				};
				l3=0;
				for (i=0; i<n; i++) l3=l3+getAtom(aList[i])->nc;
				if (l3 != 0) l2=-l2;
				for (i=0; i<n; i++) l2=l2 + a[aList[i]];  //Addition from 04.04.94}
				wSphere[i1-1]=l2;
			} else wSphere[i1-1]=0;
		};
		if ((secAt == -1) || (secAt == 65700)) result=0; else result=a[secAt];
		//freeing resources
		for (j=0; j<3; j++) delete(aPrev[j]);


		return result;
};


void TSimpleMolecule::makeEquivalentList(std::vector<int>& equivalenceList, bool isTopologyOnly) const {
	//AllAboutCycles must be called early

	std::vector<int> a(listarSize());
	std::vector<int> b(listarSize());
	TSimpleMolecule * em;
	int i,j,n;
	std::vector<adjustedlist> bk(listarSize());
	std::vector<std::vector<int> *> aeqList(0);
	std::vector<int> * lL;


	em= new TSimpleMolecule();
	em->moleculeCopy(*this);

	//Only topology-nothing more!!!!
	if (isTopologyOnly) {
		for (i=0; i<em->nAtoms(); i++) {
			em->getAtom(i)->na=6;
			em->getAtom(i)->nv=4;
			em->getAtom(i)->rl=0;
			em->getAtom(i)->nc=0;
		};
		for (i=0; i<em->nBonds(); i++) em->getBond(i)->tb=1;
		em->defineAtomConn();
	} else {
		for (i=0; i<em->nBonds(); i++) if ((em->getBond(i)->tb >= 9) && (em->getBond(i)->tb <= 11)) em->getBond(i)->tb=1;
	};

	em->defineBondConn(bk);
	//
	for (i=0; i<em->nAtoms(); i++) a[i]=em->getAtom(i)->allAtAtom();
	for (i=0; i<em->nBonds(); i++) b[i]=em->getBond(i)->bondConversion();
	for (i=0; i<em->nAtoms(); i++) {
		lL=new vector<int>(em->listarSize());
		for (j=0; j<lL->size(); j++) (*lL)[i]=0;
		em->fragmentSecond(em->nAtoms(),i,-1,a,b,bk,*lL);
		aeqList.push_back(lL);
	};
	equivalenceList.resize(em->nAtoms());


	n=0;
	for (i=0; i<equivalenceList.size(); i++) equivalenceList[i]=0;
	for (i=0; i<em->nAtoms(); i++) if (equivalenceList[i] == 0) {
		n++;
		equivalenceList[i]=n;
		if (i < (em->nAtoms()-1)) for (j=i+1; j<em->nAtoms(); j++) if (equivalenceList[j] == 0) {
			if (compareAtoms(j,i,aeqList)) equivalenceList[j]=n;
		}
	};

	//freeing resorces
	for (int i=0; i<aeqList.size(); i++) delete(aeqList.at(i));
	aeqList.clear();
	delete(em);
};

int TSimpleMolecule::correctOverlapped() {
	double r;
	TSimpleMolecule * smCopy=new TSimpleMolecule();
	TSimpleMolecule * bestStore=new TSimpleMolecule();//NULL;
	int i, j, n, k, kk;
	bool test;
	int result;
	std::vector<int> currentValues(nBonds()+nAtoms());
	std::vector<int> maxValues(nBonds() + nAtoms());
	std::vector<int> rotBondList(nBonds() + nAtoms());
	int at1,at2;
	int nn;


	r=averageBondLength()/blDenominator;
	k=hasOverlapped(r,false);
	result=k;
	smCopy->moleculeCopy(*this);
	test=true;
	while (test) {
		test=false;
		for (i=0; i<nBonds(); i++) if (getBond(i)->db == 0) {
			flipSmall(i);
			kk=hasOverlapped(r,false);
			if (kk<k) {
				test=true;
				k=kk;
				smCopy->moleculeCopy(*this);
			} else moleculeCopy(*smCopy);
		};
		if (k == 0) test=false;
	};
	result=k;
	if (result > 0)  {  //multiply bonds rotation...
		smCopy->moleculeCopy(*this);

		for (i=0; i<rotBondList.size(); i++) rotBondList[i]=0;
		for (i=0; i<nBonds(); i++) if (getBond(i)->db == 0) {
			at1=getBond(i)->at[0];
			at2=getBond(i)->at[1];
			if ((getAtom(at1)->nb >= 2) && (getAtom(at2)->nb >=2)) {
				//checking if greter 1
				test=false;
				if (getAtom(at1)->nb == 3) {
					test=true;
					for (j=0; j<getAtom(at1)->nb; j++) {
						k=getAtom(at1)->ac[j];
						if ((k != at2) && (getAtom(k)->nb > 1)) test=false;
					};
				};
				if ((! test) && (getAtom(at2)->nb == 3)) {
					test=true;
					for (j=0; j<getAtom(at2)->nb; j++) {
						k=getAtom(at2)->ac[j];
						if ((k != at1) && (getAtom(k)->nb > 1)) test=false;
					};
				};

				if (! test) rotBondList[i]=1;
			};
		};

		nn=0;
		for (i=0; i<rotBondList.size(); i++) if (rotBondList[i] == 1) nn++;
		if ((nn <= nRotBondsMax) && (nn>0)) {  //All possible combinations of rotations...
			currentValues.resize(nn);
			maxValues.resize(nn);
			for (i=0; i<nn; i++) {
				currentValues[i]=0;
				maxValues[i]=1;
			};
			bestStore->moleculeCopy(*smCopy);
			test=true;
			while (test && (k>0)) {
				moleculeCopy(*smCopy);
				n=0;
				for (i=0; i<rotBondList.size(); i++) if (rotBondList[i] != 0) {
					n++;
					if (currentValues[n-1] == 1) flipSmall(i);
				};
				kk=hasOverlapped(r,false);
				if (kk < k) {
					k=kk;
					bestStore->moleculeCopy(*this);
				};
				test=incrementValues(currentValues,maxValues);
			};// (k=0) or (not test);
			moleculeCopy(*bestStore);
			result=k;
		};
	};
	if (result>0) { //Bond enlarge....]
		k=hasOverlapped(r,false);
		n=1;
		while (n != 0) {
			n=0;
			for (i=0; i<nBonds(); i++) if (getBond(i)->db == 0) {
				smCopy->moleculeCopy(*this);
				smCopy->bondEnlarge(i);
				kk=smCopy->hasOverlapped(r,false);
				if (kk < k) {
					k=kk;
					n=i;
				};
				smCopy->flipSmall(i);
				kk=smCopy->hasOverlapped(r,false);
				if (kk < k) {
					k=kk;
					n=-i;
				};
			};
			if (abs(n) != 0) {
				bondEnlarge(abs(n));
				if (n < 0) flipSmall(abs(n));
			};
		};
		result=k;
	};
	//freeing resources
	if (smCopy != NULL) delete(smCopy);
	if (bestStore != NULL) delete(bestStore);

	return result;
};

std::string TSimpleMolecule::getAtomSymbol(int atAtom)
{
	if(atAtom >= 0 && atAtom < fAtom.size() && fAtom[atAtom]->na <= NELEMMCDL)
		return aSymb[fAtom[atAtom]->na];
	else
		return std::string("");
}

double TSimpleMolecule::getAtomMass(int atAtom) {
	if(atAtom >= 0 && atAtom < fAtom.size() && fAtom[atAtom]->na <= NELEMMCDL)
		return aMass[fAtom[atAtom]->na];
	else
		return 0.0;
};

void TSimpleMolecule::removeExplicitHydrogens(void)
{
	int				i, n;
	bool				test;
	//vector<int>		aList, NH, eraseFlag;
	vector<TSingleAtom *> tempAtom;
	vector<TSingleBond *> tempBond;
	TSingleBond * sb;

	test = false;
	for(i = 0; i < fAtom.size(); ++i)
		if(fAtom[i]->na == 1) {
			test = true;
			break;
		};
	
	if(!test) return;
	defineAtomConn();
	//determineFormula();
	//aList.clear();
	//NH.clear();
    //eraseFlag.clear();
	try {
	  while (test && (fAtom.size() > 1)) {
        test=false;
		tempAtom.clear();
		tempBond.clear();
		n=-1;
		for (i=fAtom.size()-1; i>=0; i--) if (getAtom(i)->na == 1) {
		  n=i;
          break;
		};
		if (n >= 0) {
          test=true;
		  for (i=0; i<nAtoms(); i++) if (i != n) tempAtom.push_back(fAtom[i]);
		  delete(fAtom[n]);
		  fAtom.resize(tempAtom.size());
		  for (i=0; i<tempAtom.size(); i++) fAtom[i]=tempAtom[i];
		  for (i=0; i< nBonds(); i++) {
			sb=getBond(i);
			if ((sb->at[0] == n) || (sb->at[1] == n)) delete sb; else {
			  if (sb->at[0] > n) sb->at[0]--;
			  if (sb->at[1] > n) sb->at[1]--;
			  tempBond.push_back(sb);
			};
		  };
		  fBond.resize(tempBond.size());
		  for (i=0; i<tempBond.size(); i++) fBond[i]=tempBond[i];
		};
	  };
	  defineAtomConn();

		/*
		aList.resize(fAtom.size(), 0);
		NH.resize(fAtom.size(), 0);
		eraseFlag.resize(fAtom.size(), 0);
		for(i = 0; i < fAtom.size(); ++i) {
   		  fAtom[i]->enumerator = i;
		  aList[i] = -1;
		  NH[i] = getNH(i);
		  test=((fAtom[i]->nb == 1) && (fAtom[fAtom[i]->ac[0]]->na != 1));	   
		  if (test) eraseFlag[i]=1;
		};
		for(i = fAtom.size(); i > 0; --i) if(eraseFlag[i-1] == 1) {
		  fAtom.erase(fAtom.begin() + (i-1));
		};
		for(i = 0; i < fAtom.size(); ++i) {
			n = fAtom[i]->enumerator;
			aList[n] = i;
		};

		for(i = fBond.size(); i > 0; --i) {
			at1 = fBond[i-1]->at[0];
			at2 = fBond[i-1]->at[1];
			if(aList[at1] == -1 || aList[at2] == -1) {
				fBond.erase(fBond.begin() + (i-1));
			} else {
				fBond[i-1]->at[0] = aList[at1];
				fBond[i-1]->at[1] = aList[at2];
			};
		};
		defineAtomConn();

    //determineFormula();
	 
		for(i = 0; i < NH.size(); ++i) { //valencies corrections
			n = aList[i];
			if(n >= 0) {
				k = getNH(n);
				if(k < NH[i]) 
					fAtom[n]->nv = fAtom[n]->nv + (NH[i] - k);
			}
		}
		*/
	} catch(...) {
	};
}





void TSimpleMolecule::cb_removeExplicitHydrogens(void)
{
	int				i, n;
	bool				test;
	vector<TSingleAtom *> tempAtom;
	vector<TSingleBond *> tempBond;
	TSingleBond * sb;

	test = false; for(i = 0; i < fAtom.size(); ++i)if(fAtom[i]->na == 1) {test = true;break;}; if(!test) return;
	try {
	  while (test && (fAtom.size() > 1)) {
        test=false;
		tempAtom.clear();
		tempBond.clear();
		n=-1;
		for (i=fAtom.size()-1; i>=0; i--) if (getAtom(i)->na == 1) {
		  n=i;
          break;
		};
		if (n >= 0) {
          test=true;
		  for (i=0; i<nAtoms(); i++) if (i != n) tempAtom.push_back(fAtom[i]);
		  delete(fAtom[n]);
		  fAtom.resize(tempAtom.size());
		  for (i=0; i<tempAtom.size(); i++) fAtom[i]=tempAtom[i];
		  for (i=0; i< nBonds(); i++) {
			sb=getBond(i);
			if ((sb->at[0] == n) || (sb->at[1] == n)) delete sb; else {
			  if (sb->at[0] > n) sb->at[0]--;
			  if (sb->at[1] > n) sb->at[1]--;
			  tempBond.push_back(sb);
			};
		  };
		  fBond.resize(tempBond.size());
		  for (i=0; i<tempBond.size(); i++) fBond[i]=tempBond[i];
		};
	  };
	  defineAtomConn();
	} catch(...) {};
}






void TSimpleMolecule::rescaleToLength(double newAvgBondLength){
  double oldBondLength=this->averageBondLength();
  double xMin, yMin;

  if (oldBondLength <= 0.000001) return;

  xMin=getAtom(0)->rx;
  yMin=getAtom(0)->ry;
  for (int i=0; i<nAtoms(); i++) {
	getAtom(i)->rx=newAvgBondLength*getAtom(i)->rx/oldBondLength;
	getAtom(i)->ry=newAvgBondLength*getAtom(i)->ry/oldBondLength;
	if (getAtom(i)->rx < xMin) xMin=getAtom(i)->rx;
	if (getAtom(i)->ry < yMin) yMin=getAtom(i)->ry;
  };
  for (int i=0; i<nAtoms(); i++) {
	getAtom(i)->rx=getAtom(i)->rx-xMin;
	getAtom(i)->ry=getAtom(i)->ry-yMin;
  };

};

TSimpleMolecule& TSimpleMolecule::operator=(const TSimpleMolecule& sm)
{
	if(this != &sm) {
		TSimpleMolecule&		smn = const_cast<TSimpleMolecule&>(sm);
		moleculeCopy(smn);
	};
	return *this;
}


TSimpleMolecule::TSimpleMolecule(const TSimpleMolecule& sm)
{
	(*this) = sm;
}

//Ring calculation routines

/*
bool listPresent(std::vector<std::vector<int> *> ringData, std::vector<int> tempList) {
  bool result=false;
  std::vector<int> * lr;
  int  i,j;
  bool test;
  
  if (ringData.size() == 0) return result;
  
  for (i=0; i<ringData.size(); i++) {
    lr=ringData.at(i);
	
	if ((lr->size() == tempList.size()) && (lr->size() > 0)) {
      test=true;
	  for (j=0; j<lr->size(); j++) if ((*lr)[j] != tempList[j]) {
        test=false;
        break;
      }
      result=test;
    }
    if (result) break;
  };
  return result;
}


bool bondPresent(std::vector<int> list, int bn) {
  int i;
  bool result=false;
  for (i=0; i< list.size(); i++) if (list[i] == bn) {
    result=true;
    break;
  };
  return result;
};
*/

int TSimpleMolecule::cyclesStatistic(std::vector<int>& nCyclesSize, std::vector<int>& cycleDetails, bool testEnumerator) {
//It is assumed, that AllAboutCycles routine has been called early}
//nCycleSize[0] contains number of aromatic 5-membered cycles
//nCycleSize[1] contains number of aromatic 6-membered cycles
//nCycleSize[2] contains number of cycles, greater high bound of nCycleSize
//nCycleSize[3]..nCycleSize[n] - number of cycles of given size.

  int result=0;
  std::vector<std::vector<int> *> ringData;
  std::vector<int> bondList,ar,bondTypes,atomNoList;
  std::vector<int> * lRef;
  std::vector<int> * tempList;
  std::vector<adjustedlist> bk(listarSize());
  int rs;
  int i,j,k,n,m;
  std::vector<int> bondCount(listarSize());
  bool testUniq;
  int at1,at2;
  bool test;

  
  defineBondConn(bk);
  bondList.resize(listarSize());
  ar.resize(listarSize());

 //addition
  bondTypes.resize(listarSize());
  atomNoList.resize(listarSize());   //+1 -because of Delphi remove in C++ code

   
  if (nCyclesSize.size() < 3) return result;
  for (i=0; i<nCyclesSize.size(); i++) nCyclesSize[i]=0;
 
  if (fBond.size() > 0) for (i=0; i<fBond.size(); i++) {
      //Addition
	n=getBond(i)->db;
    if ((n==2) || (n==3)) n=4; else n=getBond(i)->tb;
    bondTypes[i]=n;
      //end addition

    vaweBond(i,bk,rs,bondList);
	if (rs > 0) {
	  tempList=new vector<int>(0);
	  for (j=0; j<rs; j++) {
        m=bondList[j];
		tempList->push_back(m);
	  };
	  for (j=0; j<(rs-1); j++) for (k=j+1; k<rs; k++) if ((*tempList)[j] > (*tempList)[k]) {
		m=(*tempList)[j];
		(*tempList)[j]=(*tempList)[k];
		(*tempList)[k]=m;
	  };
     // tempList.Sort(sortList);
	  if (! listPresent(ringData,*tempList)) ringData.push_back(tempList);
	};
  };

//    {06042012}
    //removing with non-zero enumerators....
  if (testEnumerator) for (i=ringData.size()-1; i>=0; i--) {
    lRef=(std::vector<int> *)ringData.at(i);
    test=false;
	for (j=0; j<lRef->size(); j++) {
      n=(*lRef)[j];
	  at1=getBond(n)->at[0];
	  at2=getBond(n)->at[1];
	  if ((getAtom(at1)->enumerator != 0) || (getAtom(at2)->enumerator != 0)) {
        test=true;
        break;
	  };
	};
	if (test) {
  	  ringData[i]=NULL;
      delete(lRef);
    };
  };
//    {/06042012}


  //removing depended rings
  /*
  for i:=0 to NBONDSMAX-1 do bondCount[i]:=0;
  for i:=0 to ringData.Count-1 do begin
    lRef:=TList(RingData[i]);
    for j:=0 to lRef.Count-1 do begin
      m:=integer(lRef[j]);
      bondCount[m]:=bondCount[m]+1;
    end;
  end;
  for i:=ringData.count-1 downto 0 do begin
    lRef:=TList(RingData[i]);
    testUniq:=true;
    for j:=0 to lRef.Count-1 do begin
      m:=integer(lRef[j]);
      if bondCount[m]=1 then begin
        testUniq:=false;
        break;
      end;
    end;
    if testUniq then begin
      for j:=0 to lRef.Count-1 do begin
        m:=integer(lRef[j]);
        bondCount[m]:=bondCount[m]-1;
      end;
      freeAndNil(lRef);
      ringData.Delete(i);
    end;
  end;
  */

  //finalyze result
  result=ringData.size();
  if (ringData.size() > 0) for (i=0; i<ringData.size(); i++) if (ringData[i] != NULL) {
	lRef=(std::vector<int> *)ringData.at(i);
	n=lRef->size();
	if (n < nCyclesSize.size()) nCyclesSize[n]=nCyclesSize[n]+1; else nCyclesSize[2]=nCyclesSize[2]+1;
	if (aromatic(n,*lRef,ar)) {
        //    if (aromatic(cycleSize,bondList,ar)) {       //store, if aromatic cycle}
      if (n == 5) nCyclesSize[0]=nCyclesSize[0]+1; else nCyclesSize[1]=nCyclesSize[1]+1;
      for (j=0; j<n; j++) ar[(*lRef)[j]]=1;
    };

        //size
     cycleDetails[3*i]=n;
     n=0;
        //number of in double (or aromatic) bonds
	 for (j=0; j<lRef->size(); j++) {
       k=(*lRef)[j];
       k=bondTypes[k];
       if ((k==2) || (k==4)) n++;
	 };
     cycleDetails[3*i+1]=n;
     //Number of attached double bonds (methylene)
     for (j=0; j<atomNoList.size(); j++) atomNoList[j]=0;
	 for (j=0; j<lRef->size(); j++) {
       k=(*lRef)[j];
	   at1=getBond(k)->at[0]; 
	   at2=getBond(k)->at[1];
       atomNoList[at1]=1;
       atomNoList[at2]=1;
     };
     n=0;
	 for (j=0; j<fBond.size(); j++) {
       at1=getBond(j)->at[0];
	   at2=getBond(j)->at[1];
       if (((atomNoList[at1]==0) && (atomNoList[at2]!=0)) ||
		   ((atomNoList[at2]==0) && (atomNoList[at1]!=0))) {
         k=bondTypes[j];
         if ((k==2) || (k==4)) n++;
	   };
	 };
     cycleDetails[3*i+2]=n;
  };
  //cleanup resources
  if (ringData.size() > 0) for (i=0; i<ringData.size(); i++) if (ringData[i] != NULL) {
	lRef=(std::vector<int> *)ringData.at(i);
	delete(lRef);
  }
  return result;
};

int TSimpleMolecule::getNGauscheCarbons(){

  int i,j,n,n1,n2,nC1,nC2;
  bool test;
  int result=0;

  for (i=0; i<nBonds(); i++) if ((getBond(i)->tb == 1) && (getBond(i)->db == 0)) {  //Single acyclic bond

	n1=getBond(i)->at[0]; n2=getBond(i)->at[1];
    nC1=0; nC2=0;  //Calculation of attached carbons...
	for (j=0; j<getAtom(n1)->nb; j++) {
      n=getAtom(n1)->ac[j];
      if ((n != n2) && (getAtom(n)->na == 6)) nC1++;    //Only carbons-other atoms might have other Gaushe-values (halogens)
	};
	
	for (j=0; j<getAtom(n2)->nb; j++) {
      n=getAtom(n2)->ac[j];
      if ((n != n1) && (getAtom(n)->na == 6)) nC2++;
	};
    if ((nC1 >= 1) && (nC2 >= 1) && ((nC1+nC2) >= 3)) {
      //checking for the attached double-bond-alkenes have no gosche-correction
      test=true;
	  for (j=0; j<nBonds(); j++) if ((getBond(j)->tb == 2)
	  && ((getBond(j)->at[0] == n1) || (getBond(j)->at[0] == n2) || (getBond(j)->at[1] == n1) || (getBond(j)->at[1] == n2))) {
        test=false;
        break;
	  };
	  if (test) {
        if (nC1 == 2) result=result+1; else if (nC1 == 3) result=result+3;
        if (nC2 == 2) result=result+1; else if (nC2 == 3) result=result+3;
	  };
    };
  };
  return result;
};

int TSimpleMolecule::getNOrtho() {
  int i,j,n1,n2;
  bool test1,test2;
  int result=0;

  for (i=0; i<nBonds(); i++)  if ((getBond(i)->db == 2) || (getBond(i)->db == 3)) {  //Single acyclic bond
    n1=getBond(i)->at[0]; n2=getBond(i)->at[1];
    test1=false; test2=false;
	for (j=0; j<nBonds(); j++) if ((getBond(j)->tb == 1) && (getBond(j)->db == 0)) {
      if ((getBond(j)->at[0] == n1) || (getBond(j)->at[1]== n1)) test1=true;
      if ((getBond(j)->at[0] == n2) || (getBond(j)->at[1]== n2)) test2=true;
      if (test1 && test2) break;
	};
    if (test1 && test2) result++;
  };
  return result;
};

void TSimpleMolecule::atomBondChange() {
	//substitutes the semipolar bond with double. NBONDS-total number of bonds,
	//BOND-bond's attributes, ATOM-atom's attributes}
	int i;
	int ca,cb;

	if (nBonds()==0) return;
	for (i=0; i<nBonds(); i++) {
		ca=getAtom(getBond(i)->at[0])->nc;
		cb=getAtom(getBond(i)->at[1])->nc;
		if ((((ca<0) && (cb>0)) || ((ca>0) && (cb<0))) && ((getBond(i)->tb<3) || (getBond(i)->tb>8))) {
			if (ca<0) {
				getAtom(getBond(i)->at[0])->nc=ca+1;
				getAtom(getBond(i)->at[1])->nc=cb-1;
			};
			if (ca>0) {
				getAtom(getBond(i)->at[0])->nc=ca-1;
				getAtom(getBond(i)->at[1])->nc=cb+1;
			};
			if (getBond(i)->tb<3) getBond(i)->tb=getBond(i)->tb+1; else getBond(i)->tb=2;
		};
	};
};

bool TSimpleMolecule::stereoBondChange() {
	//Subtitute stereo bonds with single. On output TEST has TRUE value, if at
	//least one stereo bond was detected in structure with bond attributes BOND
	//nd total number of bonds NBONDS
	int i;
	bool result;

	result=false;
	if ((nBonds()==0) || (fIOPT12==2)) return result;
	//If Exact Stereo Search-no substitution
	for (i=0; i<nBonds(); i++)  {
		if (fIOPT12==3) {
			//If Replace Stereo Search-only EITHER bond substitution}
			if (getBond(i)->tb==11) getBond(i)->tb=1;
			if (getBond(i)->tb>=9) result=true;
		};
		if (fIOPT12==1) if (getBond(i)->tb>=9) getBond(i)->tb=1;
		//If No Stereo search-all bonds substitution
	};
	return result;
};


void TSimpleMolecule::correctValencies() {
  int i,n,k;

  for (i=0; i<nAtoms(); i++) {
	if ((getAtom(i)->nv > 0) && (getAtom(i)->nv < hVal[getAtom(i)->na])) {
	  n=getAtom(i)->nv;
      k=getNH(i);
      getAtom(i)->nv=hVal[getAtom(i)->na];
      if (getNH(i) != k) getAtom(i)->nv=n;
	};
  };
};

bool TSimpleMolecule::deleteHydrogenDecCharge(int an, int anInc) {
  int nH,i,k;
  bool result=false;

  nH=getAtom(an)->nv;
  nH=nH-getAtom(an)->currvalence-abs(getAtom(an)->nc)-getAtom(an)->rl;
  if (nH <= 0) {
    for (i=0; i<getAtom(an)->nb; i++) {
      k=getAtom(an)->ac[i];
      if (getAtom(k)->na == 1) {
        getAtom(anInc)->nc=getAtom(anInc)->nc+1;
        getAtom(an)->nc=getAtom(an)->nc-1;
        deleteAtom(k);
        result=true;
        break;
	  };
	};
  } else if ((getAtom(an)->na == 7) && (getAtom(an)->nv == 5)) {
    getAtom(an)->nv=3;
    getAtom(anInc)->nc=getAtom(anInc)->nc+1;
    getAtom(an)->nc=getAtom(an)->nc-1;
    result=true;
  };
  return result;
};


void TSimpleMolecule::makeStandardStructure() {
//Not necessary call AllAboutCycles. returns total charge (important for ionic compound savings!!!). Performs chemical structure normalization
  int i,j,k,m,n,nH,nP;
  bool testRemove;

  fIOPT12=1;
  fIOPT13=true;
  atomBondChange(); //Semipolar to double
  defineAtomConn();
  correctValencies();
  //Five-valent nitrogen normalization and Five-valent boron normalization
  testRemove=true;
  while (testRemove) {
    testRemove=false;
	for (i=0; i<nAtoms(); i++) if ((getAtom(i)->na == 7) || (getAtom(i)->na == 5)) {
      //Five-valent
      if (getAtom(i)->nv == 5) {
        getAtom(i)->nv=3;
		n=getAtom(i)->currvalence+getAtom(i)->nc;
        if (n < 5) for (j=n+1; j<=5; j++) addDefaultHydrogen(i,false);
        testRemove=true;
        defineAtomConn();
	  };
      //Five-coordinational
      if (! testRemove) if (((getAtom(i)->nb >= 4) && (getAtom(i)->nc == 0)) || (getAtom(i)->nb == 5)) {
        nH=0;
        nP=0;
        for (j=0; j<getAtom(i)->nb; j++) {
          k=getAtom(i)->ac[j];
          if (getAtom(k)->na == 1) nH++;
          if ((getAtom(k)->na == 8) || isHalogen(getAtom(k)->na)) nP++;
		};
        if ((nP== 1) && (nH > 0)) {
          //Search for polar...
          for (j=0; j<getAtom(i)->nb; j++) {
            k=getAtom(i)->ac[j];
            if ((getAtom(k)->na == 8) || isHalogen(getAtom(k)->na)) {
              //Search for bond...
              for (m=0; m<nBonds(); m++) if (((getBond(m)->at[0] == i) && (getBond(m)->at[1] == k))
              || ((getBond(m)->at[0] == k) && (getBond(m)->at[1] == i))) {
                deleteBond(m);
                testRemove=true;
                break;
			  };
		    };
            if (testRemove) break;
		  };
          //Search for hydrogen
          for (j=0; j<getAtom(i)->nb; j++) {
            k=getAtom(i)->ac[j];
            if (getAtom(k)->na == 1) {
              deleteAtom(k);
              break;
			};
		  };
          testRemove=true;
          defineAtomConn();
	    };
	  };
      if (testRemove) break;
    };
  };  //until not TestRemove;
  //Positive carbon processing   - wait for better time...
  /*
  for I:=1 to MainMol.FAtom.NAtoms do if (MainMol.FAtom[I].NA=6) and (MainMol.FAtom[I].NC=1) then begin
    NH:=0;
    //Nitrogen
    for J:=1 to MainMol.FAtom[I].NB do begin
      K:=MainMol.FAtom[I].AC[J];
      if
    end;
  end;
  */
  //Salts
  n=0;
  for (i=0; i<nAtoms(); i++) n=n+getAtom(i)->nc;
  if (n == 0) {  //The normalization can be perfomed in single way only for balanced charges
	testRemove=true;
    while (testRemove) { //Ammonium
      testRemove=false;
      //Appropriate proton was found...
      for (i=0; i<nAtoms(); i++) if ((getAtom(i)->nc > 0) && (getNH(i) > 0)) {
        //Halogen
        for (j=0; j<nAtoms(); j++) if ((getAtom(j)->nc < 0) && isHalogen(getAtom(j)->na)) {
          if (deleteHydrogenDecCharge(i,j)) {
            testRemove=true;
            defineAtomConn();
		  };
          if (testRemove) break;
		};
       //Oxygen, non semipolar
        if (! testRemove) {
          for (j=0; j<nAtoms(); j++) if ((getAtom(j)->nc < 0) && (getAtom(j)->na == 8)
          && (getAtom(j)->nb == 1) && (getAtom(getAtom(j)->ac[0])->nc == 0)) {
            if (deleteHydrogenDecCharge(i,j)) {
              testRemove=true;
              defineAtomConn();
			};
            if (testRemove) break;
		  };
        };
        //Oxygen, semipolar
        if (! testRemove) {
          for (j=0; j<nAtoms(); j++) if ((getAtom(j)->nc < 0) && (getAtom(j)->na == 8)
          && (getAtom(j)->nb == 1) && (getAtom(getAtom(j)->ac[0])->nc > 0)) {
            n=0;
            m=getAtom(j)->ac[0];
            for (k=0; k<getAtom(m)->nb; k++) n=n+getAtom(getAtom(m)->ac[k])->nc;
            if (abs(n) > 1) if (deleteHydrogenDecCharge(i,j)) {
              testRemove=true;
              defineAtomConn();
			};
            if (testRemove) break;
		  };
		};
        //Any other anions (sulfides). Non-semipolar
        if (!testRemove) {
          for (j=0; j<nAtoms(); j++) if (getAtom(j)->nc < 0) {
            n=0;
            if (getAtom(j)->nb > 0) for (0; k<getAtom(j)->nb; k++) n=n+getAtom(getAtom(j)->ac[k])->nc;
            if (n == 0) if (deleteHydrogenDecCharge(i,j)) {
              testRemove=true;
              defineAtomConn();
			};
            if (testRemove) break;
		  };
		};
        if (testRemove) break;
	  };
	};//until not TestRemove;
  };
  //Metal normalization - ionic structures where possible
  testRemove=true;
  while (testRemove) {
    testRemove=false;
    for (i=0; i<nAtoms(); i++) if ((getAtom(i)->nb > 0) && isMetall(getAtom(i)->na)) {
      for (m=0; m<getAtom(i)->nb; m++) {
        k=getAtom(i)->ac[m];
        for (j=0; j<nBonds(); j++) if (getBond(j)->tb == 1) if (((getBond(j)->at[0] == i) && (getBond(j)->at[1] == k)) 
			|| ((getBond(j)->at[0] == k) && (getBond(j)->at[1] == i))) {
          deleteBond(j);
          getAtom(i)->nc=getAtom(i)->nc+1;
          getAtom(k)->nc=getAtom(k)->nv-1;
          testRemove=true;
          defineAtomConn();
          break;
		};
        if (testRemove) break;
	  };
      if (testRemove) break;
	};
  }; //until not TestRemove;
};

void TSimpleMolecule::correctAzulenes(const std::vector<adjustedlist> & bk) {
  std::vector<int> bL7(listarSize());
  std::vector<int> bL5(listarSize());
  int i,j,n,k;
  int w7,w5;
  bool test;

  //By unknown reasom 7-membered aromatic is bond type=8 is used
  for (i=0; i<nBonds(); i++) if (getBond(i)->tb == 8) getBond(i)->tb=4;

  for (i=0; i<nBonds(); i++) if (getBond(i)->db == 8) {
	vaweBond(i,bk,w7,bL7); //If I-th bond belongs to cycle
    n=0;
    k=0;
    for (j=0; j<w7; j++) {
      if (getBond(bL7[j])->tb == 2) n++;
      if (getBond(bL7[j])->tb == 4) k++;
	};
    test=false;
    if ((n == 3) || (k == 7)) for (j=0; j<w7; j++) if (getBond(bL7[j])->db == 6) {
      vaweBond(bL7[j],bk,w5,bL5); 
      n=0;
      for (k=0; k<w5; k++) if (getBond(bL5[k])->tb == 2) n++;
      if (n == 2) {
		for (k=0; k<w7; k++) getBond(bL7[k])->db=3;
		for (k=0; k<w5; k++) getBond(bL5[k])->db=2;
        test=true;
	  };
      if (test) break;
	};
  };
};


void TSimpleMolecule::tautomerAnalize() {
//AllAboutCycles MUST be called prior the procedure calling
//!!! Enumerator will be destroyed in the TSimpleMolecule !!!!
  std::vector<adjustedlist> bk(listarSize());
  bool testAdded,test,test1;
  int i,j,k,nL,m,n,s,d,e,m1,m2;
  std::vector<int> chargeStore(nAtoms());
  std::vector<int> nH(nAtoms());
  bool arom1,arom2;
  int an1,an2;
  int k1,k2,k3;
  bool testArom1,testArom2;
  std::vector<int> possibleDouble(nBonds() + nAtoms());
  std::vector<int> hasPIElectrons(nAtoms()); 
  std::vector<int> hasAromatic(nAtoms());
  std::vector<int> possibleH(nAtoms());
  std::vector<int> atomDefine;
  //bool firstRun;
  int w;
//  FF:TextFile;

  //Addition from 11 September 2006 - negative charges have to be processed here
  for (i=0; i<nAtoms(); i++) {
	nH[i]=getNH(i)+getAtom(i)->nc;
	chargeStore[i]=getAtom(i)->nc;
	getAtom(i)->nc=0;
  };

  defineBondConn(bk);
  correctAzulenes(bk);
  //firstRun:=true;
  for (int runCounter=0; runCounter<2; runCounter++) {
    //20:
	for (i=0; i<possibleH.size(); i++) possibleH[i]=0;
    for (i=0; i<nBonds(); i++) {
	  getBond(i)->enumerator=0;
      if ((getBond(i)->tb == 2) || (getBond(i)->tb == 3) || (getBond(i)->db == 2) || (getBond(i)->db == 3) || (getBond(i)->tb == 4)) possibleDouble[i]=1; else possibleDouble[i]=0;
	};
    //possible double processing
    e=0;
	testAdded=true;
    while (testAdded) {
      e++;
      testAdded=false;
      for (i=0; i<nBonds(); i++) if (possibleDouble[i] == e) {
        k1=getBond(i)->at[0];
        k2=getBond(i)->at[1];
		for (j=0; j<bk[k2].nb; j++) {
		  k=bk[k2].adjusted[j];  //bond number attached
          if (k != i) {
            n=getBond(k)->at[0];
            if (n == k2) n=getBond(k)->at[1];  //K2- central atom, K1 - atom at end double bond N - atom at end single bond
            if (isHetero(getAtom(k1)->na)) {
              test=nH[n] > 0;
			} else {
              test=(nH[n] > 0) && isHetero(getAtom(n)->na);
			};
            if (test) {
              possibleH[k1]=1;
              getBond(i)->enumerator=1;
              getBond(k)->enumerator=1;
              if (possibleDouble[k] == 0) {
                possibleDouble[k]=e+1;
                testAdded=true;
			  };
			};
		  };
		};
        k1=getBond(i)->at[1];
        k2=getBond(i)->at[0];
		for (j=0; j<bk[k2].nb; j++) {
		  k=bk[k2].adjusted[j];  //bond number attached
          if (k != i) {
            n=getBond(k)->at[0];
            if (n == k2) n=getBond(k)->at[1];  //K2- central atom, K1 - atom at end double bond N - atom at end single bond
            if (isHetero(getAtom(k1)->na)) {
              test=nH[n] > 0;
			} else {
              test=(nH[n] > 0) && isHetero(getAtom(n)->na);
			};
            if (test) {
              possibleH[k1]=1;
              getBond(i)->enumerator=1;
              getBond(k)->enumerator=1;
              if (possibleDouble[k] == 0) {
                possibleDouble[k]=e+1;
                testAdded=true;
			  };
			};
		  };
		};
	  };
	};//until not testAdded;
    //End possible double processing-process again because of some H atoms might arisen...
	

    //possible H processing
    for (i=0; i<nBonds(); i++) if (possibleDouble[i] > 0) {
      k1=getBond(i)->at[0];
      k2=getBond(i)->at[1];
	  for (j=0; j<bk[k2].nb; j++) {
		k=bk[k2].adjusted[j];  //bond number attached
        if (k != i) {
          n=getBond(k)->at[0];
          if (n == k2) n=getBond(k)->at[1];  //K2- central atom, K1 - atom at end double bond N - atom at end single bond
          if (isHetero(getAtom(k1)->na)) {
            test= possibleH[n] > 0;
		  } else {
            test=(possibleH[n] > 0) && isHetero(getAtom(n)->na);
		  };
          if (test) {
			getBond(i)->enumerator=1;
            getBond(k)->enumerator=1;
		  };
		};
	  };
      k1=getBond(i)->at[1];
      k2=getBond(i)->at[0];
	  for (j=0; j<bk[k2].nb; j++) {
		k=bk[k2].adjusted[j];  //bond number attached
        if (k != i) {
          n=getBond(k)->at[0];
          if (n == k2) n=getBond(k)->at[1];  //K2- central atom, K1 - atom at end double bond N - atom at end single bond
          if (isHetero(getAtom(k1)->na)) {
            test= possibleH[n] > 0;
		  } else {
            test=(possibleH[n] > 0) && isHetero(getAtom(n)->na);
		  };
          if (test) {
			getBond(i)->enumerator=1;
            getBond(k)->enumerator=1;
		  };
		};
	  };
	};
    //end possibleH processing
	for (i=0; i<nBonds(); i++) if ((getBond(i)->enumerator != 0) || (getBond(i)->db == 2) || (getBond(i)->db == 3)) getBond(i)->tb=4;

    //More extend-determine number of Pi-electrons and move them to conjugated system
	for (i=0; i<hasPIElectrons.size(); i++) hasPIElectrons[i]=0;
	for (i=0; i<hasAromatic.size(); i++) hasAromatic[i]=0;
    for (i=0; i<nAtoms(); i++) if (isHetero(getAtom(i)->na)) hasPIElectrons[i]=1;
    for (i=0; i<nBonds(); i++) if ((getBond(i)->tb == 2) || (getBond(i)->tb == 3) || (getBond(i)->tb == 4) || (getBond(i)->db == 2) || (getBond(i)->db == 3)) {
      k=getBond(i)->at[0];
      hasPIElectrons[k]=1;
      if ((getBond(i)->tb == 4) || (getBond(i)->db == 2) || (getBond(i)->db == 3)) hasAromatic[k]=1;
      k=getBond(i)->at[1];
      hasPIElectrons[k]=1;
	  if ((getBond(i)->tb == 4) || (getBond(i)->db == 2) || (getBond(i)->db == 3)) hasAromatic[k]=1;
	};
	test=true;
	while (test) {
      test=false;
      for (i=0; i<nBonds(); i++) if (getBond(i)->tb != 4) {
        k1=getBond(i)->at[0];
        k2=getBond(i)->at[1];
        if ((hasPIElectrons[k1] > 0) && (hasPIElectrons[k2] > 0)) {
          if ((hasAromatic[k1] > 0) || (hasAromatic[k2] > 0)) {
            getBond(i)->tb=4;
            hasAromatic[k1]=1;
            hasAromatic[k2]=1;
            test=true;
	      };
	    };
	  };
	};//until not test;

    //Non-aromatic bond connects pair of atoms
	for (i=0; i<nBonds(); i++) getBond(i)->enumerator=0;
    for (i=0; i<nBonds(); i++) if (getBond(i)->tb != 4) {
      arom1=false;
      arom2=false;
      for (j=0; j<nBonds(); j++) if (getBond(j)->tb == 4) {
        if (getBond(j)->at[0] == getBond(i)->at[0]) arom1=true;
        if (getBond(j)->at[0] == getBond(i)->at[1]) arom2=true;
        if (getBond(j)->at[1] == getBond(i)->at[0]) arom1=true;
        if (getBond(j)->at[1] == getBond(i)->at[1]) arom2=true;
        if (arom1 && arom2) break;
	  };
      if (arom1 && arom2) getBond(i)->enumerator=1;
	};
	for (i=0; i<nBonds(); i++) if (getBond(i)->enumerator != 0) getBond(i)->tb=4;

    //Checking in new aromatic rings formed and execute again if it is.
  };
  //restore charges
  for (i=0; i<chargeStore.size(); i++) getAtom(i)->nc=chargeStore[i];

    //Search for special groups - nitro and nitrilo, connected to aromatic
  for (i=0; i<nAtoms(); i++) {
    n=0;
    if ((getAtom(i)->na == 7) && (getAtom(i)->nb == 3)) {
      n=1;
	} else if ((getAtom(i)->na == 6) && (getAtom(i)->nb == 2)) {
      n=2;
	};
	
    if (n > 0) {
      k1=0; k2=0; k3=0;
	  for (j=0; j<bk[i].nb; j++) {
        k=bk[i].adjusted[j];
        if ((getBond(k)->tb == 1) || (getBond(k)->tb == 4)) {
          if (k1 == 0) k1=k; else k1=-1;
		} else if (getBond(k)->tb == 2) k2=k2+1; else
        if (getBond(k)->tb == 3) k3=k3+1;
	  };
      if (k1 > 0) {
        test=false;
        if ((n == 1) && (k2 == 2)) {
          test=true;
		} else if ((n == 2) && (k3 == 1)) {
          test=true;
		};
        if (test) {
          n=getBond(k1)->at[0];
          if (n == i) n=getBond(k1)->at[1];
          test=false;
		  for (j=0; j<bk[n].nb; j++) {
			k=bk[n].adjusted[j];
            if (getBond(k)->tb == 4) test=true;
		  };
		};
        if (test) for (j=0; j<bk[i].nb; j++) {
		  k=bk[i].adjusted[j];
          getBond(k)->tb=4;
		};
	  };
	};
  };
};

void TSimpleMolecule::removeHydrogen(std::vector<int> * qHydr, std::vector<int> * qEnumerator) {

	//Explicitly-defined hydrogens are removed from the structure, which is characte-
	//rized of total number of atoms NA, total number of bonds NB, atom's attribute
	//array ATOM, bond's attribute array BOND, bond-connection matrix invariants
	//ONN. Number of explicitly-defined hydrogens for each atom are stored in
	//QHYDR. If QUERYLABEL points on some atom, the value is recalculated so that
	//the same atom is pointed after hydrogens have been removed from the structure

	//NB! QueryLabel should be array!}
	int i,j,i1,n,mm,kk;
	bool test,test1,test2;

	if (qEnumerator != NULL) {
		if (qEnumerator->size() != nAtoms()) qEnumerator->resize(nAtoms());
		for (i=0; i<qEnumerator->size(); i++) (*qEnumerator)[i]=i;
		//for (i=0; i<nAtoms(); i++) getAtom(i)->astereo=i;
	};
	for (i=0; i<nAtoms(); i++) {
		if (qHydr != NULL) (*qHydr)[i]=0;
		if ((! fIOPT11) && (getAtom(i)->na==104)) getAtom(i)->na=1;
		//D->H if no isotop sensitivity}
	};
	i=0;
	test1=false;
	if (nAtoms()>0) while (i<nAtoms()) {
	    mm=0;
		if ((getAtom(i)->na == 1) && (getAtom(i)->nb == 1) && (getAtom(getAtom(i)->ac[0])->na != 1))  { //hydrogen found
		  for (kk=0; kk<nBonds(); kk++) {
			if (getBond(kk)->at[0] == i) mm++;
			if (getBond(kk)->at[1] == i) mm++;
		  };
		};
        if (mm == 1) {
			test1=true;
			j=0;
			i1=-1; 
			test2=true;

			if (nBonds()>0) 
				do { //search for corresponding bond and hydrogen's neighbour detection
					test=(getBond(j)->at[0]==i) || (getBond(j)->at[1]==i);
					if (test) {
						i1=getBond(j)->at[0];
						if (i1==i) i1=getBond(j)->at[1];
					};
					j++;
					test2 = test || (j == nBonds());
				} while (!test && j < nBonds());
				//} while (!test2);

			deleteAtom(i);
			if ((qHydr != NULL) && (i<(nAtoms()-1))) deleteIntElement(qHydr,i);

			if (i1>i)  i1=i1-1;
			//shift of attribute's arrays
			
			if (qEnumerator != NULL) {
				for (j=0; j<qEnumerator->size(); j++) {
					if ((*qEnumerator)[j] == i) (*qEnumerator)[j]=-i1; else
						if ((*qEnumerator)[j] > i)  (*qEnumerator)[j]=(*qEnumerator)[j]-1; else
							if (((*qEnumerator)[j] < 0) && (abs((*qEnumerator)[j]) > i)) (*qEnumerator)[j]=(*qEnumerator)[j]+1;
				};
			};
			
			
			if ((i1 >= 0) && (qHydr != NULL)) (*qHydr)[i1]=(*qHydr)[i1]+1;
			//counter of explicitly defined hydrogens
			i--;
defineAtomConn();
		};
		i++;
	};
	if (test1) {
		defineAtomConn();
		//Inverse enumerator creation... New QA->OldQA array}
		n = 0;
		if (qEnumerator != NULL)  for (i = 0; i < nAtoms(); i++)	  if ((*qEnumerator)[i] == 0) {
			if (n == 0) n++; else (*qEnumerator)[i] = -1;
			//n=getAtom(i)->astereo;
			//(*qEnumerator)[n]=i;
		};
	};
};


int indexOfInteger(int value, const std::vector<int>& dataVector) {
	int result = -1;
	for (int i = 0; i < dataVector.size(); i++) if (dataVector[i] == value) {
		result = i;
		break;
	}
	return result;
}

/*

void TSimpleMolecule::extractFragment(int atomNo, int topoDistance, int firstAtomCharge, TSimpleMolecule & saveStructure) {
	//it is assumed that defineAtomConn() already called
	std::vector<int> enumeratedList, distanceList;
	int i, j, k, n, n1, n2, distN;
	TSingleAtom * at;
	TSingleAtom * atomTemp;
	TSingleBond * sb;

	saveStructure.clear();
	at = getAtom(atomNo)->clone();
	if (firstAtomCharge > 0) at->na = firstAtomCharge;
	saveStructure.addAtom(at);
	enumeratedList.push_back(atomNo);
	distanceList.push_back(0);

	//atom extraction
	for (i = 0; i < topoDistance; i++) {
		distN = distanceList.size();  //store old count for distance. Higher value will be filled with higher distances
		for (j = 0; j < distN; j++) if (distanceList[j] == i) {       //only for last distances
			n = enumeratedList[j];  //atomNumber
			atomTemp = getAtom(n);
			for (k = 0; k < atomTemp->nb; k++) {
				n = atomTemp->ac[k];
				if (indexOfInteger(n, enumeratedList) < 0) {
					at = getAtom(n)->clone();
					saveStructure.addAtom(at);
					enumeratedList.push_back(n);
					distanceList.push_back(i + 1);
				};
			};

		};
	};  //end atom collection
	//collect bonds
	for (i = 0; i < nBonds(); i++) {
		n1 = getBond(i)->at[0];
		n2 = getBond(i)->at[1];
		n1 = indexOfInteger(n1, enumeratedList);
		if (n1 >= 0) n2 = indexOfInteger(n2, enumeratedList);
		if ((n1 >= 0) && (n2 >= 0)) {
			sb = getBond(i)->clone();
			sb->at[0] = n1;
			sb->at[1] = n2;
			saveStructure.addBond(sb);
		};
	};
	saveStructure.defineAtomConn();
};
*/

void TSimpleMolecule::extractFragment(int atomNo, int topoDistance, int firstAtomCharge, TSimpleMolecule & saveStructure) {  //, std::vector<adjustedlist> * bk) {
	//it is assumed that defineAtomConn() already called
	std::vector<int> * atomEnumerator = new std::vector<int>(nAtoms());
	std::vector<int> atomList;
	//std::vector<int> atomOrder;
	std::vector<int> accumulatedList;

	int i, j, k, n, n1, n2, distN, na;
	TSingleAtom * at;
	TSingleAtom * atomTemp;
	TSingleBond * sb;

	saveStructure.clear();
	for (i = 0; i < atomEnumerator->size(); i++) (*atomEnumerator)[i] = -1;
	atomList.push_back(atomNo);
	na = 0;
	(*atomEnumerator)[atomNo] = na;
	at = getAtom(atomNo)->clone();
	if (firstAtomCharge > 0) at->na = firstAtomCharge;
	saveStructure.addAtom(at);
	na++;

	for (i = 0; i < topoDistance; i++) {
		for (j = 0; j < atomList.size(); j++) {
			n = atomList[j];
			for (k = 0; k < getAtom(n)->nb; k++) {
				n1 = getAtom(n)->ac[k];
				if ((*atomEnumerator)[n1] < 0) {
				  (*atomEnumerator)[n1] = na;
				  na++;
				  accumulatedList.push_back(n1);
				  at = getAtom(n1)->clone();
				  saveStructure.addAtom(at);
                }   
			}
		}
		atomList.resize(accumulatedList.size());
		for (j = 0; j < accumulatedList.size(); j++) atomList[j] = accumulatedList[j];
		accumulatedList.clear();
	};

	for (i = 0; i < nBonds(); i++) {
		n1 = getBond(i)->at[0];
		n2 = getBond(i)->at[1];
		if (((*atomEnumerator)[n1] >= 0) && ((*atomEnumerator)[n2] >= 0)) {
			sb = getBond(i)->clone();
			sb->at[0] = (*atomEnumerator)[n1];
			sb->at[1] = (*atomEnumerator)[n2];
			saveStructure.addBond(sb);
		}
	};
	
	/*
	std::vector<int> enumeratedList, distanceList;
	int i, j, k, n, n1, n2, distN;
	TSingleAtom * at;
	TSingleAtom * atomTemp;
	TSingleBond * sb;

	saveStructure.clear();

	at = getAtom(atomNo)->clone();
	if (firstAtomCharge > 0) at->na = firstAtomCharge;
	saveStructure.addAtom(at);
	enumeratedList.push_back(atomNo);
	distanceList.push_back(0);

	//atom extraction
	for (i = 0; i < topoDistance; i++) {
		distN = distanceList.size();  //store old count for distance. Higher value will be filled with higher distances
		for (j = 0; j < distN; j++) if (distanceList[j] == i) {       //only for last distances
			n = enumeratedList[j];  //atomNumber
			atomTemp = getAtom(n);
			for (k = 0; k < atomTemp->nb; k++) {
				n = atomTemp->ac[k];
				if (indexOfInteger(n, enumeratedList) < 0) {
					at = getAtom(n)->clone();
					saveStructure.addAtom(at);
					enumeratedList.push_back(n);
					distanceList.push_back(i + 1);
				};
			};

		};
	};  //end atom collection
		//collect bonds
	for (i = 0; i < nBonds(); i++) {
		n1 = getBond(i)->at[0];
		n2 = getBond(i)->at[1];
		n1 = indexOfInteger(n1, enumeratedList);
		if (n1 >= 0) n2 = indexOfInteger(n2, enumeratedList);
		if ((n1 >= 0) && (n2 >= 0)) {
			sb = getBond(i)->clone();
			sb->at[0] = n1;
			sb->at[1] = n2;
			saveStructure.addBond(sb);
		};
	};
	*/
	delete atomEnumerator;
	saveStructure.defineAtomConn();
};


TSimpleMolecule * TSimpleMolecule::extractFragment(int atomN, std::vector<int> * enumerator) {
  std::vector<int> list(listarSize());
  std::vector<int> inverseList(listarSize());
  int i,j,k,nA;
  bool test;
  TSimpleMolecule * result=NULL;
  TSingleAtom * sa;
  TSingleBond * sb;

  if ((atomN < 0) || (atomN >= nAtoms())) return result;
  if (enumerator) for (i=0; i<enumerator->size(); i++) (*enumerator)[i]=-1;
  for (i=0; i<nAtoms(); i++) inverseList[i]=-1;
  test=makeFragment(nA,list,atomN,-1);
  if (nA > 1) for (i=0; i<(nA-1); i++) for (j=i+1; j<nA; j++) if (list[i] > list[j]) {
    k=list[i];
    list[i]=list[j];
    list[j]=k;
  };
  if (nA > 0) for (i=0; i<nA; i++) inverseList[list[i]]=i;
  result=new (TSimpleMolecule);
  for (i=0; i<nA; i++) {
    sa=this->getAtom(list[i])->clone();
    result->addAtom(sa);
    if (enumerator) (*enumerator)[list[i]]=i;
  }; 
  for (i=0; i<nBonds(); i++) if ((inverseList[getBond(i)->at[0]] != -1) && (inverseList[getBond(i)->at[1]] != -1)) {
    sb=getBond(i)->clone();
    sb->at[0]=inverseList[sb->at[0]];
    sb->at[1]=inverseList[sb->at[1]];
    result->addBond(sb);
  };
  return result;
};

bool TSimpleMolecule::extractFragment(int atomNo, int acyclicBondNo, TSimpleMolecule & smExtracted) const {
	int i, j, n, nNew, n1, n2;
	std::vector<int> oldToNewEnumerator, atomFound, atomGenerated;
	TSingleAtom * sa;
	TSingleBond * sb;

	if (getBond(acyclicBondNo)->db >= 2) return false;

	oldToNewEnumerator.resize(nAtoms());
	for (i = 0; i < oldToNewEnumerator.size(); i++) oldToNewEnumerator[i] = -1;
	if (atomNo == getBond(acyclicBondNo)->at[0]) n = getBond(acyclicBondNo)->at[1]; else if (atomNo == getBond(acyclicBondNo)->at[1]) n = getBond(acyclicBondNo)->at[0]; else return false;
	oldToNewEnumerator[n] = -2;
	smExtracted.clear();
	nNew = 0;
	atomFound.push_back(atomNo);
	//try to put first atomNo

	
	while (atomFound.size()>0) {
		/*
		for (i = 0; i < atomFound.size(); i++) for (j=0; j<getAtom(atomFound[i])->nb; j++) {
			n = getAtom(atomFound[i])->ac[j];
			if (oldToNewEnumerator[n] == -1) {
				oldToNewEnumerator[n] = nNew;
				nNew++;
				sa = getAtom(n)->clone();
				smExtracted.addAtom(sa);
				atomGenerated.push_back(n);
			}
		};
		*/
		for (i = 0; i < atomFound.size(); i++) {
			n = atomFound[i];
			if (oldToNewEnumerator[n] == -1) {
				oldToNewEnumerator[n] = nNew;
				nNew++;
				sa = getAtom(n)->clone();
				smExtracted.addAtom(sa);
				//atomGenerated.push_back(n);
				for (j = 0; j < getAtom(atomFound[i])->nb; j++) {
					n = getAtom(atomFound[i])->ac[j];
					if (oldToNewEnumerator[n] == -1) atomGenerated.push_back(n);
				}
			}
		};

		atomFound.clear();
		for (i = 0; i < atomGenerated.size(); i++) atomFound.push_back(atomGenerated[i]);
		atomGenerated.clear();
	}
	
	for (i = 0; i < nBonds(); i++) if ((oldToNewEnumerator[getBond(i)->at[0]] >= 0) && (oldToNewEnumerator[getBond(i)->at[1]] >= 0)) {
		sb = getBond(i)->clone();
		sb->at[0] = oldToNewEnumerator[sb->at[0]];
		sb->at[1] = oldToNewEnumerator[sb->at[1]];
		smExtracted.addBond(sb);
	};
	smExtracted.defineAtomConn();
	return true;
}

void TSimpleMolecule::createFragmentList(std::vector<TSimpleMolecule *> & list) {
//makes a list of fragments from single structure....

  std::vector<int> enumerator(nAtoms());
  std::vector<int> enumeratorStore(nAtoms());
  TSimpleMolecule * sm;
  int i,j,n;

  for (i=0; i<enumeratorStore.size(); i++) enumeratorStore[i]=-1;
  n=0;
  while (n>=0) {
    n=-1;
    for (i=0; i<enumeratorStore.size(); i++) if (enumeratorStore[i] == -1) {
      n=i;
      break;
	};
    if (n>=0) {
      sm=extractFragment(n,&enumerator);
	  for (i=0; i<enumeratorStore.size(); i++) if (enumerator[i] >= 0) enumeratorStore[i]=1;
	  sm->defineAtomConn();
		/*
        SM.DetermineFormula;
        CV:=CMolecule(Self.ClassType);
        EM:=CV.Create;
        EM.MoleculeCopy(SM);
		*/
      list.push_back(sm);
	};
  }; //until N<0;


  if (list.size()>1) for (i=0; i<(list.size()-1); i++) for (j=i+1; j<list.size(); j++) {
	n=0;
	if (list[i]->nAtoms() > list[j]->nAtoms()) n= 1; else
	if (list[i]->nAtoms() < list[j]->nAtoms()) n=-1; else
	if (list[i]->nBonds() > list[j]->nBonds()) n= 1; else
	if (list[i]->nBonds() < list[j]->nBonds()) n=-1; else
	if (list[i]->getMolWeight() > list[j]->getMolWeight()) n=-1; else
	if (list[i]->getMolWeight() < list[j]->getMolWeight()) n= 1;
    //Molweight - inverse to remove heavy metals....
    //if TSimpleMolecule(List[I]).MolWeight>TSimpleMolecule(List[J]).MolWeight then N:=-1 else
    //if TSimpleMolecule(List[I]).MolWeight<TSimpleMolecule(List[J]).MolWeight then N:= 1 else N:=0;
    if (n < 0) {
	  sm=list[i];
	  list[i]=list[j];
	  list[j]=sm;
	};
  };

};


int TSimpleMolecule::processCharges() {
//assumed single-connected molecule...
  int i,j,k;
  TSimpleMolecule em;
  int result=0;

  em.moleculeCopy(*this);
  em.defineAtomConn();
  for (i=0; i<em.nAtoms(); i++) {
    result=result+em.getAtom(i)->nc;
    if ((em.getAtom(i)->nc == 1) && isHetero(em.getAtom(i)->na)) {
      //try to search for hydrogen for deletion...
	  for (j=0; j<em.nBonds(); j++) if ((em.getBond(j)->at[0] == i) || (em.getBond(j)->at[1] == i)) {
        k=em.getBond(j)->at[0];
        if (k == i) k=em.getBond(j)->at[1];
        if ((em.getAtom(k)->na == 1) && (em.getAtom(k)->nb == 1)) {
          em.getAtom(k)->enumerator=65535;
          em.getAtom(k)->nc=0;
          break;
		};
	  };
    } else em.getAtom(i)->nc=0;
  };
  for (i=em.nAtoms()-1; i>=0; i--) if (em.getAtom(i)->enumerator == 65535) em.deleteAtom(i);
  this->moleculeCopy(em);
  this->defineAtomConn();
  return result;
};


int  TSimpleMolecule::makeUniqueSymbol(const std::vector<adjustedlist> & bondConnection) {
/*
  typedef struct liArray {
    long int data[MAXBONDS];
  } liArray;
  typedef struct wArray {
    unsigned short data[MAXBONDS];
  } wArray;
*/
  int result=0;
  int i,j,k,m,n,ii;
//  BK:TBondConnection;
  std::vector<unsigned short> define(listarSize());
  long int uk1,uk2,uk,nnn;
  bool test;
  std::vector<long int> uk3(listarSize());
  std::vector<long int> uk4(listarSize());
  float r1;
  double r,rr1,rr2;


  if (nAtoms() == 0) return result;
  if (nAtoms() == 1) {  //Partial case - single - atom. to save charge, radical
	result=getAtom(0)->rl*100+getAtom(0)->nc;
    return result;
  };

  for (i=0; i<nAtoms(); i++) {
    uk3[i]=0; uk4[i]=0;
    for (j=0; j<nAtoms(); j++) define[j]=0; //all atoms have not been took into consideration yet}
    define[i]=1;                       //exept I-th}
    k=0;
	test=false;
    while (! test) {
      k++;                           //K-counter of currently hot neighbour sphere
      test=true;
      for (j=0; j<nAtoms(); j++) if (define[j] == k) { //if atom is hot...}
        uk1=getAtom(j)->na;
        uk=getAtom(j)->nc;
        uk1=uk1 ^ uk;
        uk=getAtom(j)->iz;
        uk1=uk1 ^ uk;
        //attributes, stored in byte
        uk1=(uk1*16384) / k;
        uk3[i]=uk3[i]+uk1;    //atom's attributes accumulation
        uk1=0;
        for (ii=0; ii<bondConnection[j].nb; ii++) { //packing of bond's attributes
		  m=bondConnection[j].adjusted[ii];;
          nnn=getBond(m)->tb;;
          if ((getBond(m)->db == 2) || (getBond(m)->db == 3)) nnn=8;  //!!!! MUST BE 4!!!!
          n=getBond(m)->at[0]; if (n == j) n=getBond(m)->at[1];
          uk2=getAtom(n)->na;
          uk=getAtom(n)->nc;
          uk2=uk2 ^ uk;
          uk=getAtom(n)->iz;
          uk2=uk2 ^ uk;
          uk1=uk1+nnn*uk2;
		};
        uk1=uk1 / k;
        uk4[i]=uk4[i]+uk1;    //bond's attributes accumulation
        for (ii=0; ii<getAtom(j)->nb; ii++) { //labelling of next neighbour sphere
		  m=getAtom(j)->ac[ii];
		  if (define[m] == 0) { define[m]=k+1; test=false; };
		};
	  };
	}; // while until Test; {until all atoms will be defined}
  };
  r=0;
  for (i=0; i<nAtoms(); i++) {
    rr1=uk3[i];
    rr2=uk4[i];
    r=r+rr1*rr2; //{formation of final code}
  };
  r1=r;
  memcpy(&uk,&r1,4);
  result=uk;


  return result;
};


int TSimpleMolecule::chargeConversion(int atn) {

/*for atom's number ATN in array ATOM return some connected with charge values:
 3 - if radical label present.
 2 - if charge <0
 1 - if charge >0
 0 - charge=0*/
  int result=0;

  if (getAtom(atn)->rl != 0) result=3; else //[ATN].NC>=9 then begin ChargeConversion:=3; Exit; end;
  if (getAtom(atn)->nc < 0) result=2; else
  if (getAtom(atn)->nc > 0) result=1;
  return result;
};

int TSimpleMolecule::valencyConversion(int atn) {

/*for atom's number ATN in array ATOM returns some connected with valency value:
 =2-Valence of ATN is less, then usual
 =1-Valence of ATN is more, then usual
 =0-usual valency*/

int k,k1,k2;
int result=0;

    //Default hydrogen
  k1=getAtom(atn)->nv;
  k1=k1-getAtom(atn)->currvalence-abs(getAtom(atn)->nc)-getAtom(atn)->rl;
  if (k1 < 0) k1=0;
  k2=hVal[getAtom(atn)->na];
  k2=k2-getAtom(atn)->currvalence-abs(getAtom(atn)->nc)-getAtom(atn)->rl;
  if (k2 < 0) k2=0;
  if (k1 == k2) result=0; else if (k1 < k2) result=1; else result=2;
  return result;
};

int  TSimpleMolecule::allAtAtom(int atn) {
/*Define a digital representation of atom, they include:
 a) Position of atom in the Periodic System
 b) Its charge
 c) Its valency*/
  int b1,b2,b3,w;

  b1=encoder(getAtom(atn)->na);
  b2=chargeConversion(atn);
  b3=valencyConversion(atn);
  w=b3;
  w=w << 2;
  w=w+b2;
  w=w << 5;
  w=w+b1;
  //Addition from 1 June 2001}
  if (getAtom(atn)->rl != 0) w= ~ w;
  return w;
};

int  TSimpleMolecule::bondConversion(int bnb) {
//generate a code for bond's nuber BNB in array BOND. It is took into consi deration bond type and cycle size}
  int b1,b2;

  if (getBond(bnb)->tb >= 9) b1=1; else if (getBond(bnb)->tb > 4)  b1=0; else b1=getBond(bnb)->tb;
  b2=7;
  switch (getBond(bnb)->db) {
    case 0: {
      b2=5;
	  break;
      };
	case 1: {
      b2=0;
	  break;
	  };
	case 2: {
      b1=4; b2=1;
	  break;
	};
	case 3: {
      b1=4; b2=2;
	  break;
	};
	case 4:{
      b2=3;
	  break;
    };
	case 5: {
	  b2=4;
	  break;
	};
    case 6: {
      b2=6;
	  break;
	};
	b2=7;
  };
  b2=b2 << 2;
  b2=b2+b1;
  return b2;
};

int compareFragments(const twoSphereRecord s2, int itemIndex, const std::vector<twoSphereRecord> fragments) {
  twoSphereRecord s1;
  int result=0;

  s1=fragments[itemIndex];
  if (s1.c1 > s2.c1) result= 1; else
  if (s1.c1 < s2.c1) result=-1; else 
  if (s1.c2 > s2.c2) result= 1; else
  if (s1.c2 < s2.c2) result=-1; else result=0;
  return result;
};


void quickSortFragments(int iLo, int iHi, std::vector<twoSphereRecord> & fragments) {
  int lo, hi;
  twoSphereRecord mid,t1,t2;

  lo=iLo;
  hi=iHi;
  mid=fragments[(lo + hi) /2];
  while (lo<=hi) {
    while (compareFragments(mid,lo,fragments) < 0) lo++;
    while (compareFragments(mid,hi,fragments) > 0) hi--;
    if (lo <= hi) {
      t1=fragments[lo];
	  t2=fragments[hi];
	  fragments[lo]=t2;
	  fragments[hi]=t1;
	  lo++;
	  hi--;
	};
  };
  if (hi > iLo) quickSortFragments(iLo,hi,fragments);
  if (lo < iHi) quickSortFragments(lo,iHi,fragments);
};


void TSimpleMolecule::processCodeList() {
  int i;
  twoSphereRecord d1,d2;
  std::vector<twoSphereRecord> temp;

  if (fragmentCodeList.size() <= 1) return;

  quickSortFragments(0,fragmentCodeList.size()-1,fragmentCodeList);
  temp.push_back(fragmentCodeList[0]);
  for (i=1; i<fragmentCodeList.size(); i++) {
    d1=fragmentCodeList[i];
	if (compareFragments(d1,i-1,fragmentCodeList) != 0) temp.push_back(fragmentCodeList[i]);
  };
  fragmentCodeList.resize(temp.size());
  for (i=0; i<temp.size(); i++) fragmentCodeList[i]=temp[i];
};


int  TSimpleMolecule::makeSecondSymbol(const std::vector<adjustedlist> & bondConnection){
  int result=0;
  int i,j,n,k;
  std::vector<int>a(listarSize());
  std::vector<int>b(listarSize());
  std::vector<int>ws;
  std::vector<long int> ws1(listarSize());
  std::vector<long int> ws2(listarSize());
  long int un,ii;
  bool test;
  twoSphereRecord sr;

  if (nAtoms() == 0) return result;

  n=listarSize();
  if (n < 10) n=10;
  ws.resize(n);

  
  for (j=0; j<nAtoms(); j++) a[j]=allAtAtom(j);
  for (j=0; j<nBonds(); j++) b[j]=bondConversion(j);
  fragmentCodeList.resize(nAtoms());
  for (j=0; j<nAtoms(); j++) {
    fragmentSecond(2,j,0,a,b,bondConnection,ws);
    ws1[j]=ws[0];
    ws2[j]=ws[1];
    sr.c1=ws[0];
    sr.c2=ws[1];
    fragmentCodeList[j]=sr;
  };
  processCodeList();
  if (nAtoms() > 1) for (i=0; i<(nAtoms()-1); i++) for (j=i+1; j<nAtoms(); j++) if ((ws2[i] > ws2[j]) || ((ws2[i] == ws2[j]) && (ws1[i] > ws1[j]))) {
    n=ws2[j];
    ws2[j]=ws2[i];
    ws2[i]=n;
    n=ws1[j];
    ws1[j]=ws1[i];
    ws1[i]=n;
  };
  n=nAtoms(); //N:=FAtom.NAtoms;
  i=0;  //I:=0;
  if (nAtoms() > 1) while (i < (n-1)) { //if FAtom.NAtoms>1 then repeat
    //I:=I+1;
    test=(ws2[i] == ws2[i+1]);
    if (test) {
      for (j=i; j<n; j++) ws2[j]=ws2[j+1];
      for (j=i; j<n; j++) ws1[j]=ws1[j+1];
      i--; //I:=I-1;
      n--;  //N:=N-1;
	};
	i++;
  }; //until I>=(N-1);
  

  un=0;
  for (i=0; i<n; i++) {
    k=ws1[i];
    un=un ^ (k << 16);
    k=ws2[i];
    un=un ^ k;
  };
  result=un;
  return result;
};


void TSimpleMolecule::calculateAllIndeces(){
  TSimpleMolecule smCalc;
  int i,j,k,nH,chargeTotal;
  std::vector<int> qhydr;
  std::vector <TSimpleMolecule *> fragList;
  std::vector<adjustedlist> bk(listarSize());
  int uc1, uc2, uc3, ch;
  double dw, mwTotal;
  float  r;
  bool test,testSingleH;
  fragmentCode data;
  int k1,k2,k3;
  //bool testQuery;

  defineAtomConn();
  allAboutCycles();
  smCalc.moleculeCopy(*this);
  for (i=0; i<nBonds(); i++) if (smCalc.getBond(i)->tb == 4) smCalc.getBond(i)->db=3;
  smCalc.unC1=0;
  smCalc.unC2=0;
  smCalc.unC3=0;
  smCalc.inChIKey="";
  smCalc.makeStandardStructure();
  for (i=0; i<smCalc.nAtoms(); i++) smCalc.getAtom(i)->astereo=smCalc.getNH(i);  //VERY important !!!!!!
  smCalc.allAboutCycles();
  smCalc.tautomerAnalize();
  smCalc.allAboutCycles();
  for (i=0; i<smCalc.nBonds(); i++) if ((smCalc.getBond(i)->db == 2) || (smCalc.getBond(i)->db == 3)) smCalc.getBond(i)->tb=4;
  smCalc.defineAtomConn();
  qhydr.resize(smCalc.listarSize());
  for (i=0; i<smCalc.nAtoms(); i++) if (smCalc.getAtom(i)->na == 1) {
    if (smCalc.getAtom(i)->iz == 1) {
       smCalc.getAtom(i)->na=DEUTERIUM_ATOM;
       smCalc.getAtom(i)->iz=0;
	};
    if (smCalc.getAtom(i)->iz == 2) {
      smCalc.getAtom(i)->na=TRITIUM_ATOM;
      smCalc.getAtom(i)->iz=0;
	};
  };
  smCalc.fIOPT11=true;
  smCalc.fIOPT12=1;
  smCalc.fIOPT13=true;
  //testSingleH=((smCalc.nAtoms() == 1) && (smCalc.getAtom(0)->na == 1));
  smCalc.removeHydrogen(&qhydr,NULL); //{Connection is calculated inside RemoveHydr}
  smCalc.atomBondChange();          //semipolar bond conversion}
  for (i=0; i<smCalc.nBonds(); i++) if (smCalc.getBond(i)->tb >= 9) smCalc.getBond(i)->tb=1;
  smCalc.defineAtomConn();
  mwTotal=0;
  for (i=0; i<smCalc.nAtoms(); i++) {
    nH=smCalc.getAtom(i)->nv;
	nH=nH-smCalc.getAtom(i)->currvalence-abs(smCalc.getAtom(i)->nc)-smCalc.getAtom(i)->rl;
    if (nH < 0) nH=0;
	mwTotal=mwTotal+aMass[smCalc.getAtom(i)->na]+nH*aMass[1];;
	if (smCalc.getAtom(i)->astereo > nH) smCalc.getAtom(i)->nv=smCalc.getAtom(i)->astereo+smCalc.getAtom(i)->currvalence+abs(smCalc.getAtom(i)->nc)+smCalc.getAtom(i)->rl;
  };
  smCalc.createFragmentList(fragList);


  chargeTotal=0;
  mixtureCodeList.clear();
  for (i=0; i<fragList.size(); i++) {
    k=fragList[i]->processCharges();  //removes all charges and make molecule neutral, total charge is stored in K
    chargeTotal=chargeTotal+k;
    //Ne goditsja for Salt Part. Bad Molweight, also as
	fragList[i]->defineAtomConn();
	fragList[i]->defineBondConn(bk);

    uc1=fragList[i]->makeUniqueSymbol(bk)-k;
    uc2=fragList[i]->makeSecondSymbol(bk)+k;
    ch=k;
    dw=0;
	for (j=0; j<fragList[i]->nAtoms(); j++) {
	  test=true;
	  if (fragList[i]->getAtom(j)->na == 1) {
		int kk=0;
		for (int nn=0; nn<fragList[i]->nBonds(); nn++) {
		  if (fragList[i]->getBond(nn)->at[0] == j) kk++;
		  if (fragList[i]->getBond(nn)->at[1] == j) kk++;
		};
		test=(kk == 0);
	  };
	  if (test) dw=dw+aMass[fragList[i]->getAtom(j)->na]+fragList[i]->getAtom(j)->astereo*aMass[1];
	};
    if (k < 0) dw=dw-k*aMass[1];
	//mwTotal=mwTotal+dw;
	r=dw;
	memcpy(&uc3,&r,4);
//	inChIKey=convertMolToInChIKey(*fragList[i]);
    test=false;
    for (j=0; j<i; j++) if ((fragList[j]->unC1 == uc1) && (fragList[j]->unC2 == uc2) && (fragList[j]->unC3 == uc3) && (fragList[j]->inChIKey == inChIKey)) {
	  test=true;
	  break;
	};
	if (! test) {  //unique fragment
  	  fragList[i]->unC1=uc1;
	  fragList[i]->unC2=uc2;
	  fragList[i]->unC3=uc3;
	  fragList[i]->inChIKey=inChIKey;
	  data.unc1=uc1;
	  data.unc2=uc2;
	  data.unc3=uc3;
	  setInChIKeyString(inChIKey,data);
	  mixtureCodeList.push_back(data);
	};
  };

  k1=0; k2=0; k3=0;
  for (i=0; i<mixtureCodeList.size(); i++) {
	data=mixtureCodeList[i];
    if (i == 0) {
	  k1=data.unc1;
	  k2=data.unc2;
	} else {
      k1=k1 ^ data.unc1;
      k2=k2 ^ data.unc2;
	};
  };
  unC1=k1;
  unC2=k2;
  r=mwTotal;
  memcpy(&unC3,&r,4);
//  inChIKey=convertMolToInChIKey(*this);
  for (i=0; i<fragList.size(); i++) delete(fragList[i]);
};


double	TSimpleMolecule::getMolWeight(bool calc){
  double result=0;
  int j,k;

  if(calc){defineAtomConn();};
  for (j=0; j<nAtoms(); j++) { 
  	k=getAtom(j)->nv;
	k=k-getAtom(j)->currvalence-abs(getAtom(j)->nc)-getAtom(j)->rl;
    if (k < 0) k=0;
	result=result+aMass[getAtom(j)->na]+k*aMass[1];
  };

  return result;
}


void TSimpleMolecule::extractFragmentFast(const std::vector<int> aList, TSimpleMolecule & smExtracted) {
//Is uses AEnumerator-array will be destored on output
  int i,n,k1,k2;
  TSingleAtom * sa;
  TSingleBond * sb;

  for (i=0; i<nAtoms(); i++) getAtom(i)->enumerator=0;
  smExtracted.clear();
  for (i=0; i<aList.size(); i++) {
    n=aList[i];
    sa=getAtom(n)->clone();
    smExtracted.addAtom(sa);
	getAtom(n)->enumerator=i+1;
  };
  for (i=0; i<nBonds(); i++) {
    k1=getBond(i)->at[0];
    k2=getBond(i)->at[1];
	if ((getAtom(k1)->enumerator > 0) && (getAtom(k2)->enumerator > 0)) {
	  sb=getBond(i)->clone();
      sb->at[0]=getAtom(k1)->enumerator-1;
      sb->at[1]=getAtom(k2)->enumerator-1;
      smExtracted.addBond(sb);
	};
  };
  smExtracted.defineAtomConn();
  //After return Bond Array save topology (rings...)
};


bool TSimpleMolecule::addMoleculeNoScale(TSimpleMolecule * molecule) {
  int i,naStore;
  TSingleAtom * sa;
  TSingleBond * sb;
  bool result=false;

  //if ((nAtoms()+molecule->nAtoms()) > NATOMSMAXSMBIG) return result;
  //if ((nBonds()+molecule->nBonds()) > NBONDSMAXSMBIG) return result;
  if (molecule->nAtoms() == 0) {
    result=true;
    return result;
  };

  naStore=nAtoms();
  for (i=0; i<molecule->nAtoms(); i++) {
	sa=molecule->getAtom(i)->clone();
    //FAtom.Append(SA,1);
	addAtom(sa);
  };
  if (molecule->nBonds() > 0) for (i=0; i<molecule->nBonds(); i++) {
	sb=molecule->getBond(i)->clone();
	sb->at[0]=sb->at[0]+naStore;
	sb->at[1]=sb->at[1]+naStore;
    addBond(sb);
  };
  defineAtomConn();
  result=true;
  return result;
};

int	TSimpleMolecule::getNumberOfAtoms(int na) {
  int result=0;
  int i,k;
  TSingleAtom * atom;

  for (i=0; i<nAtoms(); i++) {
	if (getAtom(i)->na == na) result++;
	if (na == 1) {
  	  atom=getAtom(i);
	  k=atom->nv;
	  k=k-(atom->currvalence)+(atom->nc*TSingleAtom::chargeDeltaValency(atom->na))-(atom->rl);
	  if (k < 0) k=0;
	  result=result+k;
	};
  };
  return result;
};

bool TSimpleMolecule::formulaIdentical(TSimpleMolecule * other){
  bool result=true;
  bool hasQuery=false;
  int thisCount[NELEMMCDL];
  int otherCount[NELEMMCDL];
  int i,n,k;
  TSingleAtom * atom;

  for (i=0; i<NELEMMCDL; i++) {
	thisCount[i]=0;
	otherCount[i]=0;
  };

  for (i=0; i<nAtoms(); i++) {
    atom=getAtom(i);
	n=atom->na;
	if (n <= NELEMMCDL) thisCount[n-1]++; else hasQuery=true;
    k=atom->nv;
//	k=k-(atom->currvalence)+(atom->nc*TSingleAtom::chargeDeltaValency(atom->na))-(atom->rl);
	k=k-(atom->currvalence)-abs(atom->nc)-(atom->rl);
	if (k < 0) k=0;
	thisCount[0]=thisCount[0]+k;
  };
  for (i=0; i<other->nAtoms(); i++) {
    atom=other->getAtom(i);
	n=atom->na;
	if (n <= NELEMMCDL) otherCount[n-1]++; else hasQuery=true;
    k=atom->nv;
//	 k=k-(atom->currvalence)+(atom->nc*TSingleAtom::chargeDeltaValency(atom->na))-(atom->rl);
	k=k-(atom->currvalence)-abs(atom->nc)-(atom->rl);
	if (k < 0) k=0;
	otherCount[0]=otherCount[0]+k;
  };
  for (i=0; i<NELEMMCDL; i++) if (thisCount[i] != otherCount[i]) {
	result=false;
	break;
  };
  return result;
};

std::string	TSimpleMolecule::getMolformula(bool isHTML) const {
  std::string result="";
  int atomCount[NELEMMCDL];
  std::string atomCode[NELEMMCDL];
  int i,j,n,k;
  const TSingleAtom * atom;
  bool hasQuery=false;
  std::string s;

  for (i=0; i<NELEMMCDL; i++) {
	atomCount[i]=0;
	atomCode[i]="";
  };

  for (i=0; i<nAtoms(); i++) {
    atom=getAtom(i);
	n=atom->na;
	if (n <= NELEMMCDL) atomCount[n-1]++; else hasQuery=true;
    k=atom->nv;
//	k=k-(atom->currvalence)+(atom->nc*TSingleAtom::chargeDeltaValency(atom->na))-(atom->rl);
	k=k-(atom->currvalence)-abs(atom->nc)-(atom->rl);
	if (k < 0) k=0;
	atomCount[0]=atomCount[0]+k;
  };

  for (i=0; i<NELEMMCDL; i++) if (atomCount[i] > 0) {
	atomCode[i]=aSymb[i+1];
	if (atomCount[i] > 1) {
	  //if (isHTML) atomCode[i]=atomCode[i]+"<sub>";
	  if(isHTML){ 	atomCode[i]=atomCode[i]+"<sub>"+intToStr(atomCount[i])+"</sub>";}
	  else{			atomCode[i]=atomCode[i]+intToStr(atomCount[i]);};
  }};
  if (atomCount[5] > 0) {   //Carbon
	result=result+atomCode[5];
	atomCode[5]="";
	atomCount[5]=0;
  };
  if (atomCount[0] > 0) {   //Hydrogen
	result=result+atomCode[0];
	atomCode[0]="";
	atomCount[0]=0;
  };
  //Compact
  n=0;
  for (i=0; i<NELEMMCDL; i++) if (atomCount[i]>0) {
	atomCount[n]=atomCount[i];
	atomCount[i]=0;
	atomCode[n]=atomCode[i];
	atomCode[i]="";
	n++;
  };
  for (i=0; i<(n-1); i++) for (j=i+1; j<n; j++) if (atomCode[i] > atomCode[j]) {
	s=atomCode[i];
	atomCode[i]=atomCode[j];
	atomCode[j]=s;
  };
  for (i=0; i<n; i++) result=result+atomCode[i];
  return result;
};


void TSimpleMolecule::init(){
  refofs=NULL;
  fIOPT10=true;
  fIOPT11=true;   //isotope
  fIOPT12=true;   //stereo bond change
  fIOPT13=true;   //semipolar bond as double
  unC1=0;
  unC2=0;
  unC3=0;
  inChIKey="";
  userData=0;
  chiral=0;
  fScreenData[0]=0;
  fScreenData[1]=0;
  fScreenData[2]=0;
  fScreenData[3]=0;
  fScreenData[4]=0;
  fAtom.clear();
  fBond.clear();
};

TSimpleMolecule::TSimpleMolecule() {
  init();
};

bool TSimpleMolecule::isSP3Atom(int atomNo, const std::vector<adjustedlist> & bondConnection) const {
  bool result=false;
  int n,i,j,k,an;
  int ns, nd, nt;

  ns=0;
  nd=0;
  nt=0;
  for (i = 0; i < bondConnection[atomNo].nb; i++) {
	k=bondConnection[atomNo].adjusted[i];
	k=getBond(k)->tb;
	if ((k == 1) || (k == 9) || (k == 10) || (k == 11)) k=1;
	if (k == 1) ns++; else if (k == 2) nd++; else if (k == 3) nt++;
  }
//  if ((nd == 0) && (nt == 0)) { //07.05.2021
	  //result = true;
	  //return result;
  //}
  n = getAtom(atomNo)->na;
  if (n == 6) result=((nd == 0) && (nt == 0)); else
  if ((n == 7) || (n==8)) {
	  result = ((nd == 0) && (nt == 0));
	  if (result) {  //checking if second-neighbours have no double bonds
		  nd = 0;
		  for (i = 0; i < bondConnection[atomNo].nb; i++) {
			  k = bondConnection[atomNo].adjusted[i];
			  an = getBond(k)->at[0];
			  if (an == atomNo) an = getBond(k)->at[1];
			  for (j = 0; j < bondConnection[an].nb; j++) {
				  k = bondConnection[an].adjusted[j];
				  k = getBond(k)->tb;
				  if (k == 2) nd++;
			  };
			  if (nd > 0) break;
		  };
		  result = (nd == 0);
	  };
  } else if (n == 15) result=((nd == 1) && (ns == 3)); else
  if (n == 16) result=((nd==2) && (ns == 2));
  return result;
};

bool TSimpleMolecule::isSP2Atom(int atomNo, const std::vector<adjustedlist> & bondConnection) const {
  bool result=false;
  int n,i,j,k,an;
  int ns, nd, nt;

  ns=0;
  nd=0;
  nt=0;
  for (i = 0; i < bondConnection[atomNo].nb; i++) {
	k=bondConnection[atomNo].adjusted[i];
	k=getBond(k)->tb;
	if ((k == 1) || (k == 9) || (k == 10) || (k == 11)) k=1;
	if (k == 1) ns++; else if (k == 2) nd++; else if (k == 3) nt++;
  }
  n = getAtom(atomNo)->na;
  if (nd > 0) {  //07.05.2021
	  result = true;
	  return result;
  }
  if (n == 6) result=(nd == 1); else
	  if ((n == 7) || (n == 8)) {
		  result = (nd > 0);
		  if (!result) {
			  nd = 0;
			  for (i = 0; i < bondConnection[atomNo].nb; i++) {
  			      k = bondConnection[atomNo].adjusted[i];
				  an = getBond(k)->at[0];
			      if (an == atomNo) an = getBond(k)->at[1];
				  for (j = 0; j < bondConnection[an].nb; j++) {
					  k = bondConnection[an].adjusted[j];
					  k = getBond(k)->tb;
					  if (k == 2) nd++;
				  };
				  if (nd > 0) break;
			  };
			  result = (nd > 0);
		  };
	  };
  return result;
};

bool TSimpleMolecule::isSPAtom(int atomNo, const std::vector<adjustedlist> & bondConnection) {
  bool result=false;
  int n,i,k;
  int ns, nd, nt;

  ns=0;
  nd=0;
  nt=0;
  for (i = 0; i < bondConnection[atomNo].nb; i++) {
	k=bondConnection[atomNo].adjusted[i];
	k=getBond(k)->tb;
	if ((k == 1) || (k == 9) || (k == 10) || (k == 11)) k=1;
	if (k == 1) ns++; else if (k == 2) nd++; else if (k == 3) nt++;
  }
  n = getAtom(atomNo)->na;
  if (n == 6) result=(nt == 1); else
  if (n == 7) result=(nt == 1); 
  return result;
};

bool TSimpleMolecule::isAromaticAtom(int atomNo, const std::vector<adjustedlist> & bondConnection) {
  //It is assumed that All AboutCycle was called early
  bool result=false;
  int i,k;

  for (i=0; i<bondConnection[atomNo].nb; i++) {
	k=bondConnection[atomNo].adjusted[i];
	k=getBond(k)->db;
	result=((k == 2) || (k == 3));
	if (result) break;
  };
  return result;
};

bool TSimpleMolecule::isPlanarRing(int cycleSize,  const std::vector<int> & bondList, const std::vector<adjustedlist> & bondConnection) {
  //It is assumed that bond list is caconized and deterine ring path
  bool result=false;
  int i;
  int nSP3,an,n;

  nSP3=0;
  if (cycleSize <= 3) return true; else {
	an=getBond(bondList[0])->at[0];
	if (isSP3Atom(an, bondConnection)) nSP3++;
	for (i=0; i<(cycleSize-1); i++) {   // last bond is not treated because of it closes the ring
      n=getBond(bondList[i])->at[0];
	  if (n == an) n=getBond(bondList[i])->at[1];
	  if (isSP3Atom(n, bondConnection)) nSP3++;
	  an=n;
	};
  };
  result=(nSP3 < 2);
  return result;
};

void TSimpleMolecule::zeroMolDynProperties() {
	for (int i = 0; i < nAtoms(); i++) {
		getAtom(i)->chirality = 0;
		getAtom(i)->aromatic = 0;
		getAtom(i)->hybridization = 0;
		getAtom(i)->rKind1 = 0;
		getAtom(i)->rSize1 = 0;
		getAtom(i)->rKind2 = 0;
		getAtom(i)->rSize2 = 0;
		getAtom(i)->rKind3 = 0;
		getAtom(i)->rSize3 = 0;
		getAtom(i)->rKind4 = 0;
		getAtom(i)->rSize4 = 0;
		getAtom(i)->rKind5 = 0;
		getAtom(i)->rSize5 = 0;
		getAtom(i)->rKind6 = 0;
		getAtom(i)->rSize6 = 0;
		getAtom(i)->rejectInList = false;
		if (getAtom(i)->atomList) getAtom(i)->atomList->clear();
	};
};

void TSimpleMolecule::fillMolDynStructureProperties() {
  //all about cycles must be called early
  int i,j,an,n,nn;
  int cycleSize;
  std::vector<adjustedlist> bk(listarSize());
  std::vector<int>bondList(listarSize());
  std::vector<bool>bondChecked(nBonds());
  bool isPlanar,isSP3,isSP2,isSP;
  int at1, at2;
  bool test;

  /*   Must be called early
  for (i=0; i< nAtoms(); i++) {
	getAtom(i)->chirality=0;
	getAtom(i)->aromatic=0;
	getAtom(i)->hybridization=0;
	getAtom(i)->rKind1=0;
	getAtom(i)->rSize1=0;
	getAtom(i)->rKind2=0;
	getAtom(i)->rSize2=0;
	getAtom(i)->rKind3=0;
	getAtom(i)->rSize3=0;
	getAtom(i)->rKind4=0;
	getAtom(i)->rSize4=0;
	getAtom(i)->rKind5=0;
	getAtom(i)->rSize5=0;
	getAtom(i)->rKind6=0;
	getAtom(i)->rSize6=0;
	getAtom(i)->rejectInList=false;
	if (getAtom(i)->atomList) getAtom(i)->atomList->clear();
  }; */
  if (nBonds() == 0) return;
  zeroMolDynProperties();
  defineAtomConn();
  
  defineBondConn(bk);
  // This code is calculated in the prior call of allabout cycles-not so
  for (i = 0; i < bondChecked.size(); i++) bondChecked[i] = false;
  for (i=0; i<nBonds(); i++) if (! bondChecked[i]) {
	  at1 = getBond(i)->at[0];
	  at2 = getBond(i)->at[1];
	  test = ((getAtom(at1)->nb == 1) || (getAtom(at2)->nb == 1));
	  if (!test) {
		  vaweBond(i, bk, cycleSize, bondList, CYCLE_MAX_SIZE); //Is I-th bond cyclic?
		  bondChecked[i] = true;
		  if (cycleSize > 0) {       //Yes, cyclic
			  canonizeCycle(cycleSize, bondList);
			  for (int j = 0; j < cycleSize; j++) {
				  bondChecked[bondList[j]] = true;
				  isPlanar = isPlanarRing(cycleSize, bondList, bk);
				  an = getBond(bondList[j])->at[0];
				  getAtom(an)->setRingProperty(cycleSize, isPlanar);
				  an = getBond(bondList[j])->at[1];
				  getAtom(an)->setRingProperty(cycleSize, isPlanar);
			  };
		  };
	  };
	  bondChecked[i]= true;
  };
  
  allAboutCycles(CYCLE_MAX_SIZE);  //changes from 17 MArch 2021
  
  for (i = 0; i < nBonds(); i++) if ((getBond(i)->db == 7) || (getBond(i)->db == 6)) {
	  //correct aromaticity. Definition differs from my
	  if (getBond(i)->db == 7) nn = 3; else nn = 2;
	  vaweBond(i, bk, cycleSize, bondList, CYCLE_MAX_SIZE); //Is I-th bond cyclic?
	  if (cycleSize > 0) {       //Yes, cyclic
		  bool testArom = true;
		  for (j = 0; j < cycleSize; j++) {
			  n = getBond(bondList[j])->at[0];
			  if (!isSP2Atom(n, bk)) {
				  if (isHetero(getAtom(n)->na) && (getAtom(n)->nb == hVal[getAtom(n)->na])) testArom = true; else {
					  testArom = false;
					  break;
				  };
			  }
			  n = getBond(bondList[j])->at[1];
			  if (!isSP2Atom(n, bk)) {
				  if (isHetero(getAtom(n)->na) && (getAtom(n)->nb == hVal[getAtom(n)->na])) testArom = true; else {
					  testArom = false;
					  break;
				  };
			  }
		  };
		  if (testArom) for (j = 0; j < cycleSize; j++) {
			  //if (getBond(bondList[j])->db == 7) getBond(bondList[j])->db = 3; else if (getBond(bondList[j])->db == 6) getBond(bondList[j])->db = 2;
			  if (getBond(bondList[j])->db == 0) getBond(bondList[j])->db = nn; else if (getBond(bondList[j])->db > nn) getBond(bondList[j])->db = nn;
		  }
	  };
  };


  for (i=0; i<nAtoms(); i++) {
	isSP3=isSP3Atom(i,bk);
	if (isSP3) getAtom(i)->hybridization=1; else {
	  isSP2=isSP2Atom(i,bk);
	  if (isSP2) getAtom(i)->hybridization=2; else {
		isSP=isSPAtom(i,bk);
		if (isSP) getAtom(i)->hybridization=3;
	  };
	};
	if (isAromaticAtom(i,bk)) getAtom(i)->aromatic=1; else getAtom(i)->aromatic=2;
  };

  //free(bk);
};

int  TSimpleMolecule::topologicalRadius(int atomIndex) {
	int result = 0;
	int i,j,n;
	bool test = true;

	if ((atomIndex < 0) || (atomIndex >= nAtoms())) return result;
	defineAtomConn();
	for (i = 1; i < nAtoms(); i++) getAtom(i)->enumerator = -1;
	getAtom(atomIndex)->enumerator = 0;
	while (test) {
		test = false;
		for (i = 0; i < nAtoms(); i++) if (getAtom(i)->enumerator == result) for (j=0; j<getAtom(i)->nb; j++) {
			n = getAtom(i)->ac[j];
			if (getAtom(n)->enumerator == -1) {
				getAtom(n)->enumerator = result + 1;
				test = true;
			};
		};
		if (test) result++;
	}

	return result;
};

int  TSimpleMolecule::topoDistance(int an1, int an2) const {
	int result = -1;
	int radius = 0;
	int i, j, k, n;
	std::vector<int> prevList, generatedList, usedList;
	

	prevList.push_back(an1);
	usedList.push_back(an1);
	while (prevList.size() > 0) {
		radius++;
		for (i = 0; i < prevList.size(); i++) {
			n = prevList[i];
			for (j = 0; j < getAtom(n)->nb; j++) {
				k = getAtom(n)->ac[j];
				if (k == an2) {
					result = radius;
					return result;
				};
				if (indexOfInteger(k, usedList) < 0) {
					usedList.push_back(k);
					generatedList.push_back(k);
				}
			}
		}
		prevList.clear();
		for (i = 0; i < generatedList.size(); i++) prevList.push_back(generatedList[i]);
		generatedList.clear();
	}
	return result;
}


void TSimpleMolecule::extractZeroBonds(int bondIndex, TSimpleMolecule & sMol, std::vector<adjustedlist> * bk) {
	//it is assumed that defineAtomConnectin was executed early
	std::vector<int> atomEnumerator;
	std::vector<int> tempList, currentSphereList;
	TSingleAtom * sa;
	TSingleBond * sb;
	bool test;
	int i, j, k, n, m;
	int n1, n2;
	adjustedlist bl;

	sMol.clear();
	sb = this->getBond(bondIndex);
	if (sb->tb != 0) return;
	sa = this->getAtom(sb->at[0])->clone();
	sMol.addAtom(sa);
	atomEnumerator.push_back(sb->at[0]);
	currentSphereList.push_back(sb->at[0]);
	sa = this->getAtom(sb->at[1])->clone();
	sMol.addAtom(sa);
	atomEnumerator.push_back(sb->at[1]);
	currentSphereList.push_back(sb->at[1]);
	while (currentSphereList.size() > 0) {
		test = false;
		tempList.clear();
		for (i = 0; i < currentSphereList.size(); i++) {
			bl = (*bk)[currentSphereList[i]];
			for (j = 0; j < bl.nb; j++) {
				sb = this->getBond(bl.adjusted[j]);
				if (sb->tb == 0) {
					n = indexOfInteger(sb->at[0], atomEnumerator);
					if (n < 0) {
						n = sb->at[0];
						sa = this->getAtom(n)->clone();
						sMol.addAtom(sa);
						atomEnumerator.push_back(n);
						tempList.push_back(n);
					};
					n = indexOfInteger(sb->at[1], atomEnumerator);
					if (n < 0){
						n = sb->at[1];
						sa = this->getAtom(n)->clone();
						sMol.addAtom(sa);
						atomEnumerator.push_back(n);
						tempList.push_back(n);
					};
				} else {  //not to add to tempList
					n = indexOfInteger(sb->at[0], atomEnumerator);
					if (n < 0) {
						n = sb->at[0];
						sa = this->getAtom(n)->clone();
						sMol.addAtom(sa);
						atomEnumerator.push_back(n);
					};
					n = indexOfInteger(sb->at[1], atomEnumerator);
					if (n < 0) {
						n = sb->at[1];
						sa = this->getAtom(n)->clone();
						sMol.addAtom(sa);
						atomEnumerator.push_back(n);
					};
				}
			}
		};
		currentSphereList.clear();
		for (i = 0; i < tempList.size(); i++) currentSphereList.push_back(tempList[i]);
	};

	for (i = 0; i < this->nBonds(); i++) {
		n1 = indexOfInteger(this->getBond(i)->at[0], atomEnumerator);
		n2 = indexOfInteger(this->getBond(i)->at[1], atomEnumerator);
		if ((n1 >= 0) && (n2 >= 0)) {
			sb = this->getBond(i)->clone();
			sb->at[0] = n1;
			sb->at[1] = n2;
			//sb->enumerator = i + 1;  
			if (sb->tb == 0) sb->enumerator = i+1;
			sMol.addBond(sb);
		}
	};
	sMol.defineAtomConn();
	sMol.allAboutCycles();
};


void TSimpleMolecule::extractSphericalEnviroment(int atomIndex, int topoDistance, TSimpleMolecule & sMol) {
	//it is assumed that defineAtomConnectin was executed early
	std::vector<int> atomEnumerator;
	std::vector<int> tempList, currentSphereList;
	TSingleAtom * sa;
	TSingleBond * sb;
	bool test;
	int i, j, k, n, m;
	int n1, n2;

	sMol.clear();
	sa = this->getAtom(atomIndex)->clone();
	sMol.addAtom(sa);
	atomEnumerator.push_back(atomIndex);
	currentSphereList.push_back(atomIndex);
	for (i = 0; i < topoDistance; i++) {
		tempList.clear();
		for (j = 0; j < currentSphereList.size(); j++) {
			n = currentSphereList[j];
			for (k = 0; k < this->getAtom(n)->nb; k++) {
				m = this->getAtom(n)->ac[k];
				if (indexOfInteger(m, atomEnumerator) < 0) {
					sa = this->getAtom(m)->clone();
					sMol.addAtom(sa);
					atomEnumerator.push_back(m);
					tempList.push_back(m);
				};
			};
		};
		//currentSphereList.assign(tempList);
		currentSphereList.clear();
		for (j = 0; j < tempList.size(); j++) currentSphereList.push_back(tempList[j]);
	};
	for (i = 0; i < this->nBonds(); i++) {
		n1 = indexOfInteger(this->getBond(i)->at[0],atomEnumerator);
		n2 = indexOfInteger(this->getBond(i)->at[1],atomEnumerator);
		if ((n1 >= 0) && (n2 >= 0)) {
			sb = this->getBond(i)->clone();
			sb->at[0] = n1;
			sb->at[1] = n2;
			if (sb->tb == 0) sb->tb = 5;
			sMol.addBond(sb);
		}
	};
	sMol.defineAtomConn();
	sMol.allAboutCycles();
};


//----------------------------------------------------------------------------------------------------


std::string getAtomSymbol(TSimpleMolecule & sm, int atAtom, int atEx, int priority, std::string ndData) {
 //!  nepravil'no
  //For given atom AtAtom returns atom symbol with given prioritate. AtEx-do not take part into calculations
  std::string result="";
  std::vector<std::string> collectedSymbols(15);
  //string collectedSymbols [] = new string [15];
  int nPrior=0;
  int i,j,n;
  bool test;

  //result=aSymb[sm.getAtom(atAtom)->na];
  //result="at="+intToStr(atAtom)+" nt="+intToStr(sm.nAtoms())+" na="+intToStr(sm.getAtom(atAtom)->na);
 
  for (i=0; i<sm.getAtom(atAtom)->nb; i++) {
    n=sm.getAtom(atAtom)->ac[i];
    if (n != atEx) {
      nPrior++;
      test=false;
      if (sm.getAtom(n)->anum.length() > 0) test=true;
      if (test) collectedSymbols[nPrior-1]=sm.getAtom(n)->anum; else collectedSymbols[nPrior-1]=aSymb[sm.getAtom(n)->na];
    };
  };
  n=sm.getNH(atAtom);
  for (i=0; i<n; i++) {
    nPrior++;
    collectedSymbols[nPrior-1]="H";
  };
  if (nPrior < 4) for (i=0; i<4; i++) {
    nPrior++;
    collectedSymbols[nPrior-1]=ndData;//"zz";  //instead of 0 - less prioritate....
    if (nPrior == 4) break;
  };
  for (i=0; i<(nPrior-1); i++) for (j=i+1; j<nPrior; j++) if (compareStringsNumbers(collectedSymbols[i],collectedSymbols[j])>0) {
    result=collectedSymbols[i];
    collectedSymbols[i]=collectedSymbols[j];
    collectedSymbols[j]=result;
  };
  if ((priority > 0) && (priority <= nPrior)) result=collectedSymbols[priority-1];
  return result;
};


std::string getAtomSymbol(TSimpleMolecule & sm, int atAtom) {
  bool test=false;
  string result="";

  if (sm.getAtom(atAtom)->anum.length() > 0) test=true;
  if (test) result=sm.getAtom(atAtom)->anum; else result=aSymb[sm.getAtom(atAtom)->na];
  return result;
};

double getAtomMass(TSimpleMolecule & sm, int atAtom) {
	double		result = 0.0;

	result = aMass[sm.getAtom(atAtom)->na];
	return result;
};

int getAtomValency(TSimpleMolecule & sm, int atAtom) {
	int		result = 0;

	result = hVal[sm.getAtom(atAtom)->na];
	return result;
};

//-------------------------------------------------------------------------------------------
#define MAX_DIM_SYMM 1000 
#define R_UNDEF -12345678

int determineSymmetry(TSimpleMolecule & sm, TSimpleMolecule & em) {
  int result=0;
  int  diffAtomsCodes[MAX_DIM_SYMM];
  int  diffAtomsNo[MAX_DIM_SYMM];
  int  diffAtomIndex[MAX_DIM_SYMM];
  double x,y,x2,y2,xy,xMin,yMin,xMax,yMax,r1,r2,z,a,b,s;
  std::vector<int>eq;
  int i,j,n,nDiffAtoms,nSym,nLyn,nLyn_2,code;

  nDiffAtoms=0;
  code=0;
  for (i=0; i<MAX_DIM_SYMM; i++) {
    diffAtomsCodes[i]=-1;
    diffAtomsNo[i]=-1;
    diffAtomIndex[i]=-1;
  };
  sm.defineAtomConn();
  sm.allAboutCycles();
  em.moleculeCopy(sm);
  em.makeEquivalentList(eq,false);
  for (i=0; i<eq.size(); i++) if (eq[i] != 0) {
    diffAtomsCodes[nDiffAtoms]=eq[i];
    diffAtomsNo[nDiffAtoms]=1;
    diffAtomIndex[nDiffAtoms]=i;
	em.getAtom(i)->enumerator=eq[i];
    eq[i]=0;
	for (j=i+1; j<eq.size(); j++) if (eq[j] ==  diffAtomsCodes[nDiffAtoms]) {
      diffAtomsNo[nDiffAtoms]=diffAtomsNo[nDiffAtoms]+1;
	  em.getAtom(j)->enumerator=eq[j];
      eq[j]=0;
	};
    nDiffAtoms++;
  };

  x=0; y=0; x2=0; y2=0; xy=0; nSym=1000; nLyn=0; xMin=R_UNDEF; yMin=R_UNDEF; xMax=R_UNDEF; yMax=R_UNDEF;

  for (i=0; i<nDiffAtoms; i++) {
	if (diffAtomsNo[i] == 1) {
      nLyn++;
      n=diffAtomIndex[i];
	  r1=em.getAtom(n)->rx;
	  r2=em.getAtom(n)->ry;
      x=x+r1;
      y=y+r2;
      x2=x2+r1*r1;
      y2=y2+r2*r2;
      xy=xy+r1*r2;
	};
    if ((diffAtomsNo[i] >= 2) && (diffAtomsNo[i] < nSym)) nSym=diffAtomsNo[i];
	if (diffAtomsNo[i] == 4) {  //TD symmetry
      if (eq.size() == 5) nSym=12;
	};
  };
  if (nLyn <= 2) {
    if (eq.size() == 2) {
	  if ((nSym == 2) && (em.getAtom(0)->na == 6)) {
        nSym=1;
		if (em.getBond(0)->tb == 1) nSym=6; else   //Ethane
		if (em.getBond(0)->tb == 2) nSym=4;       //Ethene
	  } else if (nSym == 1) {
		if ((em.getAtom(0)->na == 6) || (em.getAtom(1)->na == 6)) nSym=3; //methyl derivatives
	  };
	};
    if (nSym < 1000) result=nSym;
  } else {
    //Method of less-squares
    z=nLyn*x2-x*x;
	if (abs(z) > 0.1) {
      b=(nLyn*xy-x*y)/z;
      a=(y*x2-xy*x)/z;
      code=1;
	} else {
      //Almost vertical line
      z=nLyn*y2-y*y;
      b=(nLyn*xy-x*y)/z;
      a=(x*y2-xy*y)/z;
      code=2;
	};
    s=0;
	for (i=0; i<nDiffAtoms; i++) if (diffAtomsNo[i] == 1) {
      n=diffAtomIndex[i];
	  r1=em.getAtom(n)->rx;
	  r2=em.getAtom(n)->ry;
	  if (code == 1) {
        y=a+b*r1;
        s=s+(y-r2)*(y-r2);
	  } else {
        x=a+b*r2;
        s=s+(x-r1)*(x-r1);
	  };
	};
    s=sqrt(s/(nLyn-2));
	if (s < 0.1) {
      if (nSym < 1000) result=nSym;
    };
  };
  if (eq.size() == 1) {
	if (sm.getAtom(0)->na == 6) result=12;
	if (sm.getAtom(0)->na == 7) result=3;
  };
  if (result == 0) result=1;
  return result;
}

void analyzeRotors(TSimpleMolecule & sm, std::vector<std::string> & explanation, std::vector<int> & symnumber){
  std::vector<int> eq;
  TSimpleMolecule em;
  std::vector<adjustedlist> bk(sm.listarSize());
  std::vector<int> tempList;
  bool test;
  int n,i,j,k,nH,an1,bn1;
  std::string s;
  int occurences[16];
  int atomN[16];
  int formulaData[NELEMMAX];
  int w;
  std::vector<int>aList;

  explanation.clear();
  symnumber.clear();

  n=determineSymmetry(sm,em);
  explanation.push_back("Molecular symmetry, symmetry number: "+intToStr(n));
  symnumber.push_back(n);
  em.defineBondConn(bk);
  //Rotors analyzing
  for (i=0; i<em.nAtoms(); i++) {
	nH=em.getNH(i);
	if (em.getAtom(i)->nb == 1) {
	  k=bk[i].adjusted[0];
	  if ((nH > 1) && ((em.getBond(k)->tb == 1) || (em.getBond(k)->tb >=9 ))) {
        s=getAtomSymbol(em,i)+"H"+intToStr(nH)+" rotor, symmetry number: "+intToStr(nH);
		explanation.push_back(s);
        symnumber.push_back(nH);
	  };
	} else if ((nH == 0) && (em.getAtom(i)->nb >= 3)) {
        //explanation.AddObject('CH3: 3',pointer(3));
      tempList.clear();
	  for (j=0; j<16; j++){ 
        occurences[j]=0;
        atomN[j]=-1;
	  };
	  for (j=0; j<em.getAtom(i)->nb; j++) {
        n=em.getAtom(i)->ac[j];
		n=em.getAtom(n)->enumerator;
		k=-1;
		for (int kk=0; kk<tempList.size(); kk++) if (tempList[kk] == n) {
          k=kk;
		  break;
		};
		if (k < 0) {
          tempList.push_back(n);
		  k=tempList.size()-1;
		};
        occurences[k]=occurences[k]+1;
		atomN[k]=em.getAtom(i)->ac[j];
	  };
      test=false; an1=-1; bn1=-1;
	  
	  if (tempList.size() == 2) if ((occurences[0] == 1) || (occurences[1] == 1)) {
        //Possyble symmetric
		if (occurences[0] == 1) {
		  an1=atomN[0];
		} else {
          an1=atomN[1];
		};
        //an1 - single atom. Checking that:
        //1) to an1 connected at least one more atom and
        //2) not double bond
		if ((em.getAtom(an1)->nb > 1) || (em.getNH(an1) > 0)) for (j=0; j<bk[an1].nb; j++) {
		  n=bk[an1].adjusted[j];
		  if (((em.getBond(n)->at[0] == i) && (em.getBond(n)->at[1] == an1)) || ((em.getBond(n)->at[0] == an1) && (em.getBond(n)->at[1] == i))) {
			test=((em.getBond(n)->tb == 1) || (em.getBond(n)->tb >= 9)) && (em.getBond(n)->db == 0);
            if (test) bn1=n;
		  };
          if (test) break;
	    };
	  };
	  
	  if (test) {
        //Fine! Found! determine molecular formula
		aList.resize(em.nAtoms());
		em.makeFragment(w,aList,i,an1);
        for (j=0; j<NELEMMAX; j++) formulaData[j]=0;
		for (j=0; j<w; j++) {
          n=aList[j];
		  k=em.getAtom(n)->na;
          formulaData[k-1]=formulaData[k-1]+1;
          k=em.getNH(n);
          formulaData[0]=formulaData[0]+k;
		};
		
        s="";
		if (formulaData[5] > 0) {
          s=s+aSymb[5+1];
          if (formulaData[5] > 1) s=s+intToStr(formulaData[5]);
		};
		if (formulaData[0] > 0) {
          s=s+aSymb[0+1];
          if (formulaData[0] > 1) s=s+intToStr(formulaData[0]);
		};
		for (j=1; j<100; j++) if ((j != 5) && (formulaData[j] > 0)) {
          s=s+aSymb[j+1];
          if (formulaData[j] > 1) s=s+intToStr(formulaData[j]);
		};
		n=em.getAtom(i)->nb-1;
        s=s+" rotor, symmetry number: "+intToStr(n);
        explanation.push_back(s);
		symnumber.push_back(n);
		
	  };
	  
	};
  };
};


bool TSimpleMolecule::isAromatic(const std::vector<int> & bondList) const {
	bool result = false;
	int i,n;
	for (i = 0; i < bondList.size(); i++) {
		n = bondList[i];
		result = ((getBond(n)->db == 2) || (getBond(n)->db == 3));
		if (!result) break;
	}
	return result;
}

void TSimpleMolecule::bondListToAtomList(const std::vector<int> & bondList, std::vector<int> & atomList) const {
  std::vector<int>tempList;
	int i, n;
    atomList.clear();
	if (bondList.size() == 0) return;
	atomList.clear();
	atomList.reserve(bondList.size() + 2);
	tempList.reserve(2 * bondList.size());
	for (i = 0; i < bondList.size(); i++) { //if (bondList[i]<nBonds()) {
		n = getBond(bondList[i])->at[0];
		tempList.push_back(n);
		n = getBond(bondList[i])->at[1];
		tempList.push_back(n);
	};
	quickSortIntegers(tempList);
	atomList.push_back(tempList[0]);
	for (i=1; i<tempList.size(); i++) if (tempList[i - 1]!=tempList[i]) atomList.push_back(tempList[i]);
}

int TSimpleMolecule::getOppositeLabel(int atomNo) const {
	//is used for zvezda polymer search
		//must check if fragments are connected...
	int j;
	int result = -1;
	for (j = 0; j < nAtoms(); j++) if ((j != atomNo) && (getAtom(j)->na == getAtom(atomNo)->na)) {
		result = j;
		break;
	};
	return result;
};


int TSimpleMolecule::getOppositeAttachedAtom(int atomNo) const {
  //must check if fragments are connected...
	int n;
    int result = -1;
	n = getOppositeLabel(atomNo);
	if (n >= 0) result = getAtom(n)->ac[0];
	return result;
};

int TSimpleMolecule::bondNumber(int an1, int an2) const {
	int i, at1, at2;
	int result = -1;
	for (i = 0; i < nBonds(); i++) {
		at1 = getBond(i)->at[0];
		at2 = getBond(i)->at[1];
		if (((at1 == an1) && (at2 == an2)) || ((at1 == an2) && (at2 == an1))) {
			result = i;
			break;
		};
	};
	return result;
};

bool TSimpleMolecule::atomInRing(int atomNo, const std::vector<int> & bondList) const {
	int i, n;
	bool result=false;

	for (i = 0; i < bondList.size(); i++) {
		n = bondList[i];
		if ((getBond(n)->at[0] == atomNo) || (getBond(n)->at[1] == atomNo)) {
			result = true;
			break;
		};
	};
	return result;
};

bool TSimpleMolecule::hasDoubleBond(int atomNo, const  std::vector<adjustedlist> & bk) const {
	int i, n;
	bool result = false;

	for (i = 0; i < bk[atomNo].nb; i++) {
		n = bk[atomNo].adjusted[i];
		if (getBond(n)->tb == 2) {
			result = true;
			break;
		};
	};
	return result;
}


void TSimpleMolecule::cyclesCalculate(int & nTotalCycles, int & nAromFive, int & nAromSix, int & nCondenced, std::vector<std::vector<int> *> * ringList) const {
	//It is assumed, that AllAboutCycles routine has been called early
	//Ncondenced = -1 - does not performs polycycles calculations !!!!
	std::vector<std::vector<int> *> ringData;
	std::vector<int> * tempList;
	std::vector<int> * tempList1;
	std::vector<int> * tempList2;
	std::vector<adjustedlist> bk(listarSize());
	std::vector<int> bondList(listarSize());
	int rs;  //word
	int i, j, k, m;
	bool test, testAdded;

	nTotalCycles= 0;
	nAromFive = 0;
	nAromSix = 0;
	//!!!Do NOT set initial value of NCondenced because of - 1 will be analized later!
	defineBondConn(bk);
	tempList = new (std::vector<int>);
	for (i = 0; i < nBonds(); i++) {
		vaweBond(i, bk, rs, bondList, 0);
		if (rs > 0) {
			tempList->clear();
			for (j = 0; j < rs; j++) {
				m = bondList[j];
				tempList->push_back(m);
			};
			quickSortIntegers(*tempList);
			if (!listPresent(ringData, *tempList)) {
				ringData.push_back(tempList);
				tempList = new (std::vector<int>);
			};
		};
	};
	nTotalCycles = ringData.size();
	if (ringData.size() > 0) for (i = 0; i < ringData.size(); i++) if (isAromatic(*ringData[i])) {
		if (ringData[i]->size() == 5) nAromFive++; else nAromSix++;
	};

	//Condenced No.of cycles
	if (nCondenced != -1) {
		nCondenced = 0;
		test = true;

		while (test) {
			test = false;
			if (ringData.size() > 1) for (i = 0; i <= ringData.size() - 2; i++) {
				tempList1 = ringData[i];
				for (j = i + 1; j < ringData.size(); j++) {
					tempList2 = ringData[j];
					for (k = 0; k < tempList2->size(); k++) if (bondPresent(*tempList1, (*tempList2)[k])) {
						test = true;
						nCondenced++;
						break;
					};
					if (test) {
						ringData[i] = NULL;
						ringData[j] = NULL;
						for (k = 0; k<tempList2->size(); k++) tempList1->push_back((*tempList2)[k]);
						delete (tempList2);
						tempList2 = NULL;
						packRingData(ringData);
					};
					if (test) break;
				};
				if (test) break;
			};

			if (test) {
				testAdded = true;
				while (testAdded) {
					testAdded = false;
					for (i = ringData.size() - 1; i >= 0; i--) {
						tempList2 = ringData[i];
						for (j = 0; j < tempList2->size(); j++)if (bondPresent(*tempList1, (*tempList2)[j])) {
							testAdded = true;
							break;
						};
						if (testAdded) {
							ringData[i] = NULL;
							packRingData(ringData);
						};
						if (testAdded) break;
					};
					if (testAdded) {
						for (i = 0; i<tempList2->size(); i++) tempList1->push_back((*tempList2)[i]);
						delete (tempList2);
					};
				}; // until not TestAdded;
				delete (tempList1);
			};
		}; //until not Test;
		
	};  
	if ((ringList != NULL) && (ringData.size()>0)) {
		ringList->clear();
		ringList->reserve(ringData.size());
		for (i=0; i<ringData.size(); i++) ringList->push_back(ringData[i]);
		ringData.clear();
	};
	for (i = 0; i < ringData.size(); i++) delete (ringData[i]);
};

void TSimpleMolecule::removeFarPendant(const std::vector<std::vector<int> *> & ringList, const std::vector<int> & mainChain, int backboneDistance) {
	std::vector<int> addedAtoms, accumList, tempList, atomList;
	std::vector<int> * bondList;
	int i,j,k,n,nn, deepIndex;
	bool closeToBackbone;
	TSimpleMolecule tempMol;
	bool ringAdded;

	for (i = 0; i < nAtoms(); i++) if (findQuick(mainChain, mainChain.size(), i) < 0) {
		closeToBackbone = false;
		accumList.clear();
		accumList.push_back(i);
		tempList.clear();
		for (deepIndex = 0; deepIndex < backboneDistance; deepIndex++) {
			for (j = 0; j < accumList.size(); j++) {
				n = accumList[j];
				for (k = 0; k < getAtom(n)->nb; k++) {
					nn = getAtom(n)->ac[k];
					if (findQuick(mainChain, mainChain.size(), nn) >= 0) {
						closeToBackbone = true;
						break;
					} else tempList.push_back(nn);
				}
				if (closeToBackbone) break;
			}
			if (closeToBackbone) break;
			accumList.clear();
			for (j = 0; j < tempList.size(); j++) accumList.push_back(tempList[j]);
			tempList.clear();
		}

		if (closeToBackbone) addedAtoms.push_back(i);
	};
	//rings processing
	for (i = 0; i < ringList.size(); i++) {
		bondList = ringList[i];
		bondListToAtomList((*bondList), atomList);
		ringAdded=false;
		for (j = 0; j < addedAtoms.size(); j++) if (findQuick(atomList, atomList.size(), addedAtoms[j]) >= 0) {
			ringAdded = true;
			break;
		}
		if (ringAdded) for (j = 0; j < atomList.size(); j++) addedAtoms.push_back(atomList[j]);
	}
	//backbone list addition
	for (i = 0; i < mainChain.size(); i++) addedAtoms.push_back(mainChain[i]);
	//remove duplicated definitions
	quickSortIntegers(addedAtoms);
	atomList.clear();
	atomList.reserve(addedAtoms.size());
	atomList.push_back(addedAtoms[0]);
	for (i = 1; i < addedAtoms.size(); i++) if (addedAtoms[i - 1] != addedAtoms[i]) atomList.push_back(addedAtoms[i]);
	tempMol.moleculeCopy(*this);
	tempMol.defineAtomConn();
	tempMol.extractFragmentFast(atomList, *this);
	defineAtomConn();
};

/*
std::string TSimpleMolecule::calculateInChIKey() {
	return convertMolToInChIKey(*this);
}

std::string TSimpleMolecule::calculateInChIKeyEx() const {
	return convertMolfileToInchIKey(*this);
}
*/

void TSimpleMolecule::makeStandardPolymer(const std::vector<int> & mainChain) {
	TSimpleMolecule smTemp;
	TSimpleMolecule * fragMol;
	int an1, an2, i, n, an1a, an2a, backboneLength, an3, an4, bn;
	std::vector<int> eqList, chainList, atomEnumerator;
	bool test;
	int unitLength, nFrag, topoDistance;
	std::vector<adjustedlist> bk;

	bk.resize(listarSize());
	defineBondConn(bk);
	
	defineAtomConn();
	smTemp.moleculeCopy(*this);

	an1 = -1; an2 = -1; an1a = -1; an2a = -1;
	for (i = 0; i < smTemp.nAtoms(); i++) if (smTemp.getAtom(i)->na == ID_ZVEZDA) {
		if (an1 < 0) an1 = i; else {
			an2 = i;
			break;
		};
	};
	for (i = 0; i < smTemp.nBonds(); i++) if (smTemp.getBond(i)->db < 2) chainList.push_back(i);
	quickSortIntegers(chainList);
	//make cyclic
	if ((an1 >= 0) && (an2 >= 0) && (smTemp.getAtom(an1)->nb == 1) && (smTemp.getAtom(an2)->nb == 1)) {
		//make cyclic
		n = smTemp.topoDistance(an1, an2);
		an1a = smTemp.getAtom(an1)->ac[0];
		an2a = smTemp.getAtom(an2)->ac[0];
		if (n == 2) {
			clear();
			return;
		} else if (n >= 4) {
			smTemp.addBond(1, an1a, an2a);
			//smTemp.defineAtomConn();
			for (i = smTemp.nBonds() - 1; i >= 0; i--) if ((smTemp.getBond(i)->at[0] == an2) || (smTemp.getBond(i)->at[1] == an2)) {
				smTemp.deleteBond(i);
				break;
			};
			for (i = smTemp.nBonds() - 1; i >= 0; i--) if ((smTemp.getBond(i)->at[0] == an1) || (smTemp.getBond(i)->at[1] == an1)) {
				smTemp.deleteBond(i);
				break;
			};
		};
		smTemp.defineAtomConn();
		smTemp.allAboutCycles();
		smTemp.makeEquivalentList(eqList, false);
		if (n == 3) {
			if (eqList[an1a] == eqList[an2a]) {
				unitLength = 1;
				nFrag = 2;
				n = -1;
				for (i = nBonds() - 1; i >= 0; i--) if (((getBond(i)->at[0] == an1a) && (getBond(i)->at[1] == an2a)) || ((getBond(i)->at[0] == an2a) && (getBond(i)->at[1] == an1a))) {
					deleteBond(i); 
					break;
				}
				addAtom(ID_ZVEZDA, 0, getAtom(an2a)->rx, getAtom(an2a)->ry);
				addBond(1, an1a, nAtoms() - 1);
				defineAtomConn();
				fragMol = extractFragment(an1a, NULL);
				moleculeCopy(*fragMol);
				defineAtomConn();
				allAboutCycles();
				delete fragMol;
				return;
			} else {
				clear();
				return;
			}
		} else {
			unitLength = 1;
			//check if single atom
			n = -1;
			for (i = 0; i < eqList.size(); i++) {
				test = (i == an1) || (i == an2);
				if (!test) {
					if (n == -1) n = eqList[i]; else unitLength = 0;
				}
				if (unitLength == 0) break;
			};
			if (unitLength == 0) {
				//copy original molecule back to smTemp
				smTemp.clear();
				smTemp.moleculeCopy(*this);
				smTemp.atomDistanceList(an1a, an1, chainList);
				//an1 will be start it has distance=0 
				//an2-end of chain it has maximal topo distance 
				topoDistance = chainList[an2];
				if (topoDistance < 0) {  //not connected atoms
					clear();
					return;
				}
				nFrag = 0;
				for (unitLength = 2; unitLength <= (topoDistance / 2); unitLength++) if ((topoDistance % unitLength) == 0) {
					nFrag = topoDistance / unitLength;
					for (i = 1; i < nFrag; i++) {
						test = fragmentsIdentical(0, i*unitLength, unitLength, chainList, eqList);
						if (!test) {
							nFrag = 0;
							break;
						}
					};
					if (nFrag > 0) break;
				}
				if (nFrag > 0) unitLength = topoDistance / nFrag; else unitLength = 0;
			}
		};
		if (unitLength > 0) {
			an3 = -1; an4 = -1; bn = -1;
			for (i = 0; i < chainList.size(); i++) if ((chainList[i] == (unitLength - 1)) && (findQuick(mainChain, mainChain.size(), i) >= 0)) {
				if (an3 == -1) an3 = i; else an3 = -2;
			};
			if (an3 >= 0) {
				for (i = 0; i < bk[an3].nb; i++) {
					bn = bk[an3].adjusted[i];
					an4 = getBond(bn)->at[0];
					if (an4 == an3) an4 = getBond(bn)->at[1];
					if ((getBond(bn)->db < 2) && (findQuick(mainChain, mainChain.size(), an4) >= 0) && (chainList[an4] == unitLength)) break; else {
						bn = -1;
						an4 = -1;
					}
				}
			};
			if (an4 >= 0) {
				//move asteric to an3 atom
				n = bk[an2].adjusted[0];  //far asteric-single bond
				if (getBond(n)->at[0] == an2) getBond(n)->at[1] = an3; else getBond(n)->at[0] = an3;
				deleteBond(bn);
				defineAtomConn();
				fragMol = extractFragment(an3, NULL);
				//fragMol->saveMolfile("c:/temp/error.mol");
				moleculeCopy(*fragMol);
				defineAtomConn();
				allAboutCycles();
				delete fragMol;
			}
		}
		else clear();
	}
	else {
		clear();
		return;
	};



	return;
};




bool fragmentsIdentical(int fromSphere1, int fromSphere2, int sphereCount, const std::vector<int> & sphereList, const std::vector<int> & eqList) {
	//for polymers backbone chain. Check if 
	bool result = true;
	int i, j, levelNo1, levelNo2;
	std::vector<int> sphere1EQ, sphere2EQ;

	if (sphereList.size() != eqList.size()) return false;
	for (i = 0; i < sphereCount; i++) {
		levelNo1 = fromSphere1 + i;
		levelNo2 = fromSphere2 + i;
		sphere1EQ.clear();
		sphere2EQ.clear() ;
		for (j = 0; j < sphereList.size(); j++) {
			if (sphereList[j] == levelNo1) sphere1EQ.push_back(eqList[j]);
			if (sphereList[j] == levelNo2) sphere2EQ.push_back(eqList[j]);
		}
		if (sphere1EQ.size() == sphere2EQ.size()) {
			quickSortIntegers(sphere1EQ);
			quickSortIntegers(sphere2EQ);
			for (j = 0; j < sphere1EQ.size(); j++) if (sphere1EQ[j] != sphere2EQ[j]) {
				result = false;
				break;
			}
		} else result = false;

		if (!result) break;
	}

	return result;
};


} // namespace MolStruct



