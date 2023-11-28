
#define _USE_MATH_DEFINES

#include <vector>
#include <string>
#include <cmath>
#include <cstdio>

#include "edited_molecule.h"
#include "scalar.h"
//#include "maccs.h"

using namespace std;

namespace MolStruct {

//#define DEFAULTBONDLENGTH 1.44
#define DEFAULTBONDLENGTH		100

//Refacor it !!!!!!!
/*
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
*/

int TEditedMolecule::addFragment(TSimpleMolecule & eMolecule, int naDEF, int cha, int chb,
	int chb1, std::vector<int>& list, double xOldCenter, double yOldCenter, double xNewCenter,
	double yNewCenter, double scale, double cFi, double sFi, int buttonStatus, bool clearEnumerator) {

		//    The procedure is used by TEMPLATE and MAKEPOLI procedure to add large frag-
		//     ments to structure. The connection between structure and fragment may be
		//     achived both through atoms (CHA>0) and through bonds (CHB>0).
		//                              Variables description:
		//     NAtoms-number of atoms in structure
		//     NBonds-number of bonds in structure
		//     ATOM-atom's attributes in structure
		//     BOND-bond's attributes in structure
		//     CONN-some invariants of bond-connection matrix
		//     NAtoms1, NBonds1, ATOM1, BOND1, CONN1-the same for fragments (these varia-
		//       bles can contain data for more, than one fragments)
		//     LIST-atom's number in arrays ATOM1,CONN1, which should be added to structure
		//     NaDef-number of principal components in LIST array
		//     CHA-atom's number in array ATOM (structure), which was selected for addition
		//     CHB-bond's number in array BOND (structure), which was selected for addition
		//     CHB1-bond's number in array BOND1 (fragment), which was selected for addition
		//     XOldCenter,YOldCenter-old coordinates of the fragment on screen (center rotation)
		//     XNewCenter,YNewCenter-new coordinates of the fragment on screen
		//     Scale-scale variables to resize fragment
		//     CFi,SFi-cosine and sine of the angle of rotation of the fragment relative
		//      old coordinate of ths type of a fragment addition (active only for addition
		//       through atoms). If e fragment on screen
		//     ButtonStatus-regulateButtonStatus=1 then fragment will be added to structure
		//       (new bond arise). If ButtonStatus=2 then fragment will be freezed in
		//       structure (selected new atom for structure will be deleted, no new bond
		//       will be arised)}

		//    Returns 0 if addition is OK,  otherwise:
		//    =1 - Maxima number of atoms is reached
		//    =2 - Maximal number of bonds is reached

		int i,j1,j2,na,k1,k2,j,nb,n1,n2;
		int l1=0;
		int l2=0;
		bool test;
		double r1,r2,r3,r4,emBLength;
		std::vector<int> bTested(nBonds() + 4 * naDEF);  //NBONDSMAX);
		std::vector<int> bList(nBonds() + 4 * naDEF); //NBONDSMAX);
		std::vector<int> aList(nAtoms() + 4 * naDEF); //NBONDSMAX);
		TSingleAtom * sa=NULL;
		TSingleBond * sb=NULL;
		int result=0;



		//return result;
		//Analizing on boundaries}
		k2=1; k1=1;
		/*
		if ((nAtoms()+naDEF) >= NATOMSMAX) {
			result=1;
			return result;
		}
		*/
		//Analizing of bond's
		nb=0;
		emBLength=0;
		for (i=0; i<eMolecule.nAtoms(); i++) bList[i]=0;
		for (i=0; i<naDEF; i++) bList[list[i]]=1;
		for (i=0; i<eMolecule.nBonds(); i++) {
			n1=eMolecule.getBond(i)->at[0];
			n2=eMolecule.getBond(i)->at[1];
			if ((bList[n1] != 0) || (bList[n2] != 0)) {
				nb++;
				emBLength=emBLength+eMolecule.bondLength(i);
			}
		}
		if (nb > 0) emBLength=emBLength/nb;
        /*
		if ((nb+nBonds()) >= NBONDSMAX) {
			result=2;
			return result;
		}
		*/
		na=nAtoms();
		nb=nBonds();


		if (chb >= 0) {   //Bond was selected for addition

			l1=getBond(chb)->at[0];     //L1,L2,K1,K2-atoms numbers for corresponding bonds
			l2=getBond(chb)->at[1];     //in structure and in fragment
			k1=eMolecule.getBond(chb1)->at[0];
			k2=eMolecule.getBond(chb1)->at[1];
			j=0;
			for (i=0; i<naDEF; i++) {
				test=(list[i] == k1) || (list[i] == k2);
				if (! test) {
					list[j]=list[i];
					j++;
				};
			};    
			list[naDEF-2]=k1;
			list[naDEF-1]=k2;
			//It is necessary to make the pair of atoms last in array so they may be
			// deleted easy from final structure
		}
		//Addition from 13 May 2000
		if (((cha < 0) || (nBonds() == 0)) && (emBLength > 0)) {
			if (nBonds() == 0) r1=DEFAULTBONDLENGTH; else r1=averageBondLength();
			scale=r1/emBLength;
		}
		//End addition

		for (i=0; i<naDEF; i++) {
			//recalculation of atom's coordinates in fragment to new ones, compartible
			//with the structure coordinates. Some rescaling and rotation take place
			r1=eMolecule.getAtom(list[i])->rx-xOldCenter;
			r2=eMolecule.getAtom(list[i])->ry-yOldCenter;
			r3=scale*(r1*cFi-r2*sFi);
			r4=scale*(r1*sFi+r2*cFi);
			sa=eMolecule.getAtom(list[i])->clone();
			sa->rx=xNewCenter+r3;
			sa->ry=yNewCenter+r4;
			if (clearEnumerator) sa->enumerator=0;
			sa->fragIndex=0;
			addAtom(sa);
		}
		if (chb >= 0) {
			//specific for bond-to-bond addition block-determination of correspondence
			//of pair of atoms, which make bond, in fragment to the same pair of atoms in
			//structure. Determination is perfomed through internuclear distances}
			r1=getAtom(l2)->rx-getAtom(nAtoms()-1)->rx;
			r2=getAtom(l2)->ry-getAtom(nAtoms()-1)->ry;
			r3=r1*r1+r2*r2;
			r1=getAtom(l1)->rx-getAtom(nAtoms()-1)->rx;
			r2=getAtom(l1)->ry-getAtom(nAtoms()-1)->ry;
			r4=r1*r1+r2*r2;
			if (r3 < r4) {k2=nAtoms()-1; k1=nAtoms()-2;
			} else {k2=nAtoms()-2; k1=nAtoms()-1;}
		}

		for (i=0; i<eMolecule.nBonds(); i++) if ((chb < 0) || ((chb1 >= 0) && (i != chb1))) {
			//addition of bonds from fragment to structure
			test=false;
			j1=-1;
			do { //Search, if bond I from BOND1 is presented in the fragment
				j1++;
				if (eMolecule.getBond(i)->at[0] == list[j1]) test=true;
			}  while (! (test || (j1 == (naDEF-1))));
			if (test) {  //Bond is presented-addition
				j2=-1;
				do     //Search for second atom, which makes I-th bond in fragment
				j2++;
				while (eMolecule.getBond(i)->at[1] != list[j2]);
				sb=new TSingleBond();
				j1=na+j1;            //atoms, which make the bond, renumeration
				j2=na+j2;
				sb->tb=eMolecule.getBond(i)->tb;  //attributes copy
				sb->at[0]=j1;
				sb->at[1]=j2;
				if (clearEnumerator) sb->enumerator=0;
				addBond(sb);
			}
		}


		if (chb >= 0) { //bond-to-bond addition
			for (i=0; i<nBonds(); i++) { //atom renumeration
				if (getBond(i)->at[0] == k1) getBond(i)->at[0]=l1;
				if (getBond(i)->at[0] == k2) getBond(i)->at[0]=l2;
				if (getBond(i)->at[1] == k1) getBond(i)->at[1]=l1;
				if (getBond(i)->at[1] == k2) getBond(i)->at[1]=l2;
			}
			//Two atoms deleted...
			delete(fAtom[nAtoms()-1]);
			delete(fAtom[nAtoms()-2]);
			fAtom.resize(nAtoms()-2);


			//two atoms, which make the bond to be added in fragment,
			//are deleted. They occupy last positions in arrays ATOM and CONN due to
			//initial sorting
			/* Contains mistake-commented to the better time
			i=0;
			do {                 //Two or three bond connection is calculated
			j=na;
			if (j < nAtoms()) do {               //Search for concided atoms
			rr1=getAtom(i)->rx-getAtom(j)->rx;
			rr2=getAtom(i)->ry-getAtom(j)->ry;
			lR=sqrt(rr1*rr1+rr2*rr2);
			test=((i != j) & ((lR) <= 3));
			if (test) { //Concided pair of atoms found
			for (k=0; k<nBonds(); k++) {
			//Re-assignment of bonds
			if (getBond(k)->at[0] == j) getBond(k)->at[0]=i;
			if (getBond(k)->at[1] == j) getBond(k)->at[1]=i;
			}
			deleteAtom(j);                          //Deletion of atom
			k=0;
			do {      //Search for bonds, which have identical pair of atom or which
			l=k;
			do {
			test=(getBond(k)->at[0] == getBond(l)->at[0])
			&& (getBond(k)->at[1] == getBond(l)->at[1]);
			test=test || (getBond(k)->at[0] == getBond(l)->at[1])
			&& (getBond(k)->at[1] == getBond(l)->at[0]);
			test=test || (getBond(l)->at[0] == getBond(l)->at[1]);
			if (test) {         //Bad bonds was found-deletion
			deleteBond(l);
			l--;
			}
			l++;
			} while (l < nBonds());
			k++;
			} while (k < (nBonds()-1));
			j--;
			}
			j++;
			} while (j < nAtoms());
			i++;
			} while (i < na);
			*/
		}


		if ((buttonStatus>1) && (chb < 0)) {
			//if fragment 'is freezed' through atom, the atom needs to be deleted
			for (i=0; i<nBonds(); i++) {
				if (getBond(i)->at[0] == cha) getBond(i)->at[0]=na;
				if (getBond(i)->at[1] == cha) getBond(i)->at[1]=na;
			}
			defineAtomConn();
			deleteAtom(cha);                      //atom deleting
		}; //new values for attribute array

		defineAtomConn();
		//Search for badly-connected aromatic bonds
		/* commented to better time. Use alternation bond instead
		if (chb>=0) {
		for (i=0; i<nBonds(); i++) bTested[i]=0;
		do {
		i=0;
		do {                               //Search for badly connected atom
		test=(getAtom(i)->currvalence > maxVal[getAtom(i)->na]);
		i++;
		} while (! (test || (i == nAtoms())));
		if (test) {                   //Creation aromtic bond list
		n=0;
		aList[0]=i-1;
		doubleSearch=true;
		do {
		j=nb;
		test=false;
		if (j < nBonds()) do {
		test=((getBond(j)->at[0] == aList[n]) || (getBond(j)->at[1] == aList[n]));
		if (doubleSearch) test=(test && (getBond(j)->tb == 2)); else test=(test && (getBond(j)->tb == 1));
		test=(test && (bTested[j] == 0));
		j++;
		} while (! (test || (j == nBonds())));
		if (test) {              //Addition to the aromatic list
		bList[n]=j-1;
		bTested[j]=1;
		l=getBond(j-1)->at[0];
		if (l == aList[n]) l=getBond(j-1)->at[1];
		n++;
		aList[n]=l;
		doubleSearch= ! doubleSearch;
		}
		} while (test);
		if (n > 1) {                 //Alternation of aromatic chain
		for (j=1; j<n; j++) {
		if ((int)Math.IEEEremainder(J,2) == 1) fBond.setTB(BList.getValue(J),(byte)1); else fBond.setTB(BList.getValue(J),(byte)2);
		}
		fAtom.setCurrValence(AList.getValue(1),(byte)(fAtom.getCurrValence(AList.getValue(1))-1));
		fAtom.setCurrValence(AList.getValue(N),(byte)(fAtom.getCurrValence(AList.getValue(N))+1));
		Test=true;
		} else {                     //Unable to alternate
		J=Nb;
		Test=false;
		if (J < fBond.maxIndex()) do {                          //Search for double bond
		J=J+1;
		Test=((fBond.getAT(J,1) == AList.getValue(N)) || (fBond.getAT(J,2) == AList.getValue(N)));
		Test=(Test & (fBond.getTB(J) == 2));
		} while (! (Test || (J == fBond.maxIndex())));
		if (Test) {        //Conversion of found double to single bond
		fBond.setTB(J,(byte)1);
		fAtom.setCurrValence(fBond.getAT(J,1),(byte)(fAtom.getCurrValence(fBond.getAT(J,1))-1));
		fAtom.setCurrValence(fBond.getAT(J,2),(byte)(fAtom.getCurrValence(fBond.getAT(J,2))-1));
		Test=true;
		}
		}
		}
		} while (test);
		}
		*/
		return result;
}



void TEditedMolecule::addAsTemplate(TSimpleMolecule& fragmentMol, int thisAN, int smAN, int thisBN, int smBN, bool isAddition) {
	int nAtomsOld,nBondsOld,naDef1,naDef,i;
	int naNum;
	bool test,test1,test2;
	std::vector<int> list(nBonds()+fragmentMol.nBonds()+4);
	std::vector<int> listA(nBonds() + fragmentMol.nBonds() + 4);
	std::vector<int> listB(nBonds() + fragmentMol.nBonds() + 4);
	double r,scale,xOld,yOld,xNew,yNew,xCenter,yCenter;
	double r1;
	double r2;
	double xu1;
	double yu1;
	double xu2;
	double yu2;
	double cAngle;
	double rxTemp1,ryTemp1;
	int mouseButton=0;
	int n;
	double amin;



	if (fragmentMol.nAtoms() == 0) return;
	naDef1=smAN;
	test1=true;
	test2=false;
	xCenter=0;
	yCenter=0;
	if (thisBN >= 0) {
		naDef1=fragmentMol.getBond(smBN)->at[0];
		test1=false; test2=true;
	}

	scale=1;
	xOld=0;
	yOld=0;
	test=fragmentMol.makeFragment(naNum,list,naDef1,-1); //creation of template's fragment
	//Scale definitions

	if (thisAN >= 0) { //connection through atoms
		xCenter=this->getAtom(thisAN)->rx;
		yCenter=this->getAtom(thisAN)->ry;
		naDef=thisAN;
		scale=1;
		r1=0;
		r2=0;
		if (fragmentMol.getAtom(naDef1)->nb > 0) {
			amin=100000000; n=1;
			for (i=0; i<fragmentMol.getAtom(naDef1)->nb; i++) {
				rxTemp1=fragmentMol.getAtom(naDef1)->rx-fragmentMol.getAtom(fragmentMol.getAtom(naDef1)->ac[i])->rx;
				ryTemp1=fragmentMol.getAtom(naDef1)->ry-fragmentMol.getAtom(fragmentMol.getAtom(naDef1)->ac[i])->ry;
				r1=sqrt(rxTemp1*rxTemp1+ryTemp1*ryTemp1);
				if(r1 < amin) {
					amin=r1;
					n=i;
				};
			};
			rxTemp1=fragmentMol.getAtom(naDef1)->rx-fragmentMol.getAtom(fragmentMol.getAtom(naDef1)->ac[n])->rx;
			ryTemp1=fragmentMol.getAtom(naDef1)->ry-fragmentMol.getAtom(fragmentMol.getAtom(naDef1)->ac[n])->ry;
			r1=sqrt(rxTemp1*rxTemp1+ryTemp1*ryTemp1);
		} else r1=DEFAULTBONDLENGTH;

		if (this->getAtom(naDef)->nb > 0) {
			amin=100000000; n=1;
			for (i=0; i<this->getAtom(naDef)->nb; i++) {
				rxTemp1=this->getAtom(naDef)->rx-this->getAtom(this->getAtom(naDef)->ac[i])->rx;
				ryTemp1=this->getAtom(naDef)->ry-this->getAtom(this->getAtom(naDef)->ac[i])->ry;
				r2=sqrt(rxTemp1*rxTemp1+ryTemp1*ryTemp1);
				if(r2 < amin) {
					amin=r2;
					n=i;
				};
			};
			rxTemp1=this->getAtom(naDef)->rx-this->getAtom(this->getAtom(naDef)->ac[n])->rx;
			ryTemp1=this->getAtom(naDef)->ry-this->getAtom(this->getAtom(naDef)->ac[n])->ry;
			r2=sqrt(rxTemp1*rxTemp1+ryTemp1*ryTemp1);
		} else {
			if (this->nBonds() > 0) r2=this->averageBondLength(); else r2=DEFAULTBONDLENGTH;
		}
		if (r1 > 0) scale=r2/r1; else scale=1;
	} else if (thisBN >= 0) {      //connection through bonds
		rxTemp1=fragmentMol.getAtom(fragmentMol.getBond(smBN)->at[0])->rx-fragmentMol.getAtom(fragmentMol.getBond(smBN)->at[1])->rx;
		ryTemp1=fragmentMol.getAtom(fragmentMol.getBond(smBN)->at[0])->ry-fragmentMol.getAtom(fragmentMol.getBond(smBN)->at[1])->ry;
		r1=sqrt(rxTemp1*rxTemp1+ryTemp1*ryTemp1);
		rxTemp1=this->getAtom(this->getBond(thisBN)->at[0])->rx-this->getAtom(this->getBond(thisBN)->at[1])->rx;
		ryTemp1=this->getAtom(this->getBond(thisBN)->at[0])->ry-this->getAtom(this->getBond(thisBN)->at[1])->ry;
		r2=sqrt(rxTemp1*rxTemp1+ryTemp1*ryTemp1);
		if (r1 > 0) scale=r2/r1; else scale=1;
	}
	//Old coordinates definition
	if (smAN >= 0) {
		xOld=fragmentMol.getAtom(smAN)->rx;
		yOld=fragmentMol.getAtom(smAN)->ry;
	}
	if (thisBN >= 0) {
		xOld=(fragmentMol.getAtom(fragmentMol.getBond(smBN)->at[0])->rx+fragmentMol.getAtom(fragmentMol.getBond(smBN)->at[1])->rx)/2;
		yOld=(fragmentMol.getAtom(fragmentMol.getBond(smBN)->at[0])->ry+fragmentMol.getAtom(fragmentMol.getBond(smBN)->at[1])->ry)/2;
	}


	//New coordinate definition
	if (thisBN >= 0) {  //connection through bonds
		xNew=(this->getAtom(this->getBond(thisBN)->at[0])->rx+this->getAtom(this->getBond(thisBN)->at[1])->rx)/2;
		yNew=(this->getAtom(this->getBond(thisBN)->at[0])->ry+this->getAtom(this->getBond(thisBN)->at[1])->ry)/2;
	} else {        //connection through atoms}
		if (isAddition) { //addition of fragment
			this->unitVector(thisAN,xu1,yu1);
			if (xu1 > -0.999) r=yu1/(1+xu1); else if (yu1 > 0) r=1E9; else r=-1E9;
			cAngle=2*180*atan(r)/M_PI;
			if (true) {
				cAngle=(int)(cAngle/30);
				cAngle=30*cAngle;
				r=cAngle*M_PI/180;
				xu1=cos(r);
				yu1=sin(r);
			}
			xNew=r1*scale*xu1+this->getAtom(thisAN)->rx;
			yNew=r1*scale*yu1+this->getAtom(thisAN)->ry;
		} else {                     //freezing of fragment
			if (this->getAtom(thisAN)->nb == 0) {
				xu1=1; yu1=0;
			} else if (this->getAtom(thisAN)->nb == 1) {
				xu1=this->getAtom(thisAN)->rx-this->getAtom(this->getAtom(thisAN)->ac[0])->rx;
				yu1=this->getAtom(thisAN)->ry-this->getAtom(this->getAtom(thisAN)->ac[0])->ry;
				//yu1=this.fAtom.getRY(thisAN)-this.fAtom.getRY(this.fAtom.getAC(thisAN,1));
				r2=sqrt(xu1*xu1+yu1*yu1);
				if (r2 > 1.0E-4) {
					xu1=xu1/r2;
					yu1=yu1/r2;
				} else {xu1=1; yu1=0;}
			} else this->unitVector(thisAN,xu1,yu1);
			if (xu1 > -0.999) r=yu1/(1+xu1); else if (yu1 > 0) r=1E9; else r=-1E9;
			cAngle=2*180*atan(r)/M_PI;
			if (true) {
				cAngle=(int)(cAngle/30);
				cAngle=30*cAngle;
				r=cAngle*M_PI/180;
				xu1=cos(r);
				yu1=sin(r);
			}
			xNew=this->getAtom(thisAN)->rx;
			yNew=this->getAtom(thisAN)->ry;
		}
	}


	//Angle definition
	if (thisAN >= 0) {                             //connection through atoms
		fragmentMol.unitVector(smAN,xu2,yu2);
		xu2=-xu2; yu2=-yu2;
		r1=xu1*xu2+yu1*yu2;
		r2=yu1*xu2-xu1*yu2;
		nBondsOld=this->nBonds();
		nAtomsOld=this->nAtoms();
		if (isAddition) mouseButton=1; else mouseButton=2;

		this->addFragment(fragmentMol,naNum,thisAN,thisBN,smBN,list,
			xOld,yOld,xNew,yNew,scale,r1,r2,mouseButton,false);


		if (isAddition) {
			this->addBond(1,nAtomsOld,thisAN);
		} else nAtomsOld=nAtomsOld-1;
	} else if ((smBN >= 0) && (thisBN >= 0)) {       //connection through bonds
		this->bondUnitVector(thisBN,xu1,yu1);
		fragmentMol.bondUnitVector(smBN,xu2,yu2);
		xu2=-xu2; yu2=-yu2;
		r1=xu1*xu2+yu1*yu2;
		r2=yu1*xu2-xu1*yu2;
		nBondsOld=this->nBonds();
		nAtomsOld=this->nAtoms();
		addFragment(fragmentMol,naNum,thisAN,thisBN,smBN,list,
			xOld,yOld,xNew,yNew,scale,r1,r2,1,false);
	}
}


TEditedMolecule * TEditedMolecule::extractFragment(int atomN, std::vector<int> * enumerator) {

	std::vector<int> list(listarSize());
	std::vector<int>  inverseList(listarSize());
	int i,j,k;
	int na;
	bool test;
	TSingleAtom * sa;
	TSingleBond * sb;
	TEditedMolecule * result=NULL;

	if ((atomN < 0) || (atomN >= nAtoms())) return result;
	if (enumerator != NULL) for (i=0; i<nAtoms(); i++) (*enumerator)[i]=-1;
	for (i=0; i<nAtoms(); i++) inverseList[i]=-1;
	test=makeFragment(na,list,atomN,-1);


	if (na>1) for (i=0; i<(na-1); i++) for (j=i+1; j<na; j++) if (list[i]>list[j]) {
		k=list[i];
		list[i]=list[j];
		list[j]=k;
	};
	if (na > 0) for (i=0; i<na; i++) inverseList[list[i]]=i;

	result= new TEditedMolecule();
	if (na > 0) for (i=0; i<na; i++) {
		sa=this->getAtom(list[i])->clone();
		result->addAtom(sa);
		if (enumerator != NULL) (*enumerator)[list[i]]=i;
	}
	if (nBonds() > 0) for (i=0; i<nBonds(); i++)
		if (inverseList[getBond(i)->at[0]] >= 0) {
			sb=getBond(i)->clone();
			sb->at[0]=inverseList[getBond(i)->at[0]];
			sb->at[1]=inverseList[getBond(i)->at[1]];
			result->addBond(sb);
		}
		return result;
}



int TEditedMolecule::prepareQuery(TSimpleMolecule & sMol) {
	/*Returns error code:
	0  - no error
	-1 - input molecule is not assigned or number of atoms equal zero
	-2 - unconnected fragments were detected
	-3 - an exception take place*/

//	TSimpleMolecule * molecule1;
	TEditedMolecule * molecule1;
	int i,j,k,l;
	bool test,test1,test2,test3;
	int aQ1,aQ2,bNQ;
	int result;
	bool whileCondition,whileCondition1;
	int i1,i2,m,jj1;

	result=-1;
	if (sMol.nAtoms()==0) return result;
	sMol.defineAtomConn();
	sMol.allAboutCycles();
	//Three above line instead of this line: CalculateAllIndeces(SMOl);
	moleculeCopy(sMol);

	result=0;


	fIsQueryPrepare=true;
	for (i=0; i<nAtoms(); i++) { //isotops change
		if ((getAtom(i)->na==104) && (! fIOPT11)) getAtom(i)->na=1;
		if (! fIOPT11) getAtom(i)->iz=0;
	};

	//Initializinf arrays
	queryQHydr.resize(listarSize());
	queryAGer.resize(listarSize());
	queryAQTested.resize(listarSize());
	queryBQTested.resize(listarSize());
	queryBQCounter.resize(listarSize());
	queryCurrentAssignment.resize(listarSize());
	aSTested.resize(listarSize());
	bSTested.resize(listarSize());
	bSTestedStore.resize(listarSize());
	queryStas.resize(listarSize());
	queryBKNew.resize(listarSize());
	structureBKNew.resize(listarSize());
	//initializing
	for (i=0; i<listarSize(); i++) {
		queryBKNew[i].nb=0;
		queryQHydr[i]=0;
		queryAGer[i]=-1;
		queryAQTested[i]=-1;
		queryBQTested[i]=-1;
		queryBQCounter[i]=-1;
		queryCurrentAssignment[i]=-1;
		aSTested[i]=-1;
		bSTested[i]=-1;
		bSTestedStore[i]=-1;
		queryStas.setNB(i,0);
	};
	queryEnum.resize(sMol.nAtoms());
	queryInverse.resize(sMol.nAtoms());
	for (i=0; i<sMol.nAtoms(); i++) {
		queryEnum[i]=-1;
		queryInverse[i]=-1;
	};
	molecule1 = new TEditedMolecule();  //TSimpleMolecule();

	//At this step memory was not allocated to: BK1,BEQ,AEQ. Realloction will be
	//required for three above temporary allocated arrays}
	removeHydrogen(&queryQHydr,&queryEnum);  //Connection is calculated inside RemoveHydr}
	if (fIOPT13) atomBondChange();  //semipolar bond conversion
	defineBondConn(queryBKNew);
	//initial values
	for (i=0; i<nAtoms(); i++) aSTested[i]=getAtom(i)->special;
	l=1;
	if (nAtoms()>2) {
		j=queryBKNew[0].nb;
		for (i=1; i<nAtoms(); i++) {
			//search for maximally coordinated non-carbon
			k=queryBKNew[i].nb;
			test1=k>j;
			test2=(getAtom(l)->na<104) && (getAtom(l)->na!=6);
			test3=(getAtom(i)->na<104) && (getAtom(i)->na!=6);
			if ((test1 && (! test2)) || (test3 && (! test2))
				|| (test1 && test2 && test3)) {
					j=k; l=i;
			};
		};
	};

	queryAQTested[0]=0; //Sequence for query assignment begin from the atom
	aQ1=0; bNQ=-1; test=true;
	if (nBonds()>0) while (test) {  //Define sequence to match query and structure
		i=0;
		whileCondition=true;
		while (whileCondition) {
			test=(queryAQTested[i] != -1);
			if (test) {
				j=0;
				whileCondition1=true;
				if (queryBKNew[i].nb>0) while (whileCondition1) {
					k=queryBKNew[i].adjusted[j];
					j++;
					test=(bSTested[k]==-1);
					whileCondition1=test || (j>=queryBKNew[i].nb);
					whileCondition1=! whileCondition1;
				};
				if (queryBKNew[i].nb==0) test=false;
			};
			i++;
			whileCondition=test || (i==nAtoms());
			whileCondition=! whileCondition;
		};
		if (test) { //connected atom found
			i--;  //zero-index
			aQ2=i;
			whileCondition=true;
			while (whileCondition) {
				j=0;
				//search for query bond appended to AQ2, which has not been assigned yet
				whileCondition1=true;
				k=0;
				test1=true;
				if (queryBKNew[aQ2].nb>0) while (whileCondition1) {
					k=queryBKNew[aQ2].adjusted[j];
					test1=(bSTested[k]==-1);
					j++;
					whileCondition1=test1 || (j>=queryBKNew[aQ2].nb);
					whileCondition1=! whileCondition1;
				} else {
					k=0;
					test1=false;
				};
				if (test1)  {  //found
					bNQ++; bSTested[k]=bNQ; //store, that BNQ has been assigned
					l=getBond(k)->at[1];
					//atoms enumeration
					if (l==aQ2) {
						l=getBond(k)->at[0];
						getBond(k)->at[1]=l;
						getBond(k)->at[0]=aQ2;
					};
					test1=(queryAQTested[l]==-1);
					//Make AGER array for end-of cycle BNQ bond}
					if (test1) {
						queryAGer[bNQ]=l;
						aQ1=aQ1+1;
						queryAQTested[l]=aQ1;
						aQ2=l;
					} else aQ2=-1;
				};
				whileCondition=(! test1) || (aQ2==-1);
				whileCondition=! whileCondition;
			};
		};
	};
	//Detection of unconnected fragment}
	for (i=0; i<nAtoms(); i++) if (queryAQTested[i]==-1) result=-2;
	if (result>=0) {
		//Creation of inverse enumeration
		queryInverse.resize(queryAQTested.size());
		for (j=0; j<queryInverse.size(); j++) queryInverse[j]=queryAQTested[j];
		//Enumeration of atoms and bonds in query in order of assignment}
		molecule1->moleculeCopy(*this);
		for (i=0; i<nAtoms(); i++) molecule1->getAtom(queryAQTested[i])->atomCopy(this->getAtom(i));
		for (i=0; i<nAtoms(); i++) queryBQTested[queryAQTested[i]]=queryQHydr[i];
		for (i=0; i<queryQHydr.size(); i++) queryQHydr[i]=queryBQTested[i]; 
		//Do NOT enumeate QueryEnum - it will be done later throgh QueryInverse!!!}
		if (nBonds()>0) for (i=0; i<nBonds(); i++) if (queryAGer[i] != -1) queryAGer[i]=queryAQTested[queryAGer[i]];
		if (nBonds()>0) for (i=0; i<nBonds(); i++) molecule1->getBond(bSTested[i])->bondCopy(this->getBond(i));
		bSTestedStore.resize(bSTested.size());
		for (i=0; i<bSTestedStore.size(); i++) bSTestedStore[i]=bSTested[i];
		if (nBonds()>0) for (i=0; i<nBonds(); i++) {
			molecule1->getBond(i)->at[0]=queryAQTested[molecule1->getBond(i)->at[0]];
			molecule1->getBond(i)->at[1]=queryAQTested[molecule1->getBond(i)->at[1]];
		};
		//Enumeration of
		this->moleculeCopy(*molecule1);
		defineAtomConn();
		//Store NO OTHER atom attribute
		for (i=0; i<nAtoms(); i++) getAtom(queryAQTested[i])->special=aSTested[i];
		defineBondConn(queryBKNew);
		j=0;


        queryTestSpace=false;
        if (nBonds() > 0) for (i=0; i<nBonds(); i++) { //Make STAS array
 
		  if (getBond(i)->special > 2) { //CIS/TRANS label on bond
            i1=getBond(i)->at[0];
            i2=getBond(i)->at[1];
			for (j=0; j<queryBKNew[i1].nb; j++) if (queryBKNew[i1].adjusted[j] != i) {
              //for each bond, connected with first atom..}
              m=queryBKNew[i1].adjusted[j];
              for (k=0; k<queryBKNew[i2].nb; k++) if ((queryBKNew[i2].adjusted[k] != i) && (queryStas.getNB(m) < SCALAR_MAX)) {
                //..takes a bond, connected with second atom}
                jj1=queryStas.getNB(m);
                queryStas.setRef(m,jj1,queryBKNew[i2].adjusted[k]); //auxiliary bonds}
                queryStas.setRel(m,jj1,i);
                queryStas.setScalar(m,jj1,sproduct(*this,i,m,queryBKNew[i2].adjusted[k])); //Code}
                queryStas.setNB(m,jj1+1);
			  };
			};
            for (j=0; j<queryBKNew[i2].nb; j++) if (queryBKNew[i2].adjusted[j] != i) {
              //for each bond, connected with second atom..}
              m=queryBKNew[i2].adjusted[j];
			  for (k=0; k<queryBKNew[i1].nb; k++) if ((queryBKNew[i1].adjusted[k] != i) && (queryStas.getNB(m) < 50)) {
                //takes a bond, connected with the first atom}
                jj1=queryStas.getNB(m);
                queryStas.setRef(m,jj1,queryBKNew[i1].adjusted[k]);
                queryStas.setRel(m,jj1,i);
                queryStas.setScalar(m,jj1,sproduct(*this,i,m,queryBKNew[i1].adjusted[k]));
                queryStas.setNB(m,jj1+1);
			  };
			};
            queryTestSpace=true;
	      };
	    };
		//Search for last non-special atom in query}
		queryStereoQ=stereoBondChange(); //Stereo bond conversion
	};
	delete(molecule1);
	return result;
};


int TEditedMolecule::sproduct(TSimpleMolecule & molecule, int br, int i1, int i2){
/*Calculate if two bonds I1 and I2 are directed on single side of straight line
 (equation of straight line is determined by bond BR) or back. The structure
 is described by ATOM and BOND, the bonds I1 and I2 having common atom with BR.
 On output function have values:
 =0-It is impossible to calculate
 =1-at the same side
 =2-at opposite side*/
  int result=0;
  int an[4];
  double rx[4];
  double ry[4];
  int i;
  double r1,r2;

  an[0]=molecule.getBond(br)->at[0];
  an[1]=molecule.getBond(br)->at[1];
  if ((molecule.getBond(i1)->at[0] != an[0]) && (molecule.getBond(i1)->at[0] != an[1])) {
    an[2]=molecule.getBond(i1)->at[0];
  } else {
    an[2]=molecule.getBond(i1)->at[1];
  }
  if ((molecule.getBond(i2)->at[0] != an[0]) && (molecule.getBond(i2)->at[0] != an[1])) {
    an[3]=molecule.getBond(i2)->at[0];
  } else {
    an[3]=molecule.getBond(i2)->at[1];
  }
  for (i=0; i<4; i++) {
	rx[i]=molecule.getAtom(an[i])->rx;
	ry[i]=molecule.getAtom(an[i])->ry;
  };
  //to the coordinate system centered at CR[1].AN
  for (i=1; i<4; i++) {
	rx[i]=rx[i]-rx[0];
	ry[i]=ry[i]-ry[0];
  };
  //distance (including sign!) calculation
  r1=rx[1]*ry[2]-ry[1]*rx[2];
  r2=rx[1]*ry[3]-ry[1]*rx[3];
  if ((r1 == 0) || (r2 == 0)) return result;
  if ((r1*r2) > 0) result=1; else result=2;

  return result;
};

/*
    if (FBond[I1].AT[1]<>CR[1].AN) and (FBond[I1].AT[1]<>CR[2].AN) then
      CR[3].AN:=FBond[I1].AT[1] else CR[3].AN:=FBond[I1].AT[2];
    if (FBond[I2].AT[1]<>CR[1].AN) and (FBond[I2].AT[1]<>CR[2].AN) then
      CR[4].AN:=FBond[I2].AT[1] else CR[4].AN:=FBond[I2].AT[2];
    for I:=1 to 4 do begin {coordinates}
      CR[I].Rx:=FAtom[CR[I].AN].Rx;
      CR[I].Ry:=FAtom[CR[I].AN].Ry;
    end;
    for I:=2 to 4 do begin {to the coordinate system centered at CR[1].AN}
      CR[I].Rx:=CR[I].Rx-CR[1].Rx;
      CR[I].Ry:=CR[I].Ry-CR[1].Ry;
    end;
    {distance (including sign!) calculation}
    R1:=CR[2].Rx*CR[3].Ry-CR[2].Ry*CR[3].Rx;
    R2:=CR[2].Rx*CR[4].Ry-CR[2].Ry*CR[4].Rx;
    if (R1=0) or (R2=0) then begin Result:=0; Exit; end;
    if ((R1>0) and (R2>0)) or ((R1<0) and (R2<0)) then Result:=1 else Result:=2;
  end;
end;
*/

bool TEditedMolecule::getEquivalentValue(int indexQuery, int indexStructure, int nQueryElements, const std::vector<bool>&matrix) {
	bool result = false;
	int index = indexStructure*nQueryElements + indexQuery;

	result = matrix[index];
	return result;
};

void TEditedMolecule::setEquivalentValue(int indexQuery, int indexStructure, int nQueryElements, bool value, std::vector<bool>&matrix) {
	int index = indexStructure*nQueryElements + indexQuery;

	matrix[index]=value;
};


void TEditedMolecule::directBondAss(int& bnq, bool& test, bool& test1, /*const bool beq[NBONDSMAX][NQUERYMAX], 
	const bool aeq[NATOMSMAX][NQUERYMAX],*/ const std::vector<bool> &beq, const std::vector<bool> &aeq, std::vector<int>& bqcounter, std::vector<int>& aqtested,
	std::vector<int>& bstested, std::vector<int>& bqtested, std::vector<int>& astested, 
	const std::vector<int> ager, const std::vector<adjustedlist> & bqconn, const std::vector<adjustedlist>& bsconn, TSimpleMolecule * smol,
	TScalar& stas) {
		//Assign a structure atom to query atom. Variables:
		//BNQ-query bond's number, which needs to be assigned
		//BEQ,AEQ-boolean matrix, contains list of equivalent bonds
		//BQCONN,BSCONN-for each atom containes list of connected bonds in query and
		//              structure respectively.
		//QBOND,SBOND-bond's attributes for query and structure respectively
		//ASTESTED- I-th element of the array contains the query atom's number, which
		//          has been assigned to I-th atom of structure (=0-no assignment)
		//AQTESTED- I-th element of the array contains the structure atom's number,
		//          which has been assigned to I-th atom of query (=0-no assignment)
		//BSTESTED- I-th element of the array contains the query bond's number, which
		//          has been assigned to I-th bond of structure (=0-no assignment)
		//BQTESTED- I-th element of the array contains the structure bond's number,
		//          which has been assigned to I-th bond of query (=0-no assignment)
		//BQCOUNTER-for each query bond contains bond's number in array BSCONN(!)-con-
		//          nected to corresponding structure atom bonds. The structure bonds
		//          in array BSCONN with less than or equal number have been tested on
		//          assignment, with above number-no.
		//AGER-for each query bond contains query atom's number, which must be generated
		//     at assignment of this bond to structure or zero, if no atom may be gene-
		//     rated. This array is used for cyclic conditions definition (last bond in
		//     the cycle must be created between already-defined query atoms-no genera-
		//     tion).
		//TEST-on output contains TRUE, if assignment was successfull, FALSE otherwise
		//TEST1-on output contains TRUE if all bonds, connected to last assigned atom
		//      in structure will be tested, FALSE otherwise. If TEST1=TRUE and TEST=
		//      FALSE it means, that last atom in structure had been unproperly assigned
		//      -backstep is required to reassign the atom.*/
		int bns,bs,as1,as2,aq2,aq1;
		bool whiletest;
		int i,j,i1,i2;

		test=false;
		aq1=getBond(bnq)->at[0];      //first query atom, already assigned
		as1=aqtested[aq1];  //corresponding first structure atom
		if (as1<0) return;
		aq2=getBond(bnq)->at[1];      //second query atom}
		bs=bqcounter[bnq]; 
		bns=0; as2=0;
		whiletest=true;
		if (bs<bsconn[as1].nb) while (whiletest) {
			bns=bsconn[as1].adjusted[bs];  //structure bond, assigned to query BNQ}
			//if structure bond hasn't already been assigned, and query and structure
			//bonds are equivalent, then...
			if (bstested[bns]<0) if (/*beq[bns][bnq]*/getEquivalentValue(bnq,bns,this->nBonds(),beq)) {
				as2=smol->getBond(bns)->at[0];                  //second structure atom}
				if (as2==as1) as2=smol->getBond(bns)->at[1];
				test=false;
				//if no new structure atom assignment must be generated-checking, if AS2
				// equal already defined atom for query atom AQ2}
				if ((ager[bnq]<0) && (astested[as2]>=0)) test=(aqtested[aq2]==as2);
				//if new atom in structure must be assigned-testing equivalence with
				//corresponding query
				if ((ager[bnq] >= 0) && (astested[as2] < 0)) test = getEquivalentValue(aq2, as2, this->nAtoms(), aeq);  //aeq[as2][aq2];
				i=0;
				if (test && (stas.getNB(bnq) > 0)) while (test && (i<stas.getNB(bnq))) {
                  i1=stas.getRef(bnq,i); //query bond numbers
                  i2=stas.getRel(bnq,i);
                  i1=bqtested[i1];     //corresponding structure bond numbers
                  i2=bqtested[i2];
                  if ((i1 >= 0) && (i2 >= 0)) {
                    j=sproduct(*smol,i2,i1,bns);
                    test=(j == stas.getScalar(bnq,i));
				  };
                  i++;
				};
/*
      if Test and SterT then if STAS[BNQ].NB>0 then begin
        {stereo check, CIS/TRANS query bond attribute}
        I:=0;
        repeat    {for all pairs in structure}
          I:=I+1;
          I1:=STAS[BNQ].Ref[I]; {query bond numbers}
          I2:=STAS[BNQ].Rel[I];
          I1:=BQTESTED[I1];     {corresponding structure bond numbers}
          I2:=BQTESTED[I2];
          if (I1>0) and (I2>0) then begin
            J:=Sproduct(SMol,I2,I1,BNS);
            TEST:=J=STAS[BNQ].Scalar[I];
          end;
        until (I=STAS[BNQ].NB) or (not TEST);
      end;
*/


			};
			bs=bs+1;
			whiletest=(bs==bsconn[as1].nb) || test;
			whiletest=! whiletest;
		}; //until(BS=BSCONN[AS1].NB) or TEST;
		//until all list of structure bonds will be exausted or success in assignment
		if (test) { //Success
			if (ager[bnq]>=0) {    //if new atom has been assigned-store assignment}
				aqtested[aq2]=as2; astested[as2]=aq2;
			};
			bstested[bns]=bnq; //store assignment of bonds}
			bqtested[bnq]=bns;
		};
		bqcounter[bnq]=bs;    //last used structure bond's number in SCONN array}
		test1=(bs==bsconn[as1].nb); //test, if SCONN for atom given is exausted}
};

bool TEditedMolecule::allQueryPresent(const std::vector<int> qA, const std::vector<int> qB,
	int nA, int nB) {
		//the function returns TRUE if for all query atoms and bonds the equivalent
		//structure atom (bond) can be found, FALSE otherwise. QA-for each query atom
		//contains '1' if a structure atom may be associated with the query, '0' other-
		//wise. The same is true for bond's list QB. NA-total number of query atoms,
		//NB-total number of query bonds
		int i;
		bool test,whiletest;

		i=0;
		if (nA<0) return false;
		whiletest=true;
		test=false;
		while (whiletest) {
			test=qA[i]==1;
			i++;
			whiletest=(! test) || (i==nA);
			whiletest=! whiletest;
		};
		if (test && (nB>=0)) {
			i=0;
			whiletest=false;
			while (whiletest) {
				test=qB[i]==1;
				i++;
				whiletest=(! test) || (i==nB);
				whiletest=! whiletest;
			};
		};
		return test;
};

#define RMG_LABEL_MAX 8

bool TEditedMolecule::fragmentSearch(TSimpleMolecule * molecule1, std::vector<int>* bondLabel, bool useFragmentNo, bool useAddInfo) {

	int labelsThis[RMG_LABEL_MAX];
	int labelsMol[RMG_LABEL_MAX];

	int cycleNumber;
	int j,k,l,m,mm;
	bool test;
	bool test1;
	bool test2,test3,stereoS;
	int aq1,aq2,as1,as2,j1,i1;
	int ii;
	bool result=false;
	bool whiletest1,whiletest2;
	std::vector<bool> unsaturated(nAtoms()+molecule1->nAtoms());
	//bool			aEQ[NATOMSMAX][NQUERYMAX];
	//bool			bEQ[NBONDSMAX][NQUERYMAX];
	std::vector<bool> aEQ(nAtoms() * molecule1->nAtoms()+1);
	std::vector<bool> bEQ(nBonds() * molecule1->nBonds()+1);


	if ((molecule1==NULL) || (! fIsQueryPrepare)) return result;
	if (molecule1->nAtoms()==0) return result;

	for (j=0; j<RMG_LABEL_MAX; j++) {
	  labelsThis[j]=0;
	  labelsMol[j]=0;
	};

	if (useFragmentNo) {
	  for (j=0; j<nAtoms(); j++) {
		k=getAtom(j)->fragIndex;
		if ((k > 0) && (k <= RMG_LABEL_MAX)) labelsThis[k-1]=j+1;  //zero index == 1
	  };
	  for (j=0; j<molecule1->nAtoms(); j++) {
		k=molecule1->getAtom(j)->fragIndex;
		if ((k > 0) && (k <= RMG_LABEL_MAX)) labelsMol[k-1]=j+1;
	  };
	  //checking if all present...
	  for (j=0; j<RMG_LABEL_MAX; j++) if ((labelsThis[j] > 0) && (labelsMol[j] == 0)) return result;
	};

	if (molecule1->listarSize()>aSTested.size()) aSTested.resize(molecule1->listarSize());
	if (molecule1->listarSize()>bSTested.size()) bSTested.resize(molecule1->listarSize());
	if (molecule1->listarSize()>queryQHydr.size()) queryQHydr.resize(molecule1->listarSize());
	if (molecule1->listarSize()>queryAQTested.size()) queryAQTested.resize(molecule1->listarSize());
	if (molecule1->listarSize()>queryBQTested.size()) queryBQTested.resize(molecule1->listarSize());
	if (molecule1->listarSize()>queryCurrentAssignment.size()) queryCurrentAssignment.resize(molecule1->listarSize());
	if (this->listarSize()>queryAQTested.size()) queryAQTested.resize(this->listarSize());
	if (molecule1->listarSize()>structureBKNew.size()) structureBKNew.resize(molecule1->listarSize());
	if (this->listarSize()>queryBKNew.size()) queryBKNew.resize(this->listarSize());

	//for (int i=0; i<molecule1->nAtoms(); i++) for (j=0; j<nAtoms(); j++) aEQ[i][j]=false;
	//for (int i=0; i<molecule1->nBonds(); i++) for (j=0; j<nBonds(); j++) bEQ[i][j]=false;
	for (int i = 0; i < aEQ.size(); i++) aEQ[i] = false;
	for (int i = 0; i < bEQ.size(); i++) bEQ[i] = false;
	cycleNumber=0;
	//Formula filter
	test2=false;
	/*
	if (fBruttoFormula.nD>0) {
	j=0;
	if (molecule1.fBruttoFormula.nD>0) {
	whiletest1=true;
	while (whiletest1) {
	j++; k=0;
	test1.value=(fBruttoFormula.eList[j]==104) && (! fIOPT11);
	whiletest2=true;
	if ((fBruttoFormula.eList[j] != 1) && (! test1.value)) while (whiletest2) {
	k++;
	test2=((fBruttoFormula.eList[j]==molecule1.fBruttoFormula.eList[k])
	&& (fBruttoFormula.eNumber[j]<=molecule1.fBruttoFormula.eNumber[k]));
	whiletest2=test2 || (k==molecule1.fBruttoFormula.nD);
	whiletest2=! whiletest2;
	} else test2=true;
	whiletest1=(! test2) || (j==fBruttoFormula.nD);
	whiletest1=! whiletest1;
	};
	} else test2=false;
	} else test2=true;
	*/ test2=true;   //turn of molecular formula filter

	molecule1->fIOPT11=fIOPT11;
	molecule1->fIOPT12=fIOPT13;
	molecule1->fIOPT13=fIOPT13;
	if (test2) {

		//{R/S/Z/E description are removed: they are not used in substructure search}
		if (fIOPT13) molecule1->atomBondChange(); //Semipolar bond conversion
		molecule1->defineAtomConn();
		molecule1->defineBondConn(structureBKNew);
		stereoS=molecule1->stereoBondChange(); //Stereo bond conversion}
		for (j=0; j<molecule1->nAtoms(); j++) {
			//Isotops conversion}
			if ((molecule1->getAtom(j)->na==104) && (! fIOPT11)) molecule1->getAtom(j)->na=1;
			if (! fIOPT11) molecule1->getAtom(j)->iz=0;
			bSTested[j]=0;  //bsTested contains nu. explicit hydrogens - so 0 is fine!
		};
		for (j=0; j<molecule1->nAtoms(); j++) if (molecule1->getAtom(j)->na==1)
			if (molecule1->getAtom(j)->nb>0) for (k=0; k<molecule1->getAtom(j)->nb; k++) {
				//Formation of list explicitly-defined hydrogens}
				l=molecule1->getAtom(j)->ac[k];
				bSTested[l]=bSTested[l]+1;
			};
		for (k=0; k<nAtoms(); k++) aSTested[k]=-1;
		if (nAtoms()==1) {
			//Partial case - if single atom was defined}
			test2=false;
			j=0;
			whiletest1=true;
			if (molecule1->nAtoms()>0) while (whiletest1) {
				test2=TSingleAtom::atomEquivalent(molecule1->getAtom(j),getAtom(0),bSTested[j],queryQHydr[0],fIOPT10,fIOPT11,
				compareHard,rejectAnyInStructure,false,useAddInfo);
				//Addition for aromatic search}
				if (test2 && ((getAtom(0)->special & AROMATIC_MASK)!=0) ){
					test2=false;
					for (m=0; m<=structureBKNew[j].nb; m++) {
						mm=structureBKNew[j].adjusted[m];
						if ((molecule1->getBond(mm)->db==2) || (molecule1->getBond(mm)->db==3)) {
							test2=true;
							break;
						};
					};
				};
				//End addition}
				if (test2) {
					if (queryEnum.size()<1) queryEnum.resize(1);
					queryEnum[0]=0;
					if (queryAQTested.size()<1) queryAQTested.resize(1);
					queryAQTested[0]=j;
				};
				j++;
				whiletest1=test2 || (j==molecule1->nAtoms());
				whiletest1=! whiletest1;
			};
		} else if (nBonds()==0) test2=false; else {

  	      if (useFragmentNo) {   //atoms numbers enumeration - after possible removing of default hydrogens
	        for (j=0; j<RMG_LABEL_MAX; j++) labelsMol[j]=0;
	        for (j=0; j<molecule1->nAtoms(); j++) {
		      k=molecule1->getAtom(j)->fragIndex;
		      if ((k > 0) && (k <= RMG_LABEL_MAX)) labelsMol[k-1]=j+1;
	        };
		  };

			//General case - substructure search
			if (fIncludedList != NULL) {
				for (j = 0; j < molecule1->nAtoms(); j++) for (k = 0; k < nAtoms(); k++) setEquivalentValue(k, j, nAtoms(), false, aEQ); //aEQ[j][k]=false;
				for (k=0; k<nAtoms(); k++) {
					j=(*fIncludedList)[k];
					setEquivalentValue(k, j, nAtoms(),true, aEQ);//aEQ[j][k]=true;
					aSTested[k]=1;
				};
			} else {
                for (j=0; j<molecule1->nAtoms(); j++) unsaturated[j]=false;
                for (j=0; j<molecule1->nBonds(); j++) {
				  if ((molecule1->getBond(j)->tb == 2) || (molecule1->getBond(j)->tb == 3)) {
					unsaturated[molecule1->getBond(j)->at[0]]=true;
                    unsaturated[molecule1->getBond(j)->at[1]]=true;
			      };
		        };

				for (j=0; j<molecule1->nAtoms(); j++) for (k=0; k<nAtoms(); k++) {
					//Creation of atom-equivalent matrix}

					test3=TSingleAtom::atomEquivalent(molecule1->getAtom(j),getAtom(k),bSTested[j],
						queryQHydr[k],fIOPT10,fIOPT11,compareHard,rejectAnyInStructure,unsaturated[j],useAddInfo) && (getAtom(k)->nb<=molecule1->getAtom(j)->nb);

					if (test3 && useFragmentNo) {
					  m=getAtom(k)->fragIndex;
					  if ((m > 0) && (m <= RMG_LABEL_MAX)) {
						m=labelsMol[m-1]-1;
						test3=(m == j);
					  };
					};
					//Addition for aromatic search
					if (test3 && ((getAtom(k)->special & AROMATIC_MASK) !=0)) {
						test3=false;
						for (m=0; m<structureBKNew[j].nb; m++) {
							mm=structureBKNew[j].adjusted[m];
							if ((molecule1->getBond(mm)->db==2) || (molecule1->getBond(mm)->db==3)) {
								test3=true;
								break;
							};
						};
					};
					//aEQ[j][k]=test3;
					setEquivalentValue(k, j, nAtoms(),test3, aEQ);
					if (test3) aSTested[k]=1;
				};
			};

			for (j=0; j<nBonds(); j++) queryAQTested[j]=0;
			for (j=0; j<molecule1->nBonds(); j++) for (k=0; k<nBonds(); k++) {
				//Creation of bond-equivalent matrix
				aq1=getBond(k)->at[0];   aq2=getBond(k)->at[1];
				as1=molecule1->getBond(j)->at[0]; as2=molecule1->getBond(j)->at[1];
				test2=(getEquivalentValue(aq1,as1,nAtoms(),aEQ) && getEquivalentValue(aq2,as2,nAtoms(),aEQ)) || (getEquivalentValue(aq2, as1, nAtoms(), aEQ) && getEquivalentValue(aq1, as2, nAtoms(), aEQ));  
				//(aEQ[as1][aq1] && aEQ[as2][aq2]) || (aEQ[as1][aq2] && aEQ[as2][aq1]);
				if (test2) {
					test3=TSingleBond::bondEquivalent(molecule1->getBond(j),getBond(k),compareHard);
					setEquivalentValue(k, j, nBonds(), test3, bEQ);  // bEQ[j][k]=test3;
					if (test3) queryAQTested[k]=1;
				} else setEquivalentValue(k, j, nBonds(), false, bEQ); //bEQ[j][k]=false;
			};
			//Check, if all query atoms has a partner in structure
			test2=allQueryPresent(aSTested,queryAQTested,nAtoms(),nBonds());
			if (test2) {
				j1=0;

				for (k=0; k<molecule1->nAtoms(); k++) if (getEquivalentValue(0, k, nAtoms(),aEQ)) { //(aEQ[k][0]) {
					//Collection of structure atoms, which may be associated with first query atom
					j=0;
					whiletest1=true;
					while (whiletest1) {
						//test=bEQ[structureBKNew[k].adjusted[j]][0];
						test=getEquivalentValue(0, structureBKNew[k].adjusted[j], nBonds(), bEQ);
						if (test) {
							queryCurrentAssignment[j1]=k;
							j1++;
						};
						j++;
						whiletest1=test || (j==structureBKNew[k].nb);
						whiletest1=! whiletest1;
					};
				};
				i1=0;
				whiletest1=true;
				if (j1>0) while (whiletest1) { //for each structure, which may be assigned to 1-st query
					for (j=0; j<nAtoms(); j++) queryAQTested[j]=-1; //array initializing
					for (j=0; j<nBonds(); j++) queryBQCounter[j]=0;
					for (j=0; j<nBonds(); j++) queryBQTested[j]=-1;
					for (j=0; j<molecule1->nBonds(); j++) bSTested[j]=-1;
					for (j=0; j<molecule1->nAtoms(); j++) aSTested[j]=-1;
					queryAQTested[0]=queryCurrentAssignment[i1];
					aSTested[queryCurrentAssignment[i1]]=0;
					ii=0;
					whiletest2=true;
					while (whiletest2) { //start recursion
						directBondAss(ii,test,test1,bEQ,aEQ,queryBQCounter,queryAQTested,bSTested,
							queryBQTested,aSTested,queryAGer,queryBKNew,structureBKNew,molecule1,queryStas);
						if ((! test) && (ii>=1)) {
							//previous atom were badly assigned-backstep
							queryBQCounter[ii]=0;
							k=queryAGer[ii-1];
							if (k>=0) {
								l=queryAQTested[k]; queryAQTested[k]=-1; aSTested[l]=-1;
							};
							bSTested[queryBQTested[ii-1]]=-1;
							queryBQTested[ii-1]=-1;
							ii=ii-2;
							test=true;
						};
						ii=ii+1; //query bond counter
						whiletest2=(ii==nBonds()) || ((! test) && (ii==1));
						whiletest2=! whiletest2;
					};
					i1++;
					test2=(ii==nBonds()); //Checking, if success has been reached
					whiletest1=((j1==i1) || test2);
					whiletest1=! whiletest1;
				} else test2=false;
			};
		};
	};

	if ((bondLabel != NULL) && (nBonds()>0)) {
		bondLabel->resize(molecule1->nBonds());
		if (test2) for (j=0; j<molecule1->nBonds(); j++) if (bSTested[j]>=0) (*bondLabel)[j]=1; else (*bondLabel)[j]=0;
	};
	result=test2;
	return result;
};

double  TEditedMolecule::distance(int atom1, int atom2) const
{
	double		dist = 0.0;
	dist = sqrt(pow(fAtom[atom1]->rx - fAtom[atom2]->rx, 2) +
					pow(fAtom[atom1]->ry - fAtom[atom2]->ry, 2));
	return dist;
}

int TEditedMolecule::makePoly(int size, int cha1, int chb1, const double& xC, const double& yC) {
	//Draw a cycle in structure. The procedure is called by icon the Draw, DrawCycle command
	int				KCycle = 0, AtStart = 0, AtEnd = 0;
	int				Na,nAt;
	double			C, S, X1, Y1, X2, Y2;
	double			XCenter = 0.0, YCenter = 0.0;
	double			Radius, Vx, Vy, R, RL;
	double			aBondLength;
	bool				Test;
	int				I, J, K, L;
	TSingleAtom		*SA;
	TSingleBond		*SB;
	Rect				Rect;

	nAt = nAtoms();
	aBondLength = averageBondLength();
	if(averageBondLength() < 1.0)
		aBondLength = DEFAULTBONDLENGTH;
	
	//!!!!!!!!!Add check!
	//if (NAtomsMax-Sel.TransferNumber)<FAtom.NAtoms then begin
	//  MessageDlg(LoadResourceString(2081),mtError,[mbOK],0);
	//  Exit;
	//end;
	//if (NBondsMax-Sel.TransferNumber)<FBond.NBonds then begin
	//  MessageDlg(LoadResourceString(2082),mtError,[mbOK],0);
	//  Exit;
	//end;

	if (cha1 < 0 && chb1 < 0) {  //{no atom or bond selection-simple put}
		XCenter = xC;
		YCenter = yC;

		//Radius = DEFAULTBONDLENGTH / (2 * sin(M_PI / size));
		Radius = aBondLength / (2 * sin(M_PI / size));
		Vx = XCenter - Radius;
		Vy = YCenter;
		KCycle = size;
	};

	if(cha1 >= 0) {  //{atom was selected-connection through atom (spiro)}
		if(fAtom[cha1]->nb > CONNMAX - 2) {
			//MessageDlg(LoadResourceString(2083),mtError,[mbOK],0);
			return 1;
		};
		unitVector(cha1,Vx,Vy); //{best direction for cycle formation}
		if(fAtom[cha1]->nb > 0) 
			R = distance(cha1, fAtom[cha1]->ac[0]);
		else
			R = aBondLength;
			//R = DEFAULTBONDLENGTH; //!!!!!!!!!!!!!

		Radius = R / (2*sin(M_PI / size));
		XCenter = fAtom[cha1]->rx + Vx * Radius;
		YCenter = fAtom[cha1]->ry + Vy * Radius;
		Vx = fAtom[cha1]->rx;
		Vy = fAtom[cha1]->ry;
		KCycle = size - 1;
	}; 

	if(chb1 >= 0) { // {bond was selected-connection through bond (condenced)}
		if(fAtom[fBond[chb1]->at[0]]->nb >= CONNMAX || fAtom[fBond[chb1]->at[1]]->nb >= CONNMAX) {
			//MessageDlg(LoadResourceString(2083),mtError,[mbOK],0);
			return 1; //{maximal coordinational number was exceeded}
		};
		bondUnitVector(chb1,Vx,Vy); //{optimal direction}
		R = distance(fBond[chb1]->at[0], fBond[chb1]->at[1]);

		Radius = R / (2 * sin(M_PI / size) / cos(M_PI / size));
		XCenter = (fAtom[fBond[chb1]->at[0]]->rx + fAtom[fBond[chb1]->at[1]]->rx) / 2 + Radius * Vx;
		YCenter = (fAtom[fBond[chb1]->at[0]]->ry + fAtom[fBond[chb1]->at[1]]->ry) / 2 + Radius * Vy;
		AtStart = fBond[chb1]->at[0];
		AtEnd   = fBond[chb1]->at[1];
		Vx = fAtom[AtStart]->rx;
		Vy = fAtom[AtStart]->ry;
		C = cos(2 * M_PI / size);
		S = sin(2 * M_PI / size);

		X1 = Vx - XCenter;
		Y1 = Vy - YCenter;
		X2 =  X1 * C + Y1 * S + XCenter;
		Y2 = -X1 * S + Y1 * C + YCenter;

		R = sqrt(pow(X2 - fAtom[AtEnd]->rx, 2) + pow(Y2 - fAtom[AtEnd]->ry, 2));
		if(R < 0.01) {
			AtStart = fBond[chb1]->at[1];
			AtEnd   = fBond[chb1]->at[0];
			Vx = fAtom[AtStart]->rx;
			Vy = fAtom[AtStart]->ry;
		};
		KCycle = size - 2;
	};

	C = cos(2 * M_PI / size);
	S = sin(2 * M_PI / size);

	Na = fAtom.size();
	for(I = 0; I < KCycle; ++I) { // {cycle formation in structure}
		//{atom addition to structure}
		
		X1 = Vx - XCenter;
		Y1 = Vy - YCenter;
		X2 =  X1 * C + Y1 * S;
		Y2 = -X1 * S + Y1 * C;
		Vx = XCenter + X2;
		Vy = YCenter + Y2;
		SA = new TSingleAtom();
		SA->rx = Vx;
		SA->ry = Vy;
		fAtom.push_back(SA);
	};

	if(KCycle > 1)
		for(I = 0; I < KCycle-1; ++I) {
			fBond.push_back(new TSingleBond(Na + I, Na + I + 1));
		};

	if (cha1 < 0 && chb1 < 0) {
		//{if simple put-one more bond should be added}
		fBond.push_back(new TSingleBond(fAtom.size() - 1, Na));
	};

	if(cha1 >= 0) { 
		//{if connection through atom-two bonds should be added}
		fBond.push_back(new TSingleBond(cha1, Na));
		fBond.push_back(new TSingleBond(cha1, fAtom.size() - 1));
	};

	if(chb1 >= 0) { 
		//{if connection through atom-two bonds should be added (other attributes)}
		fBond.push_back(new TSingleBond(AtStart, Na));
		fBond.push_back(new TSingleBond(AtEnd, fAtom.size() - 1));

		//{Two or three bond connection is calculated}
		for(I = 0; I < nAt; ++I) {
			if(nAt < fAtom.size() && I < fAtom.size()) //  {Search for concided atoms}
			for(J = nAt; J < fAtom.size();) {
				RL = distance(I, J);
				Test = (I != J) && (RL <= 3.0);//!!!!!!!!!!!!!!!!!?????????????????
				if(Test) {     // {Concided pair of atoms found}
					for(K = 0; K < fBond.size(); ++K) {			// {Re-assignment of bonds}
						if(fBond[K]->at[0] == J) fBond[K]->at[0] = I;
						if(fBond[K]->at[1] == J) fBond[K]->at[1] = I;
					};
					deleteAtom(J);
					//{Search for bonds, which have identical pair of atom or which}
					//{make a loop. Such bonds can arise during re-assignment}
					for(K = 0; K <fBond.size() - 1; ++K) {
						if(K < fBond.size() - 1)
						for(L = K + 1; L < fBond.size();) {
							Test = (fBond[K]->at[0] == fBond[L]->at[0]) && (fBond[K]->at[1] == fBond[L]->at[1]);
							Test = Test || (fBond[K]->at[0] == fBond[L]->at[1]) && (fBond[K]->at[1] == fBond[L]->at[0]);
							Test = Test || (fBond[L]->at[0] == fBond[L]->at[1]);
							if(Test) {		//{Bad bonds was found-deletion}
								deleteBond(L);
							} else {
								++L;
							};
						};
					};
				} else { //end;
					++J;
				};
			};
		};
	};

	//DefineConn; {recalculation of attributes CONN}
	//{structure output}
	defineAtomConn();

	//!!!!!!!!!!!!!!!Fix it
	//if (not InsideBoundaries) and FIOPT8 then begin
	//    Rect:=ClientRect;
	//    FMolecule.RectResize(Rect.Right*0.15,Rect.Bottom*0.1,Rect.Right*0.85,Rect.Bottom*0.9);
	//  end;
	//  Invalidate;
	//end;
	//IsModified:=True;
};

//the procedure resize major structure (Natoms, NBonds, Atom, Bond, Conn
//see global variables description) to rect
void TEditedMolecule::rectResize(double XMin, double YMin, double XMax, double YMax)
{
	double			R, Rxmin, Rxmax, Rymin, Rymax, Rzmin, Sx, Sy, CorrX, CorrY;
	int				I;
	double			Xv, Yv;

	//Result:=FFontHeight;
	xShiftResize=0;
	yShiftResize=0;
	scaleResize=0;  //Is used to store RectResize values

	if(fAtom.empty()) return;
	if(XMin >= XMax) XMin = XMax-1;
	if(YMin >= YMax) YMin = YMax-1;
	if(fAtom.size() == 1) {
		fAtom[0]->rx = XMin + (XMax - XMin) / 2;
		fAtom[0]->ry = YMin + (YMax - YMin) / 2;
		xShiftResize=XMin;
		yShiftResize=XMax;
		scaleResize=1;  //Is used to store RectResize values

		return;
	}
	Rxmin = fAtom[0]->rx;
	Rxmax = fAtom[0]->rx;
	Rymin = fAtom[0]->ry;
	Rymax = fAtom[0]->ry;

	Sx = XMax - XMin;
	Sy = YMax - YMin;
	if(Sx == 0.0) Sx = 1;
	if(Sy == 0.0) Sy = 1;

	for(I=0; I < fAtom.size(); ++I) { //    {Minimal and maximal values}
		if(fAtom[I]->rx < Rxmin) Rxmin = fAtom[I]->rx;
		if(fAtom[I]->rx > Rxmax) Rxmax = fAtom[I]->rx;
		if(fAtom[I]->ry < Rymin) Rymin = fAtom[I]->ry;
		if(fAtom[I]->ry > Rymax) Rymax = fAtom[I]->ry;
	};


  CorrX = 0; 
  CorrY = 0;

  //{Determines, what is major-X scaling or Y scaling}
  if(((Rxmax - Rxmin) / Sx > (Rymax - Rymin)/Sy) || (YMax-YMin) == 0.0) {
    R = (Rxmax - Rxmin) / Sx;
    if(R < 1e-9) R = 1;
    CorrY = (YMax - YMin - (Rymax - Rymin) / R) / 2;
  } else {
    R = (Rymax - Rymin) / Sy;
    if(R < 1e-9) R = 1;
    CorrX = (XMax - XMin - (Rxmax - Rxmin) / R) / 2;
  };

  for(I = 0; I < fAtom.size(); ++I) {   //{Rescaling structure}
    fAtom[I]->rx = CorrX + (fAtom[I]->rx - Rxmin) / R + XMin;
    fAtom[I]->ry = CorrY + (fAtom[I]->ry - Rymin) / R + YMin;
  };
  xShiftResize=CorrX-Rxmin/R+XMin;
  yShiftResize=CorrY-Rymin/R+YMin;
  scaleResize=1/R;  //Is used to store RectResize values
};


void TEditedMolecule::findAllInsertions(TSimpleMolecule * molecule1, std::vector<std::vector<int> *> & all) {
//All contains a number of vectors for each insertion
//Each insertion vector has dimension N where N-number of atoms in INITIAL molecule, used for PrepareQuery routine
//N can be greater than number of atoms in current EditedMolecule if hydrogens are present
//Each component of vector contains corresponding atom number if molecule1. Negative value means hydrogen attached 
//to the inverse (negated) value of atom number in molecule1


	int cycleNumber;
	int j,k,l,m,mm;
	bool test;
	bool test1;
	bool test2,test3,stereoS;
	int aq1,aq2,as1,as2,j1,i1;
	int ii;
	bool result=false;
	bool whiletest1,whiletest2;
	std::vector<bool> unsaturated(nAtoms()+molecule1->nAtoms());
	//bool			aEQ[NATOMSMAX][NQUERYMAX];
	//bool			bEQ[NBONDSMAX][NQUERYMAX];
	std::vector<bool> aEQ(nAtoms() * molecule1->nAtoms() + 1);
	std::vector<bool> bEQ(nBonds() * molecule1->nBonds() + 1);

	std::vector<int> * tL;
	TSingleAtom * sAS;
	TSingleAtom * sAQ;


	for (j=0; j<all.size(); j++) delete(all[j]);
	all.clear();

    if ((! fIsQueryPrepare) && (nAtoms() == 1) && (molecule1)) {
      sAQ=getAtom(0);
	  for (j=0; j<molecule1->nAtoms(); j++) {
		sAS=molecule1->getAtom(j);
		if ((sAS->na == sAQ->na) && (sAS->nc == sAQ->nc) && (sAS->rl == sAQ->rl) && (sAS->nv == sAQ->nv)) {
          test=true;
		  if ((queryQHydr.size() > 0) && (queryQHydr[0] > molecule1->getNH(j))) test=false;
          if (test) {
            tL=new std::vector<int>;
		    tL->push_back(0);
            all.push_back(tL);
		  };  
	    };
	  };
      return;
	};


	if ((molecule1==NULL) || (! fIsQueryPrepare)) return;
	if (molecule1->nAtoms()==0) return;

	if (molecule1->listarSize()>aSTested.size()) aSTested.resize(molecule1->listarSize());
	if (molecule1->listarSize()>bSTested.size()) bSTested.resize(molecule1->listarSize());
	if (molecule1->listarSize()>queryQHydr.size()) queryQHydr.resize(molecule1->listarSize());
	if (molecule1->listarSize()>queryAQTested.size()) queryAQTested.resize(molecule1->listarSize());
	if (molecule1->listarSize()>queryBQTested.size()) queryBQTested.resize(molecule1->listarSize());
	if (molecule1->listarSize()>queryCurrentAssignment.size()) queryCurrentAssignment.resize(molecule1->listarSize());
	if (this->listarSize()>queryAQTested.size()) queryAQTested.resize(this->listarSize());
	if (this->listarSize()>queryBKNew.size()) queryBKNew.resize(this->listarSize());
	if (molecule1->listarSize()>structureBKNew.size()) structureBKNew.resize(molecule1->listarSize());

	//for (int i=0; i<molecule1->nAtoms(); i++) for (j=0; j<nAtoms(); j++) aEQ[i][j]=false;
	//for (int i=0; i<molecule1->nBonds(); i++) for (j=0; j<nBonds(); j++) bEQ[i][j]=false;
	for (int i = 0; i < aEQ.size(); i++) aEQ[i] = false;
	for (int i = 0; i < bEQ.size(); i++) bEQ[i] = false;

	int n=0;
	for (int i=0; i<queryEnum.size(); i++) if (queryEnum[i] == 0) {
	  if (n == 0) n++; else queryEnum[i]=-1;
	};

	cycleNumber=0;
	//Formula filter
	test2=true;   //turn of molecular formula filter

	molecule1->fIOPT11=fIOPT11;
	molecule1->fIOPT12=fIOPT13;
	molecule1->fIOPT13=fIOPT13;
	if (test2) {

		//{R/S/Z/E description are removed: they are not used in substructure search}
		if (fIOPT13) molecule1->atomBondChange(); //Semipolar bond conversion
		molecule1->defineAtomConn();
		molecule1->defineBondConn(structureBKNew);
		stereoS=molecule1->stereoBondChange(); //Stereo bond conversion}
		for (j=0; j<molecule1->nAtoms(); j++) {
			//Isotops conversion}
			if ((molecule1->getAtom(j)->na==104) && (! fIOPT11)) molecule1->getAtom(j)->na=1;
			if (! fIOPT11) molecule1->getAtom(j)->iz=0;
			bSTested[j]=0;  //bsTested contains nu. explicit hydrogens - so 0 is fine!
		};
		for (j=0; j<molecule1->nAtoms(); j++) if (molecule1->getAtom(j)->na==1)
			if (molecule1->getAtom(j)->nb>0) for (k=0; k<molecule1->getAtom(j)->nb; k++) {
				//Formation of list explicitly-defined hydrogens}
				l=molecule1->getAtom(j)->ac[k];
				bSTested[l]=bSTested[l]+1;
			};
		for (k=0; k<nAtoms(); k++) aSTested[k]=-1;
		if (nAtoms()==1) {
			//Partial case - if single atom was defined}
			test2=false;
			j=0;
			whiletest1=true;
			if (molecule1->nAtoms()>0) while (whiletest1) {
				test2=TSingleAtom::atomEquivalent(molecule1->getAtom(j),getAtom(0),bSTested[j],queryQHydr[0],fIOPT10,fIOPT11,
				compareHard,rejectAnyInStructure,false);
				//Addition for aromatic search}
				if (test2 && ((getAtom(0)->special & AROMATIC_MASK)!=0) ){
					test2=false;
					for (m=0; m<=structureBKNew[j].nb; m++) {
						mm=structureBKNew[j].adjusted[m];
						if ((molecule1->getBond(mm)->db==2) || (molecule1->getBond(mm)->db==3)) {
							test2=true;
							break;
						};
					};
				};
				//End addition}
				if (test2) {
					if (queryEnum.size()<1) queryEnum.resize(1);
					queryEnum[0]=0;
					if (queryAQTested.size()<1) queryAQTested.resize(1);
					queryAQTested[0]=j;
                    tL=new std::vector<int>;
		            tL->push_back(j);
                    all.push_back(tL);

				};
				j++;
				whiletest1=test2 || (j==molecule1->nAtoms());
				whiletest1=! whiletest1;
			};
		} else if (nBonds()==0) test2=false; else {
			//General case - substructure search
			if (fIncludedList != NULL) {
				for (j = 0; j < molecule1->nAtoms(); j++) for (k = 0; k < nAtoms(); k++) setEquivalentValue(k, j, nAtoms(), false, aEQ);  //aEQ[j][k]=false;
				for (k=0; k<nAtoms(); k++) {
					j=(*fIncludedList)[k];
					//aEQ[j][k]=true;
					setEquivalentValue(k, j, nAtoms(), true, aEQ);
					aSTested[k]=1;
				};
			} else {
                for (j=0; j<molecule1->nAtoms(); j++) unsaturated[j]=false;
                for (j=0; j<molecule1->nBonds(); j++) {
				  if ((molecule1->getBond(j)->tb == 2) || (molecule1->getBond(j)->tb == 3)) {
					unsaturated[molecule1->getBond(j)->at[0]]=true;
                    unsaturated[molecule1->getBond(j)->at[1]]=true;
			      };
		        };

				for (j=0; j<molecule1->nAtoms(); j++) for (k=0; k<nAtoms(); k++) {
					//Creation of atom-equivalent matrix}

					test3=TSingleAtom::atomEquivalent(molecule1->getAtom(j),getAtom(k),bSTested[j],
						queryQHydr[k],fIOPT10,fIOPT11,compareHard,rejectAnyInStructure,unsaturated[j]) && (getAtom(k)->nb<=molecule1->getAtom(j)->nb);

					//Addition for aromatic search
					if (test3 && ((getAtom(k)->special & AROMATIC_MASK) !=0)) {
						test3=false;
						for (m=0; m<structureBKNew[j].nb; m++) {
							mm=structureBKNew[j].adjusted[m];
							if ((molecule1->getBond(mm)->db==2) || (molecule1->getBond(mm)->db==3)) {
								test3=true;
								break;
							};
						};
					};
					//aEQ[j][k]=test3;
					setEquivalentValue(k, j, nAtoms(), test3, aEQ);
					if (test3) aSTested[k]=1;
				};
			};

			for (j=0; j<nBonds(); j++) queryAQTested[j]=0;
			for (j=0; j<molecule1->nBonds(); j++) for (k=0; k<nBonds(); k++) {
				//Creation of bond-equivalent matrix
				aq1=getBond(k)->at[0];   aq2=getBond(k)->at[1];
				as1=molecule1->getBond(j)->at[0]; as2=molecule1->getBond(j)->at[1];
				//test2=(aEQ[as1][aq1] && aEQ[as2][aq2]) || (aEQ[as1][aq2] && aEQ[as2][aq1]);
				test2 = (getEquivalentValue(aq1, as1, nAtoms(), aEQ) && getEquivalentValue(aq2, as2, nAtoms(), aEQ)) || (getEquivalentValue(aq2, as1, nAtoms(), aEQ) && getEquivalentValue(aq1, as2, nAtoms(), aEQ));

				if (test2) {
					test3=TSingleBond::bondEquivalent(molecule1->getBond(j),getBond(k),compareHard);
					//bEQ[j][k]=test3;
					setEquivalentValue(k, j, nBonds(), test3, bEQ);
					if (test3) queryAQTested[k]=1;
				}
				else setEquivalentValue(k,j,nBonds(),false,bEQ);  //bEQ[j][k]=false;
			};
			//Check, if all query atoms has a partner in structure
			test2=allQueryPresent(aSTested,queryAQTested,nAtoms(),nBonds());
			if (test2) {
				j1=0;

				for (k=0; k<molecule1->nAtoms(); k++) if (getEquivalentValue(0,k, nAtoms(), aEQ)) { //(aEQ[k][0]) {
					//Collection of structure atoms, which may be associated with first query atom
					j=0;
					whiletest1=true;
					while (whiletest1) {
						//test=bEQ[structureBKNew[k].adjusted[j]][0];
						test = getEquivalentValue(0, structureBKNew[k].adjusted[j], nBonds(), bEQ);
						if (test) {
							queryCurrentAssignment[j1]=k;
							j1++;
						};
						j++;
						whiletest1=test || (j==structureBKNew[k].nb);
						whiletest1=! whiletest1;
					};
				};
				i1=0;
				whiletest1=true;
				if (j1>0) while (whiletest1) { //for each structure, which may be assigned to 1-st query
					for (j=0; j<nAtoms(); j++) queryAQTested[j]=-1; //array initializing
					for (j=0; j<nBonds(); j++) queryBQCounter[j]=0;
					for (j=0; j<nBonds(); j++) queryBQTested[j]=-1;
					for (j=0; j<molecule1->nBonds(); j++) bSTested[j]=-1;
					for (j=0; j<molecule1->nAtoms(); j++) aSTested[j]=-1;
					queryAQTested[0]=queryCurrentAssignment[i1];
					aSTested[queryCurrentAssignment[i1]]=0;
					ii=0;
					whiletest2=true;
					while (whiletest2) { //start recursion
						directBondAss(ii,test,test1,bEQ,aEQ,queryBQCounter,queryAQTested,bSTested,
							queryBQTested,aSTested,queryAGer,queryBKNew,structureBKNew,molecule1,queryStas);
						if ((test) && (ii == (nBonds()-1))) {
                          tL=new std::vector<int>;
				 		  tL->resize(queryEnum.size(),0);
                          for (j=0; j<queryEnum.size(); j++) {  //Changes from 30 April 2006
                            //Here is NEW procedure, which is 1:1 copy QSMatch property
                            k=queryEnum[j];
                            if (k >= 0) k=queryInverse[k]; else if (k < 0) k=-(queryInverse[abs(k)]+1);
                            if ((k < 0) && (abs(k) <= queryAQTested.size())) {
                              //Search for H}
                              k=-(queryAQTested[abs(k)-1]+1);
							} else if ((k >= 0) && (k < queryAQTested.size())) k=queryAQTested[k];
                            //End NEW procedure - below - single-line of OLD procedure}
                            (*tL)[j]=k;
						  };
						  all.push_back(tL);
                          test=false;
                          //Addition from 6 May 2006
                          k=queryAGer[ii];
                          if (k >= 0) {
                            l=queryAQTested[k];
                            queryAQTested[k]=-1;
                            aSTested[l]=-1;
						  };
                //End addition

						};

						if ((! test) && (ii>=1)) {
							//previous atom were badly assigned-backstep
							queryBQCounter[ii]=0;
							k=queryAGer[ii-1];
							if (k>=0) {
								l=queryAQTested[k]; queryAQTested[k]=-1; aSTested[l]=-1;
							};
							bSTested[queryBQTested[ii-1]]=-1;
							queryBQTested[ii-1]=-1;
							ii=ii-2;
							test=true;
						};
						ii=ii+1; //query bond counter


						whiletest2=(ii==nBonds()) || ((! test) && (ii==1));
						whiletest2=! whiletest2;
					};
					i1++;
					test2=(ii==nBonds()); //Checking, if success has been reached
//					whiletest1=((j1==i1) || test2);
//					whiletest1=! whiletest1;
					whiletest1=(i1 < j1);
				} else test2=false;
			};
		};
	};

	/*
	if ((bondLabel != NULL) && (nBonds()>0)) {
		bondLabel->resize(molecule1->nBonds());
		if (test2) for (j=0; j<molecule1->nBonds(); j++) if (bSTested[j]>=0) (*bondLabel)[j]=1; else (*bondLabel)[j]=0;
	};
	*/
	result=test2;

};

void createSingleMolecule(TSimpleMolecule * largeMolecule, std::vector<TSimpleMolecule *> & molList) {
  double r,r1,xLeft,xMin,xMax;
  int i,j,n,nn,k;
  TSimpleMolecule * sm;
  TEditedMolecule em;
  double yMin,yMax;

  largeMolecule->clear();
  n=0; r=0;
  for (i=0; i<molList.size(); i++) {
    sm=molList[i];
    r1=sm->averageBondLength();
    if (r1>0) {
      r=r+r1;
      n++;
    };
  };
  if (n>0) {
	r=r/n; 
    for (i=0; i<molList.size(); i++) {
      sm=molList[i];
      r1=sm->averageBondLength();
	  if (r1>0) for (j=0; j<sm->nAtoms(); j++) {
		sm->getAtom(j)->rx=sm->getAtom(j)->rx*r/r1;
		sm->getAtom(j)->ry=sm->getAtom(j)->ry*r/r1;
	  };
	};
  } else r=1;
  xLeft=0;

  for (i=0; i<molList.size(); i++) {
    sm=molList[i];
      //Y scaling
    yMin=0; yMax=0;
    for (j=0; j<sm->nAtoms(); j++) {
      if ((j == 0) || (sm->getAtom(j)->ry < yMin)) yMin=sm->getAtom(j)->ry;
      if ((j == 0) || (sm->getAtom(j)->ry > yMax)) yMax=sm->getAtom(j)->ry;
    };
    for (j=0; j<sm->nAtoms(); j++) sm->getAtom(j)->ry=(sm->getAtom(j)->ry-yMin)-(yMax-yMin)/2;
    nn=sm->userData;
    if (nn == 0) nn=1;
    for (k=1; k<=nn; k++) {
      xMin=0; xMax=0;
      for (j=0; j<sm->nAtoms(); j++) {
        if ((j == 0) || (sm->getAtom(j)->rx < xMin)) xMin=sm->getAtom(j)->rx;
        if ((j == 0) || (sm->getAtom(j)->rx > xMax)) xMax=sm->getAtom(j)->rx;
	  };
      for (j=0; j<sm->nAtoms(); j++) {
        sm->getAtom(j)->rx=sm->getAtom(j)->rx+xLeft-xMin+r/1;
	  };
      em.addMoleculeNoScale(sm);
      for (j=0; j<sm->nAtoms(); j++) {
        sm->getAtom(j)->rx=sm->getAtom(j)->rx-xLeft+xMin+r/1;
	  };
      xLeft=xLeft+(xMax-xMin)+r+r/1;
	};
  };
  largeMolecule->moleculeCopy(em);
  largeMolecule->defineAtomConn();
  largeMolecule->allAboutCycles();
};

void decomposeSingleMoleculeNoCoordinates(TSimpleMolecule * largeMolecule, std::vector<TSimpleMolecule *> & molList, bool removeHydrogens) {
  TSimpleMolecule * sm;
  TEditedMolecule em;
  int w;
  std::vector<int> aList(largeMolecule->listarSize());
  std::vector<int> dList;
  bool test;
  int i,j,n;
  TSingleAtom * sa;

  for (i=0; i<molList.size(); i++) delete(molList.at(i));
  molList.clear();
 
  for (i=0; i<largeMolecule->nAtoms(); i++) {
	sa=largeMolecule->getAtom(i);
	sa->fragIndex=i+1;  //frag!!
	sa->enumerator=-1;
	sa->astereo=0;
  };

  for (i=0; i<largeMolecule->nAtoms(); i++) if (largeMolecule->getAtom(i)->astereo == 0) {
    largeMolecule->makeFragment(w,aList,i,0);
	if (w == largeMolecule->nAtoms()) {
	  sm=new TSimpleMolecule();
      sm->moleculeCopy(*largeMolecule);
      if (removeHydrogens) {
		j=sm->nAtoms();
        if (j > 1) sm->removeExplicitHydrogens();
		if (sm->nAtoms() != j) {
		  sm->calculateAllIndeces();
		};
	  };
      molList.push_back(sm);
      break;
	} else {
	  dList.clear();
	  for (j=0; j<w; j++)  dList.push_back(aList[j]);
      sm=new TSimpleMolecule();
      largeMolecule->extractFragmentFast(dList,*sm);
      if (removeHydrogens) {
        if (sm->nAtoms() > 1) sm->removeExplicitHydrogens();
	  };
      sm->calculateAllIndeces();;
      molList.push_back(sm);
	  for (j=0; j<w; j++) largeMolecule->getAtom(aList[j])->astereo=1;
	};
  };
  for (i=0; i<largeMolecule->nAtoms(); i++) largeMolecule->getAtom(i)->astereo=0;
}

} // namespace MolStruct



