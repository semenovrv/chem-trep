

#include <cstring>

#include "scalar.h"

TScalar::TScalar(void)
{
}

TScalar::TScalar(int itemCount){
  data.resize(itemCount);
};


TScalar::~TScalar(void)
{
}

void TScalar::putData(int index, TSingleScalar value){
  if ((index >= 0) && (index < data.size())) data[index]=value;
};
    
TSingleScalar TScalar::getData(int index){
  TSingleScalar result;
  memset(&result,0,sizeof(result));
  if ((index >= 0) && (index < data.size())) result=data[index];
  return result;
};

void TScalar::setNB(int index, short int value){
  TSingleScalar result;
  if ((index >= 0) && (index < data.size())) {
	result=data[index];
	result.nb=value;
	data[index]=result;
  };
};

short int TScalar::getNB(int index){
  TSingleScalar vdata;
  short int result=0;

  if ((index >= 0) && (index < data.size())) {
	vdata=data[index];
	result=vdata.nb;
  }
  return result;
};

void TScalar::setRef(int index1, int index2, short int value){
  TSingleScalar vdata;

  if ((index1 >= 0) && (index1 < data.size())) {
	vdata=data[index1];
	if ((index2 >= 0) && (index2 < SCALAR_MAX)) {
	  vdata.ref[index2]=value;
	  data[index1]=vdata;
	};
  }
};

short int TScalar::getRef(int index1, int index2){
  TSingleScalar vdata;
  short int result=-1;

  if ((index1 >= 0) && (index1 < data.size())) {
	vdata=data[index1];
	if ((index2 >= 0) && (index2 < SCALAR_MAX)) result=vdata.ref[index2];
  }
  return result;
};

void TScalar::setRel(int index1, int index2, short int value){
  TSingleScalar vdata;

  if ((index1 >= 0) && (index1 < data.size())) {
	vdata=data[index1];
	if ((index2 >= 0) && (index2 < SCALAR_MAX)) {
	  vdata.rel[index2]=value;
	  data[index1]=vdata;
	};
  }
};

short int TScalar::getRel(int index1, int index2){
  TSingleScalar vdata;
  short int result=-1;

  if ((index1 >= 0) && (index1 < data.size())) {
	vdata=data[index1];
	if ((index2 >= 0) && (index2 < SCALAR_MAX)) result=vdata.rel[index2];
  }
  return result;
};
 
void TScalar::setScalar(int index1, int index2, short int value){
  TSingleScalar vdata;

  if ((index1 >= 0) && (index1 < data.size())) {
	vdata=data[index1];
	if ((index2 >= 0) && (index2 < SCALAR_MAX)) {
	  vdata.scalar[index2]=value;
	  data[index1]=vdata;
	};
  }
};
 
short int TScalar::getScalar(int index1, int index2){
  TSingleScalar vdata;
  short int result=-1;

  if ((index1 >= 0) && (index1 < data.size())) {
	vdata=data[index1];
	if ((index2 >= 0) && (index2 < SCALAR_MAX)) result=vdata.scalar[index2];
  }
  return result;
};

void TScalar::resize(int newSize){
  TSingleScalar vdata;

  memset(&vdata,0,sizeof(TSingleScalar));
  data.resize(newSize,vdata);
};
