#ifndef _HDR_SCALAR_
#define _HDR_SCALAR_

#include <vector>

#define SCALAR_MAX 50

typedef struct tagSingleScalar {
  short int nb;          //number of different pairs of bonds
  short int ref[SCALAR_MAX];     //bond's number, connected to first atom
  short int rel[SCALAR_MAX];     //bond's number, connected to second
  short int scalar[SCALAR_MAX];   //scalar product (0, 1, 2) between above bonds
} TSingleScalar;

class TScalar {
private:
  std::vector<TSingleScalar> data;
public:
	TScalar(void);
	TScalar(int itemCount);
	~TScalar(void);
    void          putData(int index, TSingleScalar value);
    TSingleScalar getData(int index);
    void          setNB(int index, short int value);
    short int     getNB(int index);
    void          setRef(int index1, int index2,short int value);
    short int     getRef(int index1, int index2);
    void          setRel(int index1, int index2, short int value);
    short int     getRel(int index1, int index2);
    void          setScalar(int index1, int index2, short int value);
    short int     getScalar(int index1, int index2);
	void          resize(int newSize);
};

#endif