#pragma once

#include <napi.h>
#include "edited_molecule.h"

#define SVG_MODE 2
#define BOOLEAN_MODE 1

typedef struct{
    int mode=BOOLEAN_MODE;
	int width=300;
    int height=150;
    bool numerationDraw=false;
} Search_Mode;

class QueryMol : public Napi::ObjectWrap<QueryMol>
{
public:
    QueryMol(const Napi::CallbackInfo&);
	Napi::Value nAtoms(const Napi::CallbackInfo&);
	Napi::Value nBonds(const Napi::CallbackInfo&);
	Napi::Value readMolfile(const Napi::CallbackInfo&);
	Napi::Value test(const Napi::CallbackInfo&);

    static Napi::Function GetClass(Napi::Env);

private:
	MolStruct::TEditedMolecule _edmol;
	Search_Mode _mode;
};

