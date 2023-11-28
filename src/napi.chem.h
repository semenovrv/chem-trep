#pragma once

#include <napi.h>
#include "simple_molecule.h"
#include "edited_molecule.h"

#define SVG_MODE 2
#define BOOLEAN_MODE 1

typedef struct{
    int mode=BOOLEAN_MODE;
	int width=300;
    int height=150;
    bool numerationDraw=false;
} Search_Mode;

class Otrep : public Napi::ObjectWrap<Otrep>
{
public:
    Otrep(const Napi::CallbackInfo&);
    Napi::Value Greet(const Napi::CallbackInfo&);

    static Napi::Function GetClass(Napi::Env);

private:
    std::string _greeterName;
};

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

class SMol : public Napi::ObjectWrap<SMol>
{
public:
    SMol(const Napi::CallbackInfo&);
//	Napi::Value handle(const Napi::CallbackInfo&);
	Napi::Value nAtoms(const Napi::CallbackInfo&);
	Napi::Value nBonds(const Napi::CallbackInfo&);
	Napi::Value readMolfile(const Napi::CallbackInfo&);
	Napi::Value readMolfilePromise(const Napi::CallbackInfo&);
	Napi::Value prepareQuery(const Napi::CallbackInfo&);
	Napi::Value match(const Napi::CallbackInfo&);
	Napi::Value saveSVG(const Napi::CallbackInfo&);

    static Napi::Function GetClass(Napi::Env);

private:
	Search_Mode _mode;
	MolStruct::TSimpleMolecule _smol;
	MolStruct::TEditedMolecule _query;
};
