#pragma once

#include <napi.h>
#include "edited_molecule.h"
#include "svg_molecule.h"

#define SVG_MODE 2
#define BOOLEAN_MODE 1
#define SVG_CSS "@import url(css/chem-trep.svg.mol.css);"

typedef struct{
    int mode=BOOLEAN_MODE;
	int width=SVG_IMAGE_WIDTH;
    int height=SVG_IMAGE_HEIGHT;
    bool numerationDraw=false;
    std::string css = SVG_INLINE_CSS;
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

