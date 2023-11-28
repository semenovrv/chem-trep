#include "otrep.h"
#include "simple_molecule.h"
#include "svg_molecule.h"
#include "edited_molecule.h"
#include <sstream>
#include <unistd.h>

using namespace Napi;

Otrep::Otrep(const Napi::CallbackInfo& info) : ObjectWrap(info) {
    Napi::Env env = info.Env();

    if (info.Length() < 1) {
        Napi::TypeError::New(env, "Wrong number of arguments")
          .ThrowAsJavaScriptException();
        return;
    }

    if (!info[0].IsString()) {
        Napi::TypeError::New(env, "You need to name yourself")
          .ThrowAsJavaScriptException();
        return;
    }

    this->_greeterName = info[0].As<Napi::String>().Utf8Value();
}

Napi::Value Otrep::Greet(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();

    if (info.Length() < 1) {
        Napi::TypeError::New(env, "Wrong number of arguments")
          .ThrowAsJavaScriptException();
        return env.Null();
    }

    if (!info[0].IsString()) {
        Napi::TypeError::New(env, "You need to introduce yourself to greet")
          .ThrowAsJavaScriptException();
        return env.Null();
    }

    Napi::String name = info[0].As<Napi::String>();

    printf("Hello %s\n", name.Utf8Value().c_str());
    printf("I am %s\n", this->_greeterName.c_str());

    return Napi::String::New(env, this->_greeterName);
}

Napi::Function Otrep::GetClass(Napi::Env env) {
    return DefineClass(env, "Otrep", {Otrep::InstanceMethod("greet", &Otrep::Greet),});
}






class AsyncSMolTest : public Napi::AsyncWorker{
public:
static Napi::Value Create(const Napi::CallbackInfo& info,void* query,void* mode){
	AsyncSMolTest* worker = new AsyncSMolTest(info.Env(),info[0].As<Napi::String>().Utf8Value(),query,mode);
	worker->Queue();
	return worker->deferredPromise.Promise();
}
static Napi::Value Reject(Napi::Env env, const char* msg){
	Napi::Promise::Deferred failed=Napi::Promise::Deferred::New(env);
	failed.Reject(Napi::Error::New(env,msg).Value());
	return failed.Promise();
}
protected:
virtual void OnOK() override{
	finalize();
	if(result){
		if(((Search_Mode*)mode)->mode==SVG_MODE){deferredPromise.Resolve(Napi::String::New(Env(),res));}
		else{					deferredPromise.Resolve(Napi::Number::New(Env(),BOOLEAN_MODE));}
	}else{						deferredPromise.Resolve(Env().Undefined());}}
virtual void OnError(const Napi::Error& e) override{finalize();deferredPromise.Reject(e.Value());}
void Execute() override{
std::istringstream mis(mol);
MolStruct::TSimpleMolecule smol;
if(smol.readMolfile(mis)&&smol.nAtoms()>0){
	smol.defineAtomConn();
	smol.allAboutCycles();
/*	
	MolStruct::TSimpleMolecule sm;
	MolStruct::TEditedMolecule ed;
	sm = *(MolStruct::TEditedMolecule*)query;	// defineAtomConn,allAboutCycles done in operator=moleculeCopy
	ed.prepareQuery(sm); 						// defineAtomConn,allAboutCycles done in operator=moleculeCopy
*/	
	if((result=((MolStruct::TEditedMolecule*)query)->fragmentSearch(&smol,NULL))&&((Search_Mode*)mode)->mode==SVG_MODE){
		std::vector<std::string> data;
		MolStruct::svgPolymerSave(smol,data,((Search_Mode*)mode)->numerationDraw,((Search_Mode*)mode)->width,((Search_Mode*)mode)->height);
		std::string svg=data.size()?data[0]:"";
		for(size_t ii=1,i1=data.size();ii<i1;++ii){svg+="\n";svg+=data[ii];}
		res=svg;
}}}
private:
AsyncSMolTest(napi_env env,const std::string& _mol,void* _query,void* _mode) :
	Napi::AsyncWorker(env),
	mol(_mol),
	query(_query),
	mode(_mode),
	result(false),
	deferredPromise(Napi::Promise::Deferred::New(env)) {}

void finalize(){query=NULL;}


void* query;
void* mode;
std::string mol;
std::string res;
bool result;
Napi::Promise::Deferred deferredPromise;
};





QueryMol::QueryMol(const Napi::CallbackInfo& info) : ObjectWrap(info){Napi::Env env = info.Env();
	if(info.Length() < 1){				Napi::TypeError::New(env, "QueryMol: Wrong number of arguments").ThrowAsJavaScriptException();return;};
	if(!info[0].IsString()){			Napi::TypeError::New(env, "QueryMol: argument is not a string").ThrowAsJavaScriptException();return;};
    std::istringstream data(info[0].As<Napi::String>().Utf8Value());
	MolStruct::TSimpleMolecule smol;
	if(!smol.readMolfile(data)||smol.nAtoms()==0){	Napi::TypeError::New(env,"QueryMol: query is empty!").ThrowAsJavaScriptException();return;};
	if(info.Length()>1){
		if(!info[1].IsNumber()){		Napi::TypeError::New(env, "QueryMol: wrong mode type").ThrowAsJavaScriptException();return;};
		int md=info[1].As<Napi::Number>();
		if(md == SVG_MODE){
			if(info.Length() < 4){		Napi::TypeError::New(env, "QueryMol: Wrong number of mode arguments").ThrowAsJavaScriptException();return;};
			if(!info[2].IsNumber()){	Napi::TypeError::New(env, "QueryMol: svg width is not an integer").ThrowAsJavaScriptException();return;}
			else{_mode.width=info[2].As<Napi::Number>();};
			if(!info[3].IsNumber()){	Napi::TypeError::New(env, "QueryMol: svg height is not an integer").ThrowAsJavaScriptException();return;}
			else{_mode.height=info[3].As<Napi::Number>();};
			if(info.Length()>4){
				if(!info[4].IsBoolean()){Napi::TypeError::New(env, "QueryMol: svg draw mode is not a boolean").ThrowAsJavaScriptException();return;}
				else{_mode.numerationDraw=info[4].As<Napi::Boolean>();};
	}};_mode.mode=md;};
	smol.defineAtomConn();
	smol.allAboutCycles();
	_edmol.prepareQuery(smol);
}

Napi::Value QueryMol::nAtoms(const Napi::CallbackInfo& info){return Napi::Number::New(info.Env(),_edmol.nAtoms());}
Napi::Value QueryMol::nBonds(const Napi::CallbackInfo& info){return Napi::Number::New(info.Env(),_edmol.nBonds());}
Napi::Value QueryMol::readMolfile(const Napi::CallbackInfo& info){	Napi::Env env = info.Env();
	if(info.Length() < 1){	Napi::TypeError::New(env, "readMolfile: Wrong number of arguments").ThrowAsJavaScriptException();
		return env.Null();}
	if(!info[0].IsString()){Napi::TypeError::New(env, "readMolfile: argument is not a string").ThrowAsJavaScriptException();
		return env.Null();}
    std::istringstream data(info[0].As<Napi::String>().Utf8Value());
    _edmol.clear();
    return Napi::Boolean::New(env,_edmol.readMolfile(data));
}
Napi::Value QueryMol::test(const Napi::CallbackInfo& info){
	if(info.Length() < 1){	return AsyncSMolTest::Reject(info.Env(),"QueryMol.test: Wrong number of arguments");};
	if(!info[0].IsString()){return AsyncSMolTest::Reject(info.Env(),"QueryMol.test: argument is not a string");};
	return AsyncSMolTest::Create(info,&_edmol,&_mode);
}





Napi::Function QueryMol::GetClass(Napi::Env env) {
    return DefineClass(env, "QueryMol",{
			QueryMol::InstanceMethod("nAtoms",&QueryMol::nAtoms),
			QueryMol::InstanceMethod("nBonds",&QueryMol::nBonds),
			QueryMol::InstanceMethod("readMolfile",&QueryMol::readMolfile),
			QueryMol::InstanceMethod("test",&QueryMol::test),
			});
}






class AsyncSMolReader : public Napi::AsyncWorker {
public:
static Napi::Value Create(const Napi::CallbackInfo& info, void *vmol){
	if(		 info.Length()<1){		return Reject(info.Env(),"AsyncSMolReader: Wrong number of arguments");}
	else if(!info[0].IsString()){	return Reject(info.Env(),"AsyncSMolReader: argument is not a string");}
	std::string _mol = info[0].As<Napi::String>().Utf8Value();
	
	MolStruct::TSimpleMolecule* _smol = (MolStruct::TSimpleMolecule*)vmol;
	
	AsyncSMolReader* worker = new AsyncSMolReader(info.Env(),_smol,_mol);
	worker->Queue();
	return worker->deferredPromise.Promise();
}
protected:
static Napi::Value Reject(Napi::Env env, const char* msg){
	Napi::Promise::Deferred failed=Napi::Promise::Deferred::New(env);
	failed.Reject(Napi::Error::New(env,msg).Value());
	return failed.Promise();
}

void Execute() override{
	std::istringstream _mol(mol);
	smol->clear();
	result=smol->readMolfile(_mol);
	sleep(2);
    printf("Atoms/Bonds %d,%d\n",smol->nAtoms(),smol->nBonds());
}

virtual void OnOK() override{deferredPromise.Resolve(Napi::Boolean::New(Env(),result));}
virtual void OnError(const Napi::Error& e) override{deferredPromise.Reject(e.Value());}
private:
AsyncSMolReader(napi_env env, MolStruct::TSimpleMolecule *_smol, std::string& _mol) :
	Napi::AsyncWorker(env),
	mol(_mol),
	smol(_smol),
	result(false),
	deferredPromise(Napi::Promise::Deferred::New(env)) {}

std::string mol;
MolStruct::TSimpleMolecule* smol;
bool result;
Napi::Promise::Deferred deferredPromise;
};






SMol::SMol(const Napi::CallbackInfo& info) : ObjectWrap(info) {}
Napi::Value SMol::nAtoms(const Napi::CallbackInfo& info){return Napi::Number::New(info.Env(),_smol.nAtoms());}
Napi::Value SMol::nBonds(const Napi::CallbackInfo& info){return Napi::Number::New(info.Env(),_smol.nBonds());}
Napi::Value SMol::saveSVG(const Napi::CallbackInfo& info){Napi::Env env = info.Env();
    if (info.Length() < 3){		Napi::TypeError::New(env, "saveSVG: Wrong number of arguments").ThrowAsJavaScriptException();
        return env.Null();}
    if (!info[0].IsNumber()){	Napi::TypeError::New(env, "saveSVG: width is not an integer").ThrowAsJavaScriptException();
        return env.Null();}
    if (!info[1].IsNumber()){	Napi::TypeError::New(env, "saveSVG: height is not an integer").ThrowAsJavaScriptException();
        return env.Null();}
    if (!info[2].IsBoolean()){	Napi::TypeError::New(env, "saveSVG: draw mode is not a boolean").ThrowAsJavaScriptException();
        return env.Null();}

	std::vector<std::string> data;
	MolStruct::svgPolymerSave(_smol,data,info[2].As<Napi::Boolean>(),info[0].As<Napi::Number>(),info[1].As<Napi::Number>());
	std::string svg=data.size()?data[0]:"";
	for(size_t ii=1,i1=data.size();ii<i1;++ii){svg+="\n";svg+=data[ii];}
    return Napi::String::New(env,svg);
}

Napi::Value SMol::readMolfile(const Napi::CallbackInfo& info){Napi::Env env = info.Env();
	if (info.Length() < 1){		Napi::TypeError::New(env, "readMolfile: Wrong number of arguments").ThrowAsJavaScriptException();
        return env.Null();}
	if (!info[0].IsString()){	Napi::TypeError::New(env, "readMolfile: argument is not a string").ThrowAsJavaScriptException();
        return env.Null();}
    
    std::istringstream data(info[0].As<Napi::String>().Utf8Value());
	_smol.clear();
    return Napi::Boolean::New(env,_smol.readMolfile(data));
}

Napi::Value SMol::prepareQuery(const Napi::CallbackInfo& info){Napi::Env env = info.Env();
//readMolfile is completed here//
	if (info.Length() < 1){		Napi::TypeError::New(env, "prepareQuery: mode is undefined").ThrowAsJavaScriptException();
        return env.Null();}
	if (!info[0].IsNumber()){	Napi::TypeError::New(env, "prepareQuery: wrong mode type").ThrowAsJavaScriptException();
		return env.Null();}
	if(_smol.nAtoms()==0){		Napi::TypeError::New(env,"prepareQuery: query is empty!").ThrowAsJavaScriptException();
		return env.Null();}
	int md=info[0].As<Napi::Number>();
	if(md == SVG_MODE){
		if (info.Length() < 4){		Napi::TypeError::New(env, "prepareQuery: Wrong number of mode arguments").ThrowAsJavaScriptException();
		return env.Null();}
		if (!info[1].IsNumber()){	Napi::TypeError::New(env, "prepareQuery: svg width is not an integer").ThrowAsJavaScriptException();
		return env.Null();}else{_mode.width=info[1].As<Napi::Number>();}
		if (!info[2].IsNumber()){	Napi::TypeError::New(env, "prepareQuery: svg height is not an integer").ThrowAsJavaScriptException();
		return env.Null();}else{_mode.height=info[2].As<Napi::Number>();}
		if (!info[3].IsBoolean()){	Napi::TypeError::New(env, "prepareQuery: svg draw mode is not a boolean").ThrowAsJavaScriptException();
		return env.Null();}else{_mode.numerationDraw=info[3].As<Napi::Boolean>();}
	};_mode.mode=md;

	_smol.defineAtomConn();
	_smol.allAboutCycles();
	_query.prepareQuery(_smol);
}

Napi::Value SMol::match(const Napi::CallbackInfo& info){Napi::Env env = info.Env();
//prepareQuery is completed here//
//run directly after readMolfile only!!!!//
	_smol.defineAtomConn();
	_smol.allAboutCycles();
	if(_query.fragmentSearch(&_smol,NULL)){
		if(_mode.mode == SVG_MODE){
			std::vector<std::string> data;
			MolStruct::svgPolymerSave(_smol,data,_mode.numerationDraw,_mode.width,_mode.height);
			std::string svg=data.size()?data[0]:"";
			for(size_t ii=1,i1=data.size();ii<i1;++ii){svg+="\n";svg+=data[ii];}
				return Napi::String::New(env,svg);
		}else{	return Napi::Boolean::New(env,true);}
	}else{		return env.Null();}
}


Napi::Value SMol::readMolfilePromise(const Napi::CallbackInfo& info){return AsyncSMolReader::Create(info,&_smol);}

Napi::Function SMol::GetClass(Napi::Env env) {
    return DefineClass(env, "SMol",{
//			SMol::InstanceMethod("handle",&SMol::handle),
			SMol::InstanceMethod("nAtoms",&SMol::nAtoms),
			SMol::InstanceMethod("nBonds",&SMol::nBonds),
			SMol::InstanceMethod("readMolfile",&SMol::readMolfile),
			SMol::InstanceMethod("readMolfilePromise",&SMol::readMolfilePromise),
			SMol::InstanceMethod("prepareQuery",&SMol::prepareQuery),
			SMol::InstanceMethod("match",&SMol::match),
			SMol::InstanceMethod("saveSVG",&SMol::saveSVG),
			});
}



Napi::Object Init(Napi::Env env, Napi::Object exports) {
    exports.Set(Napi::String::New(env, "Otrep"),	Otrep::GetClass(env));
    exports.Set(Napi::String::New(env, "QueryMol"),	QueryMol::GetClass(env));
    exports.Set(Napi::String::New(env, "SMol"),			SMol::GetClass(env));
	exports.Set(Napi::String::New(env, "SVG_MODE"),		Napi::Number::New(env,SVG_MODE));
	exports.Set(Napi::String::New(env, "BOOLEAN_MODE"),	Napi::Number::New(env,BOOLEAN_MODE));

    return exports;
}

NODE_API_MODULE(addon, Init)

