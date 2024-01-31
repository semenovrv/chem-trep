#include "napi.chem.h"
#include "simple_molecule.h"
#include "svg_molecule.h"
#include "edited_molecule.h"
#include <sstream>
#include <unistd.h>

using namespace Napi;



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
	Napi::Env env = Env();
	if(result){
		Napi::Object rs = Napi::Object::New(env);
		rs.Set(Napi::String::New(env,"formula"),Napi::String::New(env,formula));
		rs.Set(Napi::String::New(env,"mw"),Napi::Number::New(env,mw));
		if(((Search_Mode*)mode)->mode==SVG_MODE){rs.Set(Napi::String::New(env,"svg"),Napi::String::New(env,res));}
		deferredPromise.Resolve(rs);
	}else{						deferredPromise.Resolve(env.Undefined());}}
virtual void OnError(const Napi::Error& e) override{finalize();deferredPromise.Reject(e.Value());}
void Execute() override{
std::istringstream mis(mol);
MolStruct::TSVGMolecule smol;
if(smol.readMolfile(mis)&&smol.nAtoms()>0){
	smol.defineAtomConn();
	smol.allAboutCycles();
	if(result=((MolStruct::TEditedMolecule*)query)->fragmentSearch(&smol,NULL)){
		formula=smol.getMolformula(true);
		mw=smol.getMolWeight(false);
		if(((Search_Mode*)mode)->mode==SVG_MODE){
			std::vector<std::string> data;
//			MolStruct::svgPolymerSave(smol,data,mol,((Search_Mode*)mode)->numerationDraw,((Search_Mode*)mode)->width,((Search_Mode*)mode)->height);
			//printf("cb_svgSave %s\n",((Search_Mode*)mode)->css.c_str());
			smol.cb_svgSave(data,((Search_Mode*)mode)->numerationDraw,((Search_Mode*)mode)->width,((Search_Mode*)mode)->height,"atoms",((Search_Mode*)mode)->css);
			std::string svg=data.size()?data[0]:"";
			for(size_t ii=1,i1=data.size();ii<i1;++ii){svg+="\n";svg+=data[ii];}
			res=svg;
}}}}
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
std::string formula;
double mw;
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
				if(info[4].IsBoolean()){						_mode.numerationDraw=info[4].As<Napi::Boolean>();
					if(info.Length()>5 && info[5].IsString()){	_mode.css=info[5].As<Napi::String>().Utf8Value();};
				}else if(info[4].IsString()){					_mode.css=info[4].As<Napi::String>().Utf8Value();
					if(info.Length()>5 && info[5].IsBoolean()){	_mode.numerationDraw=info[5].As<Napi::Boolean>();};
	}}};_mode.mode=md;};
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




Napi::Object Init(Napi::Env env, Napi::Object exports) {
    exports.Set(Napi::String::New(env, "QueryMol"),	QueryMol::GetClass(env));
	exports.Set(Napi::String::New(env, "SVG_MODE"),		Napi::Number::New(env,SVG_MODE));
	exports.Set(Napi::String::New(env, "BOOLEAN_MODE"),	Napi::Number::New(env,BOOLEAN_MODE));
	exports.Set(Napi::String::New(env, "SVG_CSS"),		Napi::String::New(env,SVG_CSS));

    return exports;
}

NODE_API_MODULE(addon, Init)

