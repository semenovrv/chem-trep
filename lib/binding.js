const {QueryMol,SVG_MODE,BOOLEAN_MODE,SVG_CSS} = require('../build/Release/chem-trep')
	,os=require('os')
	,{EventEmitter}=require('events');

class CQuery extends QueryMol{
constructor(mol,mode){
	super(mol,...(mode?[SVG_MODE,mode.width,mode.height].concat(mode.css?[mode.css]:[]):[BOOLEAN_MODE]))
}
//readMolfile(mol)
//test(mol)
get['nAtoms'](){return super.nAtoms()}
get['nBonds'](){return super.nBonds()}
}

class MolSearch{
#qlist;
#race;
#threads;
#emitter;
#pos;
#next;
#speed;
constructor(mol,mode){
	const ppp=mode?[SVG_MODE,mode.width,mode.height].concat(mode.css?[mode.css]:[]):[BOOLEAN_MODE];
	this.#speed=os.cpus().length;
	this.#race=new Array(this.#speed);
	this.#qlist=new Array(this.#speed);
	for(let ii=0;ii<this.#race.length;ii++){this.#qlist[ii]=new QueryMol(mol,...ppp);this.#race[ii]=Promise.resolve(ii);}
	this.#threads=Array.from(this.#qlist);
}
['next'](mol){return new Promise((rs,rj)=>Promise.race(this.#race).then(id=>
	this.#race[id]=this.#race[id].then(()=>this.#qlist[id].test(mol).then(rs,rj).then(()=>id))
))}
['on'](...ppp){this.#emitter.on(...ppp); return this}
#testNext	(molbase,qq){const mol=molbase.get(this.#pos); 
	if(mol){	if(qq||(qq=this.#threads.pop()))this.#test(mol,qq).finally(()=>this.#next(molbase,qq));
	}else		if(qq)this.#threads.push(qq)
}
#emitNext	(molbase,qq){this.#emitter.emit('next',qq,molbase.get(this.#pos),this.#pos);}
#test	(mol,qq){const id=this.#pos++;return qq.test(mol).then(res=>res?this.#emitter.emit('success',mol,res,id):res,err=>console.error('test error',id));}
#run	(molbase){this.#next=this.#emitNext; process.nextTick(()=>{const q=this.#threads;while(q.length)this.#next(molbase,q.pop())});}
['snext'](molbase){this.#next(molbase)}
['emit'](...ppp){return this.#emitter.emit(...ppp)}
['run'](molbase){const search=this.#emitter=new EventEmitter();
	this.#pos=0;
	search.on('next',(qq,mol)=>{
		if(mol)this.#test(mol,qq).finally(()=>this.#next(molbase,qq));
		else if(this.#threads.push(qq)===this.#speed)search.emit('end',this.#pos);
	});
	if(molbase.stream){	this.#next=this.#testNext;
		search.once('parsed',len=>{	this.#run(molbase)})
	}else{							this.#run(molbase)}
	return this;
}
}



class SDFParser{
#re=		new RegExp('([\\s\\S]*?)(?=\\n>|$)([\\s\\S]*$)');
#re_trim=	new RegExp('\\n+$');
['split'](rec){const mm=rec.match(this.#re);
	return{'sdata':mm[1].replace(this.#re_trim,'')+'\n','data':mm[2]
		}}
}


class SBase{
#records;
constructor(search,stream,fn){
	(this.stream=stream).on('data',data=>{
		const 	ccc=data.replace(cr,'').split('$$$$\n');
		ccc[0]=ch0+ccc[0];
		ch0=ccc[ccc.length-1];
		for(let ii=0;ii<ccc.length-1;ii++){
			this.#records.push(ccc[ii]); //parser.split(ccc[ii]);
			search.snext(this);
		}
	}).on('end',()=>{console.log(`SDF[${stream}] parsed:`,this.length);delete this.stream;process.nextTick(()=>search.emit('parsed',this.length))});
	//	const parser=new SDFParser();
	const cr=new RegExp(/\r/,'gm');
	let ch0='';
	this.#records=[];
	this.svg=[];
	search.run(this);
	if(fn)fn(this,search);
}
['get'](i){return this.#records[i]}
get['length'](){return this.#records.length}
['send'](res){this.svg.push(res);return res}
};


class SDFStreamBase{
#records;
#emitter;
constructor(stream){
	stream.on('data',data=>{
		const 	ccc=data.replace(cr,'').split('$$$$\n');
		ccc[0]=ch0+ccc[0];
		ch0=ccc[ccc.length-1];
		for(let len,cc,ii=0;ii<ccc.length-1;ii++){len=this.#records.length;cc=ccc[ii]
			this.#records.push(cc);//parser.split(ccc[ii]);
			this.#emitter.emit('record',cc,len);
	}}).on('end',()=>{this.#emitter.emit('end',this.#records.length);console.log(`SDF[${stream}] parsed:`,this.#records.length);this.#emitter=undefined;});//	const parser=new SDFParser();
	const cr=new RegExp(/\r/,'gm');
	let ch0='';
	this.#records=[];
	this.#emitter=new EventEmitter();
	this.svg=[];
}
['each'](env,mol,ii){env.search.next(mol).then(res=>{if((++env.total===this.#records.length)&&!this.#emitter)env.resolve(env.total);res&&env.each(res)})}
['eachQueryResult'](search,fn){
const env={'total':0,'each':fn,'search':search},_each=this.each.bind(this,env);
return new Promise(rs=>{env.resolve=rs;
	this.#emitter?.on('record',_each)?.on('end',len=>(env.total===len)&&rs(len))
	this.#records.forEach(_each);
})}
['send'](res){this.svg.push(res);return res}
};






module.exports = {'MolSearch':MolSearch,'SDFStreamBase':SBase,'Query':CQuery,'SVG_MODE':SVG_MODE,'BOOLEAN_MODE':BOOLEAN_MODE,'SVG_CSS_FILENAME':SVG_CSS};
