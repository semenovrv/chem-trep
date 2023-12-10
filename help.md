# chem-trep

** Input formats [http://www.mdli.com/downloads/literature/ctfile.pdf](http://www.mdli.com/downloads/literature/ctfile.pdf) **
## Original C# methods: 
```js
const mol='any MDL molfile string'; // just a string with data
const svg_query=new QueryMol(mol,SVG_MODE[,svg_width,svg_height,svg_draw_mode]); // svg_draw_mode - any boolean value to draw hydrogens
```

## Getting promise object. The promise is fullfilled with the substructure search result as svg string in case of success, or undefined 
 
```js
svg_query.test(mol).then(svg=>svg?console.log('Success:',svg):console.log('no substructure match'),err=>console.error(err));
```

```js
const query=new QueryMol(mol)
```

## Getting promise object. The promise is fullfilled with the substructure search result as svg string in case of success, or undefined 
 
```js
query.test(mol).then(result=>result===BOOLEAN_MODE?console.log('Success'):console.log('no substructure match'),err=>console.error(err));
```



```js
const {MolSearch,SDFStreamBase,Query,SVG_MODE,BOOLEAN_MODE}= equire('chem-trep'); // Extended search methods contained in binding.js
```

### Prepare new query with svg results mode; mol - any MDL molfile 
```js
let search=new MolSearch(mol,{width:400,height:300}); 
```
### Search process initialization 
```js
search
	.on('success',(_mol,svg,id)=>svg?:console.log('query mol is not _mol substructure'){})
	.on('end',count=>{console.log('Search complete!',count);}) // count - total number of compared molecules
```
### A database is a stream of text MDL SDfile. SDFStreamBase starts search process while parsing the stream the same time. All parsed records are pushed to internal array. It is possible to use the same SDFStreamBase instance as many times as it needs

```js
const fs=require('fs')
const base=new SDFStreamBase(search,fs.createReadStream('db.sdf',{'encoding':'utf8'})
```


### It is possible to use the same SDFStreamBase instance as many times as it needs  */

```js
search=new MolSearch(mol1,{width:400,height:300}); // prepare query for the next search with the next query molecule
```

### Search process initialization 
```js
search
	.on('success',(_mol,svg,id)=>svg?:console.log('query mol is not _mol substructure'){})
	.on('end',count=>{console.log('Search complete!',count);}) // count - total number of compared molecules

search.run(base);
```
