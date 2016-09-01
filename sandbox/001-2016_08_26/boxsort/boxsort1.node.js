var fs=require('fs');
var child_process=require('child_process');

var CLP=new CLParams(process.argv);

var alglist_path=CLP.unnamedParameters[0];
var datasetlist_path=CLP.unnamedParameters[1];
var outpath=CLP.namedParameters.outpath;

if ((!alglist_path)||(!datasetlist_path)||(!outpath)) {
	print_usage();
	process.exit(-1);
}

mkdir_safe(outpath);
var basepath=__dirname+'/../../../../mountainsort_experiments';

var algs=read_algs_from_text_file(alglist_path);
var datasets=read_datasets_from_text_file(datasetlist_path);

run_sorting();


function run_sorting(callback) {
	console.log(algs);
	console.log(datasets);
	for (var a in algs) {
		for (var d in datasets) {
			apply_sorting(algs[a],datasets[d]);
		}
	}
}

function apply_sorting(alg,ds) {
	var outpath0=outpath+'/'+alg.name+'-'+ds.name;
	mkdir_safe(outpath0);

	copy_file_sync(basepath+'/datasets/'+ds.folder+'/firings_true.mda',outpath0+'/firings_true.mda');
	copy_file_sync(basepath+'/datasets/'+ds.folder+'/firings_true.mda.prv',outpath0+'/firings_true.mda.prv');

	var ds_folder=basepath+'/datasets/'+ds.folder;

	var cmd='mountainprocess';
	args='queue-script '+basepath+'/algorithms/'+alg.script;
	args+=' --raw='+ds_folder+'/raw.mda.prv';
	args+=' '+ds_folder+'/params.json';
	var geom_fname=ds_folder+'/geom.csv';
	if (fs.existsSync(geom_fname)) {
		args+=' --geom='+geom_fname;
	}
	args+=' --outpath='+outpath0;
	var arguments0=alg.arguments.split(' ');
	for (var i in arguments0) {
		args+=' '+arguments0[i];
	}
	args=args.split(' ');
	make_system_call(cmd,args,function() {});
}

////////////////////////////////////////////////////////////////////////////////////////////////

function mkdir_safe(path) {
	try {
		fs.mkdirSync(path);
	}
	catch (err) {

	}
}

function read_text_file(path) {
	return fs.readFileSync(path,'utf8');
}

function CLParams(argv) {
	this.unnamedParameters=[];
	this.namedParameters={};

	var args=argv.slice(2);
	for (var i=0; i<args.length; i++) {
		var arg0=args[i];
		if (arg0.indexOf('--')==0) {
			arg0=arg0.slice(2);
			var ind=arg0.indexOf('=');
			if (ind>=0) {
				this.namedParameters[arg0.slice(0,ind)]=arg0.slice(ind+1);
			}
			else {
				this.namedParameters[arg0]=args[i+1]||'';
				i++;
			}
		}
		else {
			this.unnamedParameters.push(arg0);
		}
	}
}

function copy_file_sync(src,dst) {
	if (!fs.existsSync(src)) return;
	var data=fs.readFileSync(src);
	fs.writeFileSync(dst,data);
}

function make_system_call(cmd,args,callback) {
	console.log ('Running '+cmd+' '+args.join(' '));
	var pp=child_process.spawn(cmd,args);
	pp.stdout.setEncoding('utf8');
	pp.stderr.setEncoding('utf8');
	var done=false;
	pp.on('close', function(code) {
  		done=true;
		callback();
	});
	pp.on('error',function(err) {
		console.log ('Process error: '+cmd+' '+args.join(' '));
		console.log (err);
	});
	var all_stdout='';
	var all_stderr='';
	pp.stdout.on('data',function(data) {
		console.log ('----'+data);
		all_stdout+=data;
	});
	pp.stderr.on('data',function(data) {
		console.log ('===='+data);
		all_stderr+=data;
	});
}

function transpose_matrix(X) {
	if (X.length==0) return X;
	var Y=[];
	for (var i in X[0]) {
		Y.push([]);
	}
	for (var j in X) {
		for (var i in X[j]) {
			Y[i].push(X[j][i]);
		}
	}
	return Y;
}

function sub_confusion_matrix(CM,LM,ks) {
	var ks2=[];
	var LMb=[];
	for (var i in ks) {
		var match=LM[ks[i]-1];
		if (match) {
			ks2.push(match);
			LMb.push(ks2.length);
		}
		else {
			LMb.push(0);
		}
	}
	var K1b=ks.length;
	var K2b=ks2.length;
	var CMb=[];

	for (var i1=0; i1<K1b; i1++) {
		var tmp=[];
		for (var i2=0; i2<K2b; i2++) {
			tmp.push(CM[ks[i1]-1][ks2[i2]-1]);
		}
		tmp.push(0); //place-holder
		CMb.push(tmp);
	}
	//last row
	{
		var tmp=[];
		for (var i2=0; i2<K2b; i2++) {
			tmp.push(0); //place-holder
		}
		tmp.push(0); //place-holder
		CMb.push(tmp);
	}
	for (var i1=0; i1<K1b; i1++) {
		CMb[i1][K2b]=row_sum(CM,ks[i1]-1)-row_sum(CMb,i1);
	}
	for (var i2=0; i2<K2b; i2++) {
		CMb[K1b][i2]=col_sum(CM,ks2[i2]-1)-col_sum(CMb,i2);
	}
	return {
		CMb:CMb,
		LMb:LMb
	}
}

function read_algs_from_text_file(file_path) {
	var algs=[];
	{
		var txt=read_text_file(file_path);
		var lines=txt.split('\n');
		for (var i in lines) {
			if (lines[i].trim().slice(0,1)!='#') {
				var vals=lines[i].trim().split(' ');
				if (vals.length>=2) {
					algs.push({
						name:vals[0],
						script:vals[1],
						arguments:vals.slice(2).join(' ')
					});
				}
				else {
					if (lines[i].trim()) {
						throw 'problem in alglist file: '+lines[i].trim();
					}
				}
			}
		}
	}
	return algs;
}

function read_datasets_from_text_file(file_path) {
	var datasets=[];
	{
		var txt=read_text_file(file_path);
		var lines=txt.split('\n');
		for (var i in lines) {
			if (lines[i].trim().slice(0,1)!='#') {
				var vals=lines[i].trim().split(' ');
				if (vals.length==2) {
					datasets.push({
						name:vals[0],
						folder:vals[1]
					});
				}
				else {
					if (lines[i].trim()) {
						throw 'problem in datasetlist file: '+lines[i].trim();
					}
				}
			}
		}
	}
	return datasets;
}

function read_csv_matrix(path) {
	var ret=[];
	var txt=read_text_file(path);
	var lines=txt.split('\n');
	for (var i in lines) {
		var vals=lines[i].split(',');
		if (vals.length>0) {
			var row=[];
			for (var k=0; k<vals.length; k++) {
				row.push(Number(vals[k]));
			}
			ret.push(row);
		}
	}
	return transpose_matrix(ret); //this is because of a bad decision I made
}

function read_csv_vector(path) {
	var X=read_csv_matrix(path);
	var Y=[];
	for (var i in X) {
		for (var j in X[i]) {
			Y.push(X[i][j]);
		}
	}
	return Y;
}

function print_csv_matrix(X) {
	var txt='';
	for (var r=0; r<X.length; r++) {
		console.log (X[r].join(','));
	}
}

function row_sum(X,row) {
	var ret=0;
	for (var i in X[row]) {
		ret=ret+X[row][i];
	}
	return ret;
}
function col_sum(X,col) {
	var ret=0;
	for (var i in X) {
		ret=ret+X[i][col];
	}
	return ret;
}

function topct(num) {
	if (isNaN(num)) return '';
	if (num>1) return '>100%';
	return Math.floor(num*100)+'%';
}

function clone(X) {
	return JSON.parse(JSON.stringify(X));
}