var fs=require('fs');
var child_process=require('child_process');

var CLP=new CLParams(process.argv);

var alglist_path=CLP.unnamedParameters[0];
var datasetlist_path=CLP.unnamedParameters[1];
var basepath=CLP.namedParameters.basepath;
var outpath=CLP.namedParameters.outpath;

mkdir_safe(outpath);

var unit_numbers={
	h1:[1],
	m1:[2,3],
	m2:[2,3],
	m3:[2,3],
	m4:[2,3],
	m5:[2,3]
};

var algs=[];
{
	var txt=read_text_file(alglist_path);
	var lines=txt.split('\n');
	for (var i in lines) {
		var vals=lines[i].trim().split(' ');
		if (vals.length>=2) {
			algs.push({
				name:vals[0],
				script:vals[1],
				params:vals.slice(2).join(' ')
			});
		}
		else {
			if (lines[i].trim()) {
				throw 'problem in alglist file: '+lines[i].trim()
			}
		}
	}
}

var datasets=[];
{
	var txt=read_text_file(datasetlist_path);
	var lines=txt.split('\n');
	for (var i in lines) {
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

run_sorting(function() {
	compile_results(function() {
		process.exit(0);
	});
});

function compile_results(callback) {
	for (var d in datasets) {
		var dsname=datasets[d].name;
		var ks=unit_numbers[dsname];
		for (var ii in ks) {
			var k=ks[ii];
			for (var a in algs) {
				var algname=algs[a].name;
				var outpath0=outpath+'/'+algname+'-'+dsname;
				var CM=read_csv_matrix(outpath0+'/confusion_matrix.csv');
				var LM=read_csv_vector(outpath0+'/optimal_label_map.csv');

				var num=row_sum(CM,k-1);
				var num_fn,num_fp;
				if (LM[k-1]>0) {
					num_fn=num-CM[k-1][LM[k-1]-1];
					num_fp=col_sum(CM,LM[k-1]-1)-CM[k-1][LM[k-1]-1];
				}
				else {
					num_fn=num;
					num_fp=0;
				}

				print_csv_matrix(CM);
				console.log(algname+' '+dsname+' '+k+': num='+num+' fn='+num_fn+' fp='+num_fp);
			}
		}
	}
	callback();
}

function run_sorting(callback) {
	var num_running=0;
	for (var a in algs) {
		for (var d in datasets) {
			console.log('Applying '+algs[a].name+' to '+datasets[d].name);
			num_running++;
			apply_sorting(algs[a],datasets[d],function() {
				num_running--;
			});
		}
	}

	function on_timeout() {
		console.log('# algorithms running: '+num_running);
		if (num_running==0) {
			callback();
		}
		setTimeout(on_timeout,1000);
	}
	setTimeout(on_timeout,1000);
}

function apply_sorting(alg,ds,callback) {
	var outpath0=outpath+'/'+alg.name+'-'+ds.name;
	mkdir_safe(outpath0);

	copy_file_sync(basepath+'/datasets/'+ds.folder+'/firings_true.mda.prv',outpath0+'/firings_true.mda.prv');

	var cmd='mountainprocess';
	args=['queue-script',basepath+'/algorithms/'+alg.script];
	args.push('--raw='+basepath+'/datasets/'+ds.folder+'/raw.mda.prv');
	args.push(basepath+'/datasets/'+ds.folder+'/params.json');
	args.push('--outpath='+outpath0);
	var params0=alg.params.split(' ');
	for (var i in params0) {
		args.push(params0[i]);
	}
	console.log (cmd+' '+args.join(' '));
	
	run_process(cmd,args,function() {
		compute_confusion_matrix(outpath0,function() {
			callback();
		});
	});
}

function compute_confusion_matrix(output_path,callback) {
	var cmd='mountainprocess';
	var args=['run-process','merge_firings'];
	args.push('--firings1='+output_path+'/firings_true.mda.prv');
	args.push('--firings2='+output_path+'/firings.mda');
	args.push('--confusion_matrix='+output_path+'/confusion_matrix.csv');
	args.push('--optimal_label_map='+output_path+'/optimal_label_map.csv');
	args.push('--firings_merged='+output_path+'/firings_merged.mda');
	args.push('--max_matching_offset=4');

	run_process(cmd,args,function() {
		callback();
	});
}

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
	var data=fs.readFileSync(src);
	fs.writeFileSync(dst,data);
}

function run_process(cmd,args,callback) {
	console.log('Running '+cmd+' '+args.join(' '));
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
		console.log (data);
		all_stdout+=data;
	});
	pp.stderr.on('data',function(data) {
		console.log (data);
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
	return X[0];
}

function print_csv_matrix(X) {
	var txt='';
	for (var r=0; r<X.length; r++) {
		console.log(X[r].join(','));
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