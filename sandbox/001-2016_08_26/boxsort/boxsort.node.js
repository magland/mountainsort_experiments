var fs=require('fs');
var child_process=require('child_process');

var CLP=new CLParams(process.argv);

var alglist_path=CLP.unnamedParameters[0];
var datasetlist_path=CLP.unnamedParameters[1];
var basepath=CLP.namedParameters.basepath;
var outpath=CLP.namedParameters.outpath;

mkdir_safe(outpath);

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
	process.exit(0);
});

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
	var args=['run-process','confusion_matrix'];
	args.push('--firings1='+output_path+'/firings.mda');
	args.push('--firings2='+output_path+'/firings_true.mda.prv');
	args.push('--output='+output_path+'/confusion_matrix.csv');
	args.push('--optimal_assignments='+output_path+'/optimal_assignments.csv');
	args.push('--event_correspondence='+output_path+'/event_correspondence.mda');
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
