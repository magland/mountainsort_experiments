var fs=require('fs');
var child_process=require('child_process');

var args=process.argv.slice(2);

var alglist_path=args[0];
var datasetlist_path=args[1];

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
		exit(0);
	}
	setTimeout(on_timeout,1000);
}
setTimeout(on_timeout,1000);

function apply_sorting(alg,ds,callback) {
	var basepath='/home/magland/dev/mountainsort_experiments';
	var outpath='/home/magland/tmp/'+alg.name+'-'+ds.name;
	var cmd='mountainprocess';
	args=['queue-script',basepath+'/algorithms/'+alg.script];
	args.push('--raw='+basepath+'/datasets/'+ds.folder+'/raw.mda.prv');
	args.push(basepath+'/datasets/'+ds.folder+'/params.json');
	args.push('--outpath='+outpath);
	var params0=alg.params.split(' ');
	for (var i in params0) {
		args.push(params0[i]);
	}
	console.log(cmd+' '+args.join(' '));
	mkdir_safe(outpath);
	console.log('test 1');
	var pp=child_process.spawn(cmd,args);
	console.log('test 2');
	var done=false;
	pp.on('exit', function(code) {
  		console.log('Process finished: '+outpath+' with code '+code);
  		done=true;
		callback();
	});
	pp.on('error',function(err) {
		console.log('Process error: '+outpath);
		console.log(err);
		done=true;
		callback();
	});
	var all_stdout='';
	var all_stderr='';
	pp.stdout.on('data',function(data) {
		all_stdout+=data;
		console.log('len of stdout: '+all_stdout.length);
	});
	pp.stderr.on('data',function(data) {
		all_stderr+=data;
		console.log('len of stderr: '+all_stderr.length);
	});
	function debug1() {
		console.log(':::: '+all_stdout.slice(all_stdout.slice-30));
		console.log('++++ '+all_stderr.slice(all_stderr.slice-30));
		if (!done) {
			setTimeout(debug1,3000);
		}
	}
	setTimeout(debug1,3000);
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