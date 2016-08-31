var fs=require('fs');
var child_process=require('child_process');

var CLP=new CLParams(process.argv);

var basepath=CLP.namedParameters.basepath;
if (!basepath) {
	console.log('No basepath specified');
	return;
}
var dspath=basepath+'/datasets';

var dspath0=dspath+'/synth1';

mkdir_safe(dspath0);

cmd='mountainprocess';
args=['run-process','synthesize_timeseries_001_matlab'];
args.push('--waveforms='+dspath0+'/waveforms.mda');
args.push('--timeseries='+dspath0+'/raw.mda');
args.push('--firings_true='+dspath0+'/firings_true.mda');

params0={samplerate:30000};
write_text_file(dspath0+'/params.json',JSON.stringify(params0));

make_system_call(cmd,args,function() {
	process.exit(0);
});

function make_system_call(cmd,args,callback) {
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

function mkdir_safe(path) {
	try {
		fs.mkdirSync(path);
	}
	catch (err) {

	}
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

function read_text_file(path) {
	return fs.readFileSync(path,'utf8');
}

function write_text_file(path,txt) {
	fs.writeFileSync(path,txt,'utf8');
}
