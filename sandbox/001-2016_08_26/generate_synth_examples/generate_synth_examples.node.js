var fs=require('fs');
var child_process=require('child_process');

var CLP=new CLParams(process.argv);

var basepath=CLP.namedParameters.basepath;
if (!basepath) {
	console.log ('No basepath specified');
	return;
}
var dspath=basepath+'/datasets';

mkdir_safe(dspath+'/synth_examples');;
var dspath0=dspath+'/synth_examples/'+(CLP.namedParameters['dsname']||'default');
mkdir_safe(dspath0);

cmd='mountainprocess';
args='run-process synthesize_timeseries_001_matlab';
args+=' --waveforms='+dspath0+'/waveforms.mda';
args+=' --timeseries='+dspath0+'/raw.mda';
args+=' --firings_true='+dspath0+'/firings_true.mda';

var M0=CLP.namedParameters["M"]||5;
var K0=CLP.namedParameters["K"]||20;
var duration0=CLP.namedParameters["duration"]||600;
var amp_variation_min=CLP.namedParameters["amp_variation_min"]||1;
var amp_variation_max=CLP.namedParameters["amp_variation_max"]||1;
var firing_rate_min=CLP.namedParameters["firing_rate_min"]||1;
var firing_rate_max=CLP.namedParameters["firing_rate_max"]||1;
var noise_level=CLP.namedParameters["noise_level"]||1;

args+=" --M="+M0+" --T=800 --K="+K0+" --duration="+duration0+" --amp_variation_min="+amp_variation_min+" --amp_variation_max="+amp_variation_max;
args+=" --firing_rate_min="+firing_rate_min+" --firing_rate_max="+firing_rate_max;
args+=" --noise_level="+noise_level;

if ('_force_run' in CLP.namedParameters) args+=" --_force_run"

args=args.split(' ');

params0={samplerate:30000};
write_text_file(dspath0+'/params.json',JSON.stringify(params0));

make_system_call(cmd,args,function() {
	process.exit(0);
});

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
