#!/usr/bin/env node

var fs=require('fs');
var common=require(__dirname+'/boxsort_common.node.js');

var CLP=new common.CLParams(process.argv);

var algnames=CLP.unnamedParameters[0];
var dsnames=CLP.unnamedParameters[1];
var outpath=CLP.namedParameters.outpath;
var alglist_path=CLP.namedParameters.alglist;
var dslist_path=CLP.namedParameters.dslist;

if ((!algnames)||(!dsnames)||(!outpath)||(!alglist_path)||(!dslist_path)) {
	print_usage();
	process.exit(-1);
}

common.mkdir_safe(outpath);
var basepath=__dirname+'/../../../../mountainsort_experiments';

var algs=common.read_algs_from_text_file(alglist_path);
var datasets=common.read_datasets_from_text_file(dslist_path);

console.log(algs);
console.log(datasets);

for (var a in algs) {
	if (common.contains_alg(algnames,algs[a])) {
		for (var d in datasets) {
			if (common.contains_ds(dsnames,datasets[d])) {
				apply_sorting(algs[a],datasets[d]);
			}
		}
	}
}

common.wait_for_system_calls_to_finish(function() {

});

function print_usage() {
	console.log ('nodejs boxsort2.node.js [algnames|all] [dsnames|all] --outpath=[output] --alglist=[alglist.txt] --dslist=[dslist.txt]');
}

function apply_sorting(alg,ds,callback) {
	var outpath0=outpath+'/'+alg.name+'-'+ds.name;
	common.mkdir_safe(outpath0);

	common.copy_file_sync(basepath+'/datasets/'+ds.folder+'/firings_true.mda',outpath0+'/firings_true.mda');
	common.copy_file_sync(basepath+'/datasets/'+ds.folder+'/firings_true.mda.prv',outpath0+'/firings_true.mda.prv');

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
	common.make_system_call(cmd,args,function() {
		if (callback) callback();
	});
}


