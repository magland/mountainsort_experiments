var fs=require('fs');
var common=require(__dirname+'/boxsort_common.node.js');

//Use: npm install console.table --save
require('console.table');

var CLP=new common.CLParams(process.argv);

var algname=CLP.unnamedParameters[0];
var dsname=CLP.unnamedParameters[1];
var alglist_path=CLP.namedParameters.alglist;
var dslist_path=CLP.namedParameters.dslist;
var outpath=CLP.namedParameters.outpath;
var firings_name=CLP.namedParameters.firings_name||'firings.mda'; //or firings.annotated.mda
console.log(firings_name);

if ((!algname)||(!dsname)||(!firings_name)||(!alglist_path)||(!dslist_path)||(!outpath)) {
	print_usage();
	process.exit(-1);
}

var unit_numbers={};

var algs=common.read_algs_from_text_file(alglist_path);
var datasets=common.read_datasets_from_text_file(dslist_path);

compute_confusion_matrices(function() {
	tabulate_results();
});

function print_usage() {
	console.log('nodejs tabulate_results.node.js [algname] [dsname] --firings_name=[firings.mda] --outpath=[output] --alglist=[alglist.txt] --dslist=[dslist.txt]');
}

function compute_confusion_matrices(callback) {
	for (var a in algs) {
		var algname0=algs[a].name;
		console.log(algname0);
		if ((algname=='all')||(algname==algname0)) {
			for (var d in datasets) {
				var dsname0=datasets[d].name;
				if ((dsname=='all')||(dsname==dsname0)) {
					var outpath0=outpath+'/'+algname0+'-'+dsname0;
					compute_confusion_matrix(outpath0,firings_name,function() {});
				}
			}
		}
	}
	common.wait_for_system_calls_to_finish(callback);
}

function tabulate_results(callback) {
	for (var d in datasets) {
		var dsname0=datasets[d].name;
		if ((dsname=='all')||(dsname==dsname0)) {
			console.log ('');
			console.log ('######## DATASET: '+dsname0);
			console.log ('');
			var table0=[];
			for (var a in algs) {
				var algname0=algs[a].name;
				if ((algname=='all')||(algname==algname0)) {
					var outpath0=outpath+'/'+algname0+'-'+dsname0;
					var CM=common.read_csv_matrix(outpath0+'/confusion_matrix.csv');
					var LM=common.read_csv_vector(outpath0+'/optimal_label_map.csv');

					//var Ntrue=0,Ndetect=0,Ncorrect=0,Nincorrect=0,Nmissed=0,Nextra=0;
					var K1=CM.length-1;
					var K2=CM[0].length-1;
					
					var has_match=[];
					{
						for (var k2=1; k2<=K2; k2++) {
							has_match.push(0);
						}
						for (var k1=1; k1<=K1; k1++) {
							var k2_match=LM[k1-1];
							if (k2_match) {
								has_match[k2_match-1]=0;
							}
						}
					}

					var Ntrue=0,Ndetect=0,Ncorrect=0,Nincorrect=0,Nextra=0,Nmissed=0;
					for (var k1=1; k1<=K1; k1++) {
						for (var k2=1; k2<=K2; k2++) {
							var val=CM[k1-1][k2-1];
							Ntrue+=val;
							Ndetect+=val;
							if (LM[k1-1]==k2) {
								Ncorrect+=val;
							}
							else {
								//classified with a label that corresponds to a different true label
								Nincorrect+=val;
							}
						}
					}
					for (var k1=1; k1<=K1; k1++) {
						var val=CM[k1-1][K2];
						//true events that are not detected
						Ntrue+=val;
						Nmissed+=val;
					}
					for (var k2=1; k2<=K2; k2++) {
						var val=CM[K1][k2-1];
						//detected events that don't correspond to any true label
						Ndetect+=val;
						Nextra+=val;
					}

					var k_display='all';
					var table_row={
						DATASET:dsname0,
						ALG:algname0,
						UNIT:k_display,
						Ntrue:Ntrue,
						Ndetect:Ndetect,
						Ktrue:K1,
						Kdetect:K2,
						Correct:common.topct(Ncorrect/Ndetect),
						Incorrect:common.topct(Nincorrect/Ndetect),
						Extra:common.topct(Nextra/Ndetect),
						Missed:common.topct(Nmissed/Ntrue)
					};
					table0.push(table_row);
				}
			}
			console.table(table0);
		}
	}
	
	if (callback) callback();
}

function compute_confusion_matrix(output_path,callback) {
	var cmd='mountainprocess';
	var args=['run-process','merge_firings'];
	if (fs.existsSync(output_path+"/firings_true.mda.prv"))
		args.push('--firings1='+output_path+'/firings_true.mda.prv');
	else if (fs.existsSync(output_path+"/firings_true.mda"))
		args.push('--firings1='+output_path+'/firings_true.mda');
	else {
		if (callback) callback(); //no ground truth
		return;
	}
	args.push('--firings2='+output_path+'/'+firings_name);
	args.push('--confusion_matrix='+output_path+'/confusion_matrix.csv');
	args.push('--optimal_label_map='+output_path+'/optimal_label_map.csv');
	args.push('--firings_merged='+output_path+'/firings_merged.mda');
	args.push('--max_matching_offset=10');

	common.make_system_call(cmd,args,function() {
		if (callback) callback();
	});
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
		CMb[i1][K2b]=common.row_sum(CM,ks[i1]-1)-common.row_sum(CMb,i1);
	}
	for (var i2=0; i2<K2b; i2++) {
		CMb[K1b][i2]=common.col_sum(CM,ks2[i2]-1)-common.col_sum(CMb,i2);
	}
	return {
		CMb:CMb,
		LMb:LMb
	}
}

function compute_confusion_matrix(output_path,firings_name,callback) {
	var cmd='mountainprocess';
	var args=['run-process','merge_firings'];
	if (fs.existsSync(output_path+"/firings_true.mda.prv"))
		args.push('--firings1='+output_path+'/firings_true.mda.prv');
	else if (fs.existsSync(output_path+"/firings_true.mda"))
		args.push('--firings1='+output_path+'/firings_true.mda');
	else {
		callback(); //no ground truth file
		return;
	}
	args.push('--firings2='+output_path+'/'+firings_name);
	args.push('--confusion_matrix='+output_path+'/confusion_matrix.csv');
	args.push('--optimal_label_map='+output_path+'/optimal_label_map.csv');
	args.push('--firings_merged='+output_path+'/firings_merged.mda');
	args.push('--max_matching_offset=10');

	common.make_system_call(cmd,args,function() {
		callback();
	});
}