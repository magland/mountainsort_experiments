var fs=require('fs');

//Use: npm install console.table --save
require('console.table');

var CLP=new CLParams(process.argv);

var alglist_path=CLP.unnamedParameters[0];
var datasetlist_path=CLP.unnamedParameters[2];
var outpath=CLP.namedParameters.outpath;

if ((!alglist_path)||(!datasetlist_path)||(!outpath)) {
	print_usage();
	process.exit(-1);
}

//finish!

var basepath=__dirname+'/../../../../mountainsort_experiments';


mkdir_safe(outpath);

var alg_info={
	ms1:{clip_size:1},
	ms4:{clip_size:4},
	ms10:{clip_size:10},
	ms20:{clip_size:20},
	ms40:{clip_size:40},
	ms80:{clip_size:80},
	ms160:{clip_size:160},
	ms320:{clip_size:320},
	ms640:{clip_size:640},
};

var unit_numbers={
	h1:[1],
	m1:[2,3],
	m2:[2,3],
	m3:[2,3],
	m4:[2,3],
	m5:[2,3],
	nn:[1],
	synth1:[1,2,3,4,5]
};

var algs=[];
{
	var txt=read_text_file(alglist_path);
	var lines=txt.split('\n');
	for (var i in lines) {
		if (lines[i].trim().slice(0,1)!='#') {
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
					throw 'problem in alglist file: '+lines[i].trim();
				}
			}
		}
	}
}

var datasets=[];
{
	var txt=read_text_file(datasetlist_path);
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

if (command=="ground_truth_test") {
	if ('norun' in CLP.namedParameters) {
		compile_results(function() {
			process.exit(0);
		});
	}
	else {
		run_sorting(function() {
			compile_results(function() {
				process.exit(0);
			});
		});
	}
}
else if (command=="sort") {
	run_sorting(function() {
		process.exit(0);
	});
}
else {
	console.log('Unknown command: '+command);
	process.exit(-1);
}

function print_usage() {
	console.log('boxsort.sh [sort|ground_truth_test] [alglist].txt [datasetlist].txt --outpath=output [--norun]');
}

function compile_results(callback) {
	for (var d in datasets) {
		var dsname=datasets[d].name;
		var ks=unit_numbers[dsname]||[];
		console.log ('');
		console.log ('######## DATASET: '+dsname);
		console.log ('');
		for (var ii=0; ii<=ks.length; ii++) {
			var k=ks[ii];
			if (ii==ks.length) k='all';
			var table0=[];
			for (var a in algs) {
				var algname=algs[a].name;
				var outpath0=outpath+'/'+algname+'-'+dsname;
				var CM=read_csv_matrix(outpath0+'/confusion_matrix.csv');
				var LM=read_csv_vector(outpath0+'/optimal_label_map.csv');

				var Ntrue=0,Ncorrect=0,Nincorrect=0,Nmissed=0,Nextra=0;

								
				if (k=='all') {
					var K1=CM.length-1;
					var K2=CM[0].length-1;
					var ks_to_consider=clone(ks);
					if (ks_to_consider.length===0) {
						for (var aa=1; aa<=K1; aa++) {
							ks_to_consider.push(aa);
						}
					}

					var ret=sub_confusion_matrix(CM,LM,ks_to_consider);
					var CMb=ret.CMb;
					var LMb=ret.LMb;
					var K1b=CMb.length-1;
					var K2b=CMb[0].length-1;
					for (var k1b=1; k1b<=K1b; k1b++) {
						Ntrue+=row_sum(CMb,k1b-1);
						var k2b_match=LMb[k1b-1];
						if (k2b_match) {
							Ncorrect+=CMb[k1b-1][k2b_match-1];
						}
						for (var k2b=1; k2b<=K1b; k2b++) {
							if (k2b!=k2b_match) {
								Nincorrect+=CMb[k1b-1][k2b-1];
							}
						}
						Nmissed+=CMb[k1b-1][K2b];
					}
					Nextra+=row_sum(CMb,K1b);
				}
				else {
					var K1=CM.length-1;
					var K2=CM[0].length-1;
					var k1=k;
					Ntrue+=row_sum(CM,k1-1);
					var k2_match=LM[k1-1];
					if (k2_match) {
						Ncorrect+=CM[k1-1][k2_match-1];
					}
					for (var k2=1; k2<=K2; k2++) {
						if (k2!=k2_match) {
							Nincorrect+=CM[k1-1][k2-1];
						}
					}
					Nmissed+=CM[k1-1][K2];
					if (k2_match) {
						Nextra+=col_sum(CM,k2_match-1)-CM[k1-1][k2_match-1];
					}
				}
				var k_display=k;
				if ((k=='all')&&(ks.length>0)) k_display=ks.join(',');
				var clip_size=(alg_info[algname]||{}).clip_size||'';
				var table_row={DATASET:dsname,ALG:algname,clip_size:clip_size,UNIT:k_display,Ntrue:Ntrue,Correct:topct(Ncorrect/Ntrue),Incorrect:topct(Nincorrect/Ntrue),Missed:topct(Nmissed/Ntrue)};
				table_row.Extra=topct(Nextra/Ntrue);
				table0.push(table_row);

				/*
				{
					num1=row_sum(CM,k-1);
					num2=col_sum(CM,LM[k-1]-1);
					if (LM[k-1]>0) {
						num_fn=num1-CM[k-1][LM[k-1]-1];
						num_fp=col_sum(CM,LM[k-1]-1)-CM[k-1][LM[k-1]-1];
						
					}
					else {
						num_fn=NaN;
						num_fp=NaN;
					}	
				}
				var num_fn_frac=num_fn/num1;
				var num_fp_frac=num_fp/num2;

				//console.log('@@@@@@@@@@@@ '+algname+' '+dsname+' '+k+' '+num_fn_frac+' '+num_fp_frac);
				//print_csv_matrix(CM);
				//console.log(LM);
				//console.log(algname+' '+dsname+' '+k+': num='+num+' fn='+num_fn+' fp='+num_fp);
				//console.log('  '+algname+' '+dsname+' '+k+': '+num+'\t'+topct(num_fn/num1)+'\t'+topct(num_fp/num2));
				var clip_size=(alg_info[algname]||{}).clip_size||'';
				table0.push({ALG:algname,clip_size:clip_size,DATASET:dsname,UNIT:k,num_events:num1,false_neg:topct(num_fn_frac),false_pos:topct(num_fp_frac)});
				*/
			}
			console.table(table0);
		}
	}
	callback();
}

function run_sorting(callback) {
	var num_running=0;
	for (var a in algs) {
		for (var d in datasets) {
			console.log ('Applying '+algs[a].name+' to '+datasets[d].name);
			num_running++;
			apply_sorting(algs[a],datasets[d],function() {
				num_running--;
			});
		}
	}

	function on_timeout() {
		console.log ('# algorithms running: '+num_running);
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

	copy_file_sync(basepath+'/datasets/'+ds.folder+'/firings_true.mda',outpath0+'/firings_true.mda');
	copy_file_sync(basepath+'/datasets/'+ds.folder+'/firings_true.mda.prv',outpath0+'/firings_true.mda.prv');

	var ds_folder=basepath+'/datasets/'+ds.folder;

	var cmd='mountainprocess';
	args=['queue-script',basepath+'/algorithms/'+alg.script];
	if (fs.existsSync(ds_folder+"/raw.mda.prv"))
		args.push('--raw='+ds_folder+'/raw.mda.prv');
	else
		args.push('--raw='+ds_folder+'/raw.mda');
	args.push(ds_folder+'/params.json');
	var geom_fname=ds_folder+'/geom.csv';
	if (fs.existsSync(geom_fname)) {
		args.push('--geom='+geom_fname);
	}
	args.push('--outpath='+outpath0);
	var params0=alg.params.split(' ');
	for (var i in params0) {
		args.push(params0[i]);
	}
	console.log (cmd+' '+args.join(' '));
	
	make_system_call(cmd,args,function() {
		console.log('++++++++++++++++++++++++++++++++');
		compute_confusion_matrix(outpath0,function() {
			console.log('++++++++++++++++++++++++++++++++------');
			callback();
		});
	});
}

function compute_confusion_matrix(output_path,callback) {
	var cmd='mountainprocess';
	var args=['run-process','merge_firings'];
	if (fs.existsSync(output_path+"/firings_true.mda.prv"))
		args.push('--firings1='+output_path+'/firings_true.mda.prv');
	else if (fs.existsSync(output_path+"/firings_true.mda"))
		args.push('--firings1='+output_path+'/firings_true.mda');
	else {
		callback(); //no ground truth
		return;
	}
	args.push('--firings2='+output_path+'/firings.mda');
	args.push('--confusion_matrix='+output_path+'/confusion_matrix.csv');
	args.push('--optimal_label_map='+output_path+'/optimal_label_map.csv');
	args.push('--firings_merged='+output_path+'/firings_merged.mda');
	args.push('--max_matching_offset=10');

	make_system_call(cmd,args,function() {
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