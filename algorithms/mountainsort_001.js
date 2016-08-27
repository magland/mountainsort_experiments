function main(params) {
	run_algorithm(params);
}

function run_algorithm(params) {
	params.clip_size=params.clip_size||50;
	params.detect_threshold=params.detect_threshold||3.5;
	params.detect_interval=params.detect_interval||10;
	params.shell_increment=params.shell_increment||3;
	params.min_shell_size=params.min_shell_size||150;
	params.samplerate=params.samplerate||30000;
	params.sign=params.sign||0;
	params.freq_min=params.freq_min||300;
	params.freq_max=params.freq_max||8000;
	params.num_features=params.num_features||10;
	params.num_features2=params.num_features2||10;
	params.adj_radius=params.adj_radius||0;
	params.consolidation_factor=params.consolidation_factor||0.9;
	params.channels=params.channels||'';
	params.timerange=params.timerange||[-1,-1];
	params.use_whitening=params.use_whitening||'true';
	if (params.mask_threshold===null) params.mask_threshold=6; //to mask out artifacts. Use 0 to not mask out artifacts
	params.cluster_num_threads=params.cluster_num_threads||0;

	display_parameters(params);

	var raw=params.raw;
	var geom=params.geom;
	var outpath=params.outpath;

	var o_geom2adj={
		channels:params.channels,
		radius:params.adj_radius
	};
	var o_extract_raw={
		t1:params.timerange[0],
		t2:params.timerange[1],
		channels:params.channels	
	};
	var o_filter={
		samplerate:params.samplerate,
		freq_min:params.freq_min,
		freq_max:params.freq_max
	};
	var o_mask_out_artifacts={
		threshold:params.mask_threshold,
		interval_size:200
	};
	var o_whiten={
	};
	var o_detect={
		detect_threshold:params.detect_threshold,
		detect_interval:params.detect_interval,
		clip_size:params.clip_size,
		sign:params.sign,
		individual_channels:1
	};
	var o_branch_cluster={
		clip_size:params.clip_size,
		min_shell_size:params.min_shell_size,
		shell_increment:params.shell_increment,
		num_features:params.num_features,
		num_features2:params.num_features2,
		detect_interval:params.detect_interval,
		consolidation_factor:params.consolidation_factor,
		num_threads:params.cluster_num_threads
	};
	var o_merge_across_channels_v2={
		clip_size:params.clip_size
	};
	var o_fit_stage={
		clip_size:params.clip_size,
		min_shell_size:params.min_shell_size,
		shell_increment:params.shell_increment
	};

	initialize_pipeline();

	var adjacency_matrix='';
	if (geom) {
		adjacency_matrix='@AM';
		geom2adj(geom,adjacency_matrix,o_geom2adj);
		geom2adj(geom,outpath+'/AM.csv',o_geom2adj); //for convenience
	}

	extract_raw(raw,'@pre0',o_extract_raw);
	if (params.mask_threshold) {
		bandpass_filter('@pre0','@pre1',o_filter);
		mask_out_artifacts('@pre1','@pre1b',o_mask_out_artifacts);
	}
	else {
		bandpass_filter('@pre0','@pre1b',o_filter);
	}
	if (params.use_whitening=='true')
		whiten('@pre1b','@pre2',o_whiten);
	else
		normalize_channels('@pre1b','@pre2',{});

	detect('@pre2','@detect',o_detect);
	branch_cluster_v2('@pre2','@detect',adjacency_matrix,'@firings1',o_branch_cluster);
	merge_across_channels_v2('@pre2','@firings1','@firings2',o_merge_across_channels_v2);
	fit_stage('@pre2','@firings2','@firings3',o_fit_stage);

	save_prv('@pre0',outpath+'/pre0.mda.prv');
	save_prv('@pre1b',outpath+'/pre1b.mda.prv');
	save_prv('@pre2',outpath+'/pre2.mda.prv');
	save_prv(raw,outpath+'/raw.mda.prv');

	//copy('@firings1',outpath+'/firings_before_merge_across_channels.mda');
	//copy('@firings2',outpath+'/firings_before_fit.mda');
	copy('@firings3',outpath+'/firings.mda');

	run_pipeline();
}

function geom2adj(input,output,opts) {
	var ret={};
	ret.processor_name='geom2adj';
	ret.inputs={input:input};
	ret.outputs={output:output};
	ret.parameters=clone(opts);
	PIPELINE.addProcess(ret);
}

function extract_raw(timeseries,timeseries_out,opts) {
	var ret={};
	ret.processor_name='extract_raw';
	ret.inputs={timeseries:timeseries};
	ret.outputs={timeseries_out:timeseries_out};
	ret.parameters=clone(opts);
	PIPELINE.addProcess(ret);
}

function bandpass_filter(timeseries,timeseries_out,opts) {
	var ret={};
	ret.processor_name='bandpass_filter';
	ret.inputs={timeseries:timeseries};
	ret.outputs={timeseries_out:timeseries_out};
	ret.parameters=clone(opts);
	PIPELINE.addProcess(ret);
}

function whiten(timeseries,timeseries_out,opts) {
	var ret={};
	ret.processor_name='whiten';
	ret.inputs={timeseries:timeseries};
	ret.outputs={timeseries_out:timeseries_out};
	ret.parameters=clone(opts);
	PIPELINE.addProcess(ret);
}

function normalize_channels(timeseries,timeseries_out,opts) {
	var ret={};
	ret.processor_name='normalize_channels';
	ret.inputs={timeseries:timeseries};
	ret.outputs={timeseries_out:timeseries_out};
	ret.parameters=clone(opts);
	PIPELINE.addProcess(ret);
}

function mask_out_artifacts(timeseries,timeseries_out,opts) {
	var ret={};
	ret.processor_name='mask_out_artifacts';
	ret.inputs={timeseries:timeseries};
	ret.outputs={timeseries_out:timeseries_out};
	ret.parameters=clone(opts);
	PIPELINE.addProcess(ret);
}

function detect(timeseries,detect,opts) {
	var ret={};
	ret.processor_name='detect';
	ret.inputs={timeseries:timeseries};
	ret.outputs={detect_out:detect};
	ret.parameters=clone(opts);
	PIPELINE.addProcess(ret);
}

function branch_cluster_v2(timeseries,detect,adjacency_matrix,firings,opts) {
	var ret={};
	ret.processor_name='branch_cluster_v2';
	ret.inputs={timeseries:timeseries,detect:detect,adjacency_matrix:adjacency_matrix};
	ret.outputs={firings_out:firings};
	ret.parameters=clone(opts);
	PIPELINE.addProcess(ret);
}

function merge_across_channels(timeseries,firings,firings_out,opts) {
	var ret={};
	ret.processor_name='merge_across_channels';
	ret.inputs={timeseries:timeseries,firings:firings};
	ret.outputs={firings_out:firings_out};
	ret.parameters=clone(opts);
	PIPELINE.addProcess(ret);
}

function merge_across_channels_v2(timeseries,firings,firings_out,opts) {
	var ret={};
	ret.processor_name='merge_across_channels_v2';
	ret.inputs={timeseries:timeseries,firings:firings};
	ret.outputs={firings_out:firings_out};
	ret.parameters=clone(opts);
	PIPELINE.addProcess(ret);
}

function fit_stage(timeseries,firings,firings_out,opts) {
	var ret={};
	ret.processor_name='fit_stage';
	ret.inputs={timeseries:timeseries,firings:firings};
	ret.outputs={firings_out:firings_out};
	ret.parameters=clone(opts);
	PIPELINE.addProcess(ret);
}

function merge_stage(timeseries,firings,firings_out,opts) {
	var ret={};
	ret.processor_name='merge_stage';
	ret.inputs={timeseries:timeseries,firings:firings};
	ret.outputs={firings_out:firings_out};
	ret.parameters=clone(opts);
	PIPELINE.addProcess(ret);
}

function compute_outlier_scores(timeseries,firings,firings_out,opts) {
	var ret={};
	ret.processor_name='compute_outlier_scores';
	ret.inputs={timeseries:timeseries,firings:firings};
	ret.outputs={firings_out:firings_out};
	ret.parameters=clone(opts);
	PIPELINE.addProcess(ret);
}

function compute_detectability_scores(timeseries,firings,firings_out,opts) {
	var ret={};
	ret.processor_name='compute_detectability_scores';
	ret.inputs={timeseries:timeseries,firings:firings};
	ret.outputs={firings_out:firings_out};
	ret.parameters=clone(opts);
	PIPELINE.addProcess(ret);
}

function copy(input,output) {
	var ret={};
	ret.processor_name='copy';
	ret.inputs={input:input};
	ret.outputs={output:output};
	ret.parameters={};
	PIPELINE.addProcess(ret);
}

function save_prv(input,output) {
	if (!PIPELINE.savePrv(input,output)) {
		throw "Error in save_prv: "+input+' '+output;
	}
}

function clone(obj) {
	return JSON.parse(JSON.stringify(obj));
}

var console={
	log:function(msg) {MP.log(msg);}
};

function MSPipeline() {
	this.addProcess=function(P) {m_processes.push(clone(P));};
	this.addProcesses=function(Ps) {for (var i in Ps) {this.addProcess(Ps[i]);}};
	this.clearProcesses=function() {m_processes=[];};
	this.savePrv=function(input_path,output_path) {return _savePrv(input_path,output_path);};
	this.run=_run;

	var m_processes=[];
	var m_prvs=[];

	function sort_processes() {
		var new_processes=[];
		for (var j=0; j<m_processes.length; j++) m_processes[j]._handled=false;
		var something_changed=true;
		while (something_changed) {
			something_changed=false;
			for (var j=0; j<m_processes.length; j++) {
				var PP=m_processes[j];
				if (!PP._handled) {
					var ready_to_go=true;
					for (var pname in PP.inputs) {
						var tmp0=find_process_and_pname_for_output(m_processes,PP.inputs[pname]);
						if (tmp0.process) {
							if (!tmp0.process._handled) ready_to_go=false;
						}
					}
					if (ready_to_go) {
						something_changed=true;
						PP._handled=true;
						new_processes.push(PP);
					}
				}
			}
		}
		for (var j=0; j<m_processes.length; j++) delete(m_processes[j]._handled);
		if (new_processes.length!=m_processes.length) return false;
		m_processes=new_processes;
		return true;
	}

	function compute_process_codes() {
		//assume sorted
		for (var j=0; j<m_processes.length; j++) m_processes[j].process_code='';
		for (var j=0; j<m_processes.length; j++) {
			var PP=m_processes[j];
			var input_codes={};
			for (var pname in PP.inputs) {
				var fname=PP.inputs[pname];
				var tmp=find_process_and_pname_for_output(m_processes,fname);
				if (tmp.process) {
					if (!tmp.process.process_code) {
						console.log ('Unexpected problem in compute_process_code');
						return false;
					}
					input_codes[pname]=tmp.process.process_code+'--'+tmp.pname;
				}
				else {
					if (fname) {
						input_codes[pname]=MP.fileChecksum(fname);
						if (!input_codes[pname]) {
							console.log ('Checksum is empty, perhaps file does not exist for '+pname+': '+fname);
							return false;
						}
					}
					else {
						input_codes[pname]='<empty>';
					}
				}
			}
			var code_obj={};
			for (var pname in PP.inputs) {
				code_obj[pname]=input_codes[pname];
			}
			for (var pname in PP.parameters) {
				code_obj[pname]=clone(PP.parameters[pname]);
			}
			PP.process_code=MP.stringChecksum(JSON.stringify(code_obj));
		}
		return true;
	}

	var fname_map={};
	function _run() {

		//sort so that each process only depends on output from earlier processes
		if (!sort_processes()) {
			console.log ('Problem sorting processes. There may be cyclic dependencies.');
			return false;
		}

		//compute the process codes, which should depend on the parameters and the provenance (or checksums) of the input files
		if (!compute_process_codes()) {
			console.log ('Problem computing process codes.');
			return false;
		}

		//create the fname mapping
		fname_map={}; //maps @name to /the/corresponding/temporary/path
		for (var j=0; j<m_processes.length; j++) {
			var PP=m_processes[j];
			if (!PP.process_code) {
				console.log ('Unexpected problem. Process code is empty.');
				return false;
			}
			for (var pname in PP.outputs) {
				if (PP.outputs[pname].indexOf('@')==0) {
					fname_map[PP.outputs[pname]]=MP.createTemporaryFileName(PP.process_code+'-'+pname);
				}
			}
		}

		//apply the fname mapping
		var processes_mapped=[];
		for (var j=0; j<m_processes.length; j++) {
			var PP=clone(m_processes[j]);
			for (var pname in PP.inputs) {
				var fname=PP.inputs[pname];
				if (fname.indexOf('@')==0) {
					if (fname in fname_map) {
						PP.inputs[pname]=fname_map[fname];
					}
				}
			}
			for (var pname in PP.outputs) {
				var fname=PP.outputs[pname];
				if (fname.indexOf('@')==0) {
					if (fname in fname_map) {
						PP.outputs[pname]=fname_map[fname];
					}
				}
			}
			processes_mapped.push(PP);
		}

		var submit_the_entire_pipeline=true;
		if (submit_the_entire_pipeline) {
			var pipe={processes:[]};
			for (var j=0; j<processes_mapped.length; j++) {
				var PP=processes_mapped[j];
				var processor_name=PP.processor_name;
				var parameters=clone(PP.parameters);
				for (var pname in PP.inputs) {
					parameters[pname]=PP.inputs[pname];
				}
				for (var pname in PP.outputs) {
					parameters[pname]=PP.outputs[pname];
				}
				pipe.processes.push({processor_name:processor_name,parameters:parameters,prvs:PP.prvs||{}});
			}
			if (!MP.runPipeline(JSON.stringify(pipe))) {
				console.log ('Error running pipeline');
				return false;
			}
		}
		else {
			//run the processes one at a time (recall that they were sorted earlier)
			for (var j=0; j<processes_mapped.length; j++) {
				var PP=processes_mapped[j];
				var processor_name=PP.processor_name;
				var parameters=clone(PP.parameters);
				for (var pname in PP.inputs) {
					parameters[pname]=PP.inputs[pname];
				}
				for (var pname in PP.outputs) {
					parameters[pname]=PP.outputs[pname];
				}
				if (!MP.runProcess(processor_name,JSON.stringify(parameters))) {
					console.log ('Error running process: '+processor_name);
					return false;
				}
			}
		}

		return true;
	}	

	function find_process_and_pname_for_output(processes,fname) {
		for (var j=0; j<processes.length; j++) {
			var PP=processes[j];
			for (var pname in PP.outputs) {
				if (PP.outputs[pname]==fname) {
					return {process:processes[j],pname:pname};
				}
			}
		}
		return {process:null,pname:''};
	}

	function _savePrv(input_path,output_path) {
		var A=find_process_and_pname_for_output(m_processes,input_path);
		if (!A.process) {
			if (input_path.indexOf('@')===0) {
				console.log ('Error in _savePrv. Unable to find process for output: '+input_path);
				return false;
			}
			m_processes.push({processor_name:'',parameters:{output:input_path},prvs:{output:output_path}}); //dummy processor that will have a .prv
			return true;
		}
		if (!A.process.prvs) A.process.prvs={};
		A.process.prvs[A.pname]=output_path;
		return true;
	}
}

function display_parameters(params) {
	for (var key in params) {
		console.log (key+' = '+params[key]);
	}
}

function initialize_pipeline() {
	PIPELINE.clearProcesses();
}

function run_pipeline() {
	PIPELINE.run();
}

var PIPELINE=new MSPipeline();
