function test_handle_drift

%create_mat_file;
L=load('test_handle_drift.mat');
clips=L.data.clips;
labels=L.data.labels;
times=L.data.times;


function create_mat_file

fprintf('Reading preprocessed data...\n');
pre_fname='/disk1/tmp/mountainlab/tmp_long_term/7ab38200a26231f058b1bc6020a4dc6a31a10a29-whiten-timeseries_out.tmp';
X=readmda(pre_fname);
fprintf('Reading firings file...\n');
firings=readmda('output/ms2--drift1/firings.mda');
times=firings(2,:);
labels=firings(3,:);
K=max(labels);

fprintf('Extracting clips...\n');
data.clips=ms_extract_clips2(X,times,30);
data.times=times;
data.labels=labels;
fprintf('Saving results...\n');
save('test_handle_drift.mat','data');