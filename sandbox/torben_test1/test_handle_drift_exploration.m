%fprintf('Loading...\n');
%L=load('test_handle_drift.mat');
clips=L.data.clips;
labels=L.data.labels;
times0=L.data.times;

num_features=20;
segment_size=30000*60*10;
segment_step=30000*60;

fprintf('Computing features...\n');
features=ms_event_features(clips,num_features);

M=size(clips,1);
T=size(clips,2);
N=max(times0);
t1s=1:segment_step:N;
if (length(t1s)>1)
    t1s=t1s(1:end-1);
end;

fprintf('Defining segments...\n');
segments={};
LL=length(times0);
labels=zeros(1,LL);
for j=1:length(t1s)
    S.t1=t1s(j);
    if (j<length(t1s))
        S.t2=min(N,t1s(j)+segment_size-1);
    else
        S.t2=N;
    end;
    S.inds=find((S.t1<=times0)&(times0<=S.t2));
    segments{j}=S;
end;

S1=segments{1};
features1=features(:,S1.inds);
%labels1=isosplit5(features1);
labels(S1.inds)=labels1;
figure; ms_view_clusters(features1,labels1);
figure; ms_view_templates_from_clips(clips(:,:,S1.inds),labels1)

%for kk=2:length(segments)
for kk=2:3
    S_kk=segments{kk};
    features_kk=features(:,S_kk.inds);
    labels_kk=labels(S_kk.inds);
    opts.initial_labels=labels_kk;
    new_labels_kk=isosplit5(features_kk,opts);
    inds0=find(labels(S_kk.inds)==0);
    labels(S_kk.inds(inds0))=labels_kk(inds0);
    
    figure; ms_view_clusters(features_kk,labels_kk);
    figure; ms_view_templates_from_clips(clips(:,:,S_kk.inds),labels_kk)
end;

