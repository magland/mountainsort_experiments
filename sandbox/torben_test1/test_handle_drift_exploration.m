close all;

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
segment_durations=30000*60*[5,10,20,40,80,160,320,640];

fprintf('Defining segments...\n');
runs={};
LL=length(times0);

for rr=1:length(segment_durations)
    dur=segment_durations(rr);
    
    index_chunks={};
    tt=[1:floor(dur/2):N,inf];
    for j=1:length(tt)-1
        t1=tt(j);
        t2=tt(j+1);
        index_chunks{end+1}=find((t1<=times0)&(times0<t2));
    end;
    
    segments={};
    for j=1:length(index_chunks)-1
        S.inds1=index_chunks{j};
        S.inds2=index_chunks{j+1};
        S.inds=[S.inds1,S.inds2];
        fprintf('Segment %d of %d in run %d of %d (%d events)\n',j,length(index_chunks)-1,rr,length(segment_durations),length(S.inds));
        S.opts=struct;
        S.opts.refine_clusters=false;
        if (j>1)
            S.opts.initial_labels=[segments{j-1}.labels2,zeros(size(S.inds2))];
            S.opts.K_init=4;
        end;
        S.labels=isosplit5(features(:,S.inds));
        S.labels1=S.labels(1:length(S.inds1));
        S.labels2=S.labels(length(S.inds1)+1:end);
        segments{j}=S;
        fprintf('Found %d clusters\n',max(S.labels));
    end;
    
    R=struct;
    R.duration=dur;
    R.segments=segments;
    
    thresh=0.7;
    R.labels=zeros(1,LL);
    if (length(R.segments)>0)
        R.segments{1}.labels_matched=R.segments{1}.labels;
        R.segments{1}.labels1_matched=R.segments{1}.labels1;
        R.segments{1}.labels2_matched=R.segments{1}.labels2;
        R.labels(R.segments{1}.inds)=R.segments{1}.labels_matched;
        kk_max=max(R.segments{1}.labels_matched);
        for ss=2:length(R.segments)
            S=R.segments{ss};
            Sprev=R.segments{ss-1};
            lab0=S.labels1;
            lab0prev=Sprev.labels2_matched;
            K0=max(S.labels);
            labels_map=zeros(1,K0);
            for k=1:K0
                inds_k=find(S.labels1==k);
                k_match=0;
                if (length(inds_k)>0)
                    best_k_prev=mode(Sprev.labels2_matched(inds_k));
                    numer0=length(find(Sprev.labels2_matched(inds_k)==best_k_prev));
                    denom0=length(inds_k)+length(find(Sprev.labels2_matched==best_k_prev))-numer0;
                    if ((denom0)&&(numer0/denom0>thresh))
                        k_match=best_k_prev;
                    end;
                end;
                if (k_match==0)
                    k_match=kk_max+1;
                    kk_max=k_match;
                end;
                labels_map(k)=k_match;
            end;
            S.labels_matched=labels_map(S.labels);
            S.labels1_matched=labels_map(S.labels1);
            S.labels2_matched=labels_map(S.labels2);
            R.labels(S.inds2)=S.labels2_matched;
            R.segments{ss}=S;
        end;
    end;
    
    runs{end+1}=R;
end;
    
R=runs{end-2};

firings=zeros(3,length(R.labels));
firings(2,:)=times0;
firings(3,:)=R.labels;

% Now write out this to a firings.mda file and open in mountainview
