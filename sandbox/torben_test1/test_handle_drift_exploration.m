L=load('test_handle_drift.mat');
clips=L.data.clips;
labels=L.data.labels;
times0=L.data.times;

num_features=20;
num_neighbors=10;
segment_size=30000*60;

N=max(times0);
t1s=1:segment_size:N;
if (length(t1s)>1)
    t1s=t1s(1:end-1);
end;

segments={};
for j=1:length(t1s)
    S.t1=t1s(j);
    if (j<length(t1s))
        S.t2=t1s(j+1)-1;
    else
        S.t2=N;
    end;
    S.inds=find((S.t1<=times0)&(times0<=S.t2));
    S.clips=clips(:,:,S.inds);
    segments{j}=S;
end;

S1=segments{1};
segments{1}.labels=labels(S1.inds);
for j=2:length(segments)
    fprintf('Segment %d/%d\n',j,length(segments));
    clips1=segments{j-1}.clips;
    clips2=segments{j}.clips;
    N1=size(clips1,3);
    N2=size(clips2,3);
    FF=ms_event_features(cat(3,clips1,clips2),num_features);
    FF1=FF(:,1:N1);
    FF2=FF(:,N1+1:N1+N2);

    idx1=knnsearch(FF1',FF2','K',1);
    segments{j}.labels=segments{j-1}.labels(idx1);
end;

new_labels=zeros(size(labels));
for j=1:length(segments)
    new_labels(segments{j}.inds)=segments{j}.labels;
end;

FF=ms_event_features(clips,num_features);
figure; ms_view_clusters(FF,labels);
figure; ms_view_clusters(FF,new_labels);

