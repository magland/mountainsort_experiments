%L=load('test_handle_drift.mat');
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
    S.clips=extract_segment_of_clips(clips,times0,S.t1,S.t2);
    segments{j}=S;
end;

segments{1}.clips_traced=segments{1}.clips;
for j=2:length(segments)
    fprintf('Segment %d/%d\n',j,length(segments));
    clips1=segments{j-1}.clips;
    clips2=segments{j}.clips;
    N1=size(clips1,3);
    N2=size(clips2,3);
    FF=ms_event_features(cat(3,clips1,clips2),num_features);
    FF1=FF(:,1:N1);
    FF2=FF(:,N1+1:N1+N2);

    idx1=knnsearch(FF1',FF2','K',num_neighbors);
    clips2_1=zeros(size(clips2));
    for ii=1:num_neighbors
        clips2_1=clips2_1+clips1(:,:,idx1(:,ii));
    end;
    clips2_1=clips2_1/num_neighbors;
    
    idx2=knnsearch(FF2',FF2','K',num_neighbors+1);
    clips2_2=zeros(size(clips2));
    for ii=2:num_neighbors+1
        clips2_2=clips2_2+clips2(:,:,idx2(:,ii));
    end;
    clips2_2=clips2_2/num_neighbors;

    segments{j}.clips_traced=clips2+clips2_1-clips2_2;
end;

clips_traced=zeros(size(clips,1),size(clips,2),0);
for j=1:length(segments)
    clips_traced=cat(3,clips_traced,segments{j}.clips_traced);
end;

FF=ms_event_features(clips,num_features);
FF_traced=ms_event_features(clips_traced,num_features);

