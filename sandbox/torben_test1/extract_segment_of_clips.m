function [clips,features]=extract_segment_of_clips(clips,features,times,t1,t2)
inds=find((t1<=times)&(times<=t2));
clips=clips(:,:,inds);
features=features(:,inds);