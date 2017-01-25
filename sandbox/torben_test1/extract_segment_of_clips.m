function clips=extract_segment_of_clips(clips,times,t1,t2)
inds=find((t1<=times)&(times<=t2));
clips=clips(:,:,inds);