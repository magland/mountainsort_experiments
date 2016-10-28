function [inds,weights]=thin5(X,desired_num)

[M,N]=size(X);
inds=1:N;
while (length(inds)>desired_num)
    N0=length(inds);
    disp(N0);
    X0=X(:,inds);
    nearest0=knnsearch(X0',X0','K',2);
    nearest0=nearest0(:,2);
    dists=compute_dists(X0,X0(:,nearest0));
    dists_sorted=sort(dists);
    cutoff=dists_sorted(ceil(N0*1/10));
    
    to_keep=ones(1,N0);
    inds_to_consider_discarding=find(dists<=cutoff);
    [inds1,inds2]=split_inds(inds_to_consider_discarding);
    to_keep(inds1)=0;
    inds=inds(find(to_keep));
end

X0=X(:,inds);
nearest0=knnsearch(X0',X','K',1);
weights=accumarray(nearest0,1,[length(inds),1])';

function [inds1,inds2]=split_inds(inds)
N=length(inds);
labels=2*ones(1,N);
labels(randsample(N,ceil(N/2)))=1;
inds1=inds(find(labels==1));
inds2=inds(find(labels==2));

function ret=compute_dists(X,Y)
ret=sqrt(sum((X-Y).^2,1));