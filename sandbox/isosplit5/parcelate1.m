function [labels,parcels]=parcelate1(X,target_size)
if nargin<1, test_parcelate1; return; end;

[M,N]=size(X);
labels=zeros(1,N);

parcels={};

P.indices=1:N;
labels(P.indices)=1;
parcels{end+1}=P;

split_factor=3; % split factor around 2.71 is in a sense ideal

p_index=1;
while (p_index<=length(parcels))
    inds=parcels{p_index}.indices;
    sz=length(inds);
    if (sz>target_size*split_factor)
        iii=randsample(sz,split_factor);
        pts=X(:,inds(iii));
        dm=distance_matrix(pts,X(:,inds));
        [~,assignments]=min(dm,[],1);
        parcels{p_index}.indices=inds(find(assignments==1));
        labels(parcels{p_index}.indices)=p_index;
        for jj=2:split_factor
            PP.indices=inds(find(assignments==jj));
            parcels{end+1}=PP;
            labels(PP.indices)=length(parcels);
        end;
    else
        p_index=p_index+1;
    end;
end;

function ret=distance_matrix(X1,X2)
[M,N1]=size(X1);
[M,N2]=size(X2);
ret=zeros(N1,N2);
[aa,bb]=ndgrid(1:N1,1:N2);
for m=1:M
    ret=ret+reshape((X1(m,aa(:))-X2(m,bb(:))).^2,N1,N2);
end;
ret=sqrt(ret);

function test_parcelate1

target_size=100;
M=2;
N=1e5;
X=randn(M,N);
tic;
labels=parcelate1(X,target_size);
toc
figure; ms_view_clusters(X,labels);

