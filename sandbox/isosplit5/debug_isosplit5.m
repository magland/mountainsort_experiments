function debug_isosplit5

%rng(12);

M=10;
N0=1e5;

X=randn(M,N0);
X(1,:)=X(1,:)*1;
X2=randn(M,N0);
X2(1,:)=X2(1,:)+2.5;
X=cat(2,X,X2);

tic;
labels=isosplit5_mex(X);
toc
max(labels)

if (max(labels)>1)
    figure; ms_view_clusters(X,labels);
    inds1=find(labels==1);
    inds2=find(labels==2);
    inds12=[inds1,inds2];
    X12=X(:,inds12);
    V=mean(X(:,inds2),2)-mean(X(:,inds1),2);
    V=V/sqrt(V'*V);
    proj=V'*X(:,inds12);
    figure; hist(proj,100);
    [ds,cp]=isocut5(proj)
    error('debug');
end;