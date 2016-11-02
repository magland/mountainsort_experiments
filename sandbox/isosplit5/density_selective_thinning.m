function X_thin=density_selective_thinning(X,opts)

if nargin<1, test_density_selective_thinning; return; end;
if nargin<2, opts=struct; end;
if ~isfield(opts,'cc'), opts.cc=0.5; end;
if ~isfield(opts,'ds'), opts.ds=2; end;

[M,N]=size(X);

estdensities=kmdensity(X);
inds0=find(estdensities==0);
inds1=find(estdensities>0);
estlogdensities=zeros(1,N);
estlogdensities(inds1)=log2(estdensities(inds1));

targetlogdensities(inds1)=(estlogdensities(inds1)-mean(estlogdensities(inds1)))*opts.cc;
probs=2.^(targetlogdensities-estlogdensities);
probs=probs/(sum(probs))*N/opts.ds;
probs(inds0)=1; %always use those with estimated density of zero
inds_to_use=find(rand(size(probs))<probs);
X_thin=X(:,inds_to_use);

function ret=knn(X,K_nn)
knn=knnsearch(X',X','K',K_nn+1)';
knn=knn(K_nn+1,:);
ret=knn;

function test_density_selective_thinning

close all;
M=10;
N0=1e4;

X1=randn(M,ceil(N0/4));
X2=randn(M,ceil(N0/4));
X2(1,:)=X2(1,:)+4.5;
X3=randn(M,N0*12);
X3(1,:)=X3(1,:)+10.5;
X4=randn(M,N0*24);
X4(1,:)=X4(1,:)+15;
X5=randn(M,500);
X5(2,:)=X5(2,:)+6;
X6=randn(M,300);
X6(1,:)=X6(1,:)*1.5+15;
X6(2,:)=X6(2,:)*0.7+6;
X=cat(2,X1,X2,X3,X4,X5,X6);
N=size(X,2);
fprintf('N = %g\n',N);

oo.cc=0.3;
oo.ds=50;
ttt=tic;
X_thin=density_selective_thinning(X,oo);
fprintf('Total time for dst: %g\n',toc(ttt));

ttt=tic;
labels_thin=isosplit5(X_thin);
fprintf('Time for isosplit5 of thinned data: %g\n',toc(ttt));

figure;
subplot(4,1,1);
hist(X(1,:),100);
subplot(4,1,2);
%ms_view_clusters(X(1:2,:),labels);
ms_view_clusters(X(1:2,:));
title(sprintf('Original data: N = %d\n',size(X,2)));
subplot(4,1,3);
hist(X_thin(1,:),100);
subplot(4,1,4);
ms_view_clusters(X_thin(1:2,:),labels_thin);
title(sprintf('Thinned data: N = %d\n',size(X_thin,2)));