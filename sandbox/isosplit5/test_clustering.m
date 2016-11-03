function test_clustering

rng(6);

close all;
M=10;
N0=5e4;

% Create the data
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

ttt=tic;
labels=isosplit5(X,struct('verbose',0,'K_init',30,'refine_clusters',0));
fprintf('Time for isosplit5: %g\n',toc(ttt));

ttt=tic;
labels_mex=isosplit5_mex(X);
fprintf('Time for isosplit5_mex: %g\n',toc(ttt));

figure;
subplot(3,1,1);
hist(X(1,:),100);
subplot(3,1,2);
ms_view_clusters(X(1:2,:),labels);
title(sprintf('isosplit5: N = %d',size(X,2)));
subplot(3,1,3);
ms_view_clusters(X(1:2,:),labels_mex);
title(sprintf('isosplit5 (mex): N = %d',size(X,2)));


