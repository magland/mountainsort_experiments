function test_thinning

rng(6);

close all;
M=50;
N0=1e4;
thinning_factor=10;
thinning_param=0.3;

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

% Perform thinning
ttt=tic;
X_thin=thinning1(X,thinning_factor,thinning_param);
fprintf('Total time for thinning: %g\n',toc(ttt));

ttt=tic;
labels_original=isosplit5(X,struct('verbose',0));
fprintf('Time for isosplit5 of original data: %g\n',toc(ttt));

ttt=tic;
labels_thin=isosplit5(X_thin);
fprintf('Time for isosplit5 of thinned data: %g\n',toc(ttt));

figure;
subplot(4,1,1);
hist(X(1,:),100);
subplot(4,1,2);
ms_view_clusters(X(1:2,:),labels_original);
title(sprintf('Original data (isosplit5 mex): N = %d\n',size(X,2)));
subplot(4,1,3);
hist(X_thin(1,:),100);
subplot(4,1,4);
ms_view_clusters(X_thin(1:2,:),labels_thin);
title(sprintf('Thinned data (isosplit5 mex): N = %d\n',size(X_thin,2)));

if 0 %compare with non-mex
ttt=tic;
labels_thin_nonmex=isosplit2(X_thin);
fprintf('Time for isosplit5 of thinned data: %g\n',toc(ttt));

figure;
ms_view_clusters(X_thin(1:2,:),labels_thin_nonmex);
title(sprintf('Thinned data (isosplit5 nonmex): N = %d\n',size(X_thin,2)));
end;
