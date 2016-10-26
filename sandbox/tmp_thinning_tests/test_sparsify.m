function test_sparsify

close all;

rng(9);

N0=1e4;
%approx_num_parcels=600;
num_noise_dims=0;
eps=1;
%p_opts.max_parcel_diameter=2;

A1.N=N0*2; A1.center=[0,0]; A1.cov=[1,0;0,1];
A2.N=N0/2; A2.center=[4.5,0]; A2.cov=[1,0;0,2.5];
A3.N=N0/2; A3.center=[-5,5]; A3.cov=[2.5,1;1,2.5];
A4.N=ceil(N0/16); A4.center=[0,-8]; A4.cov=[0.5,0;0,0.5];
%AA={A1,A2,A3,A4};
AA={A1,A2,A3,A4};
for j=1:length(AA)
    M=length(AA{j}.center);
    center2=rand(1,M+num_noise_dims)*0;
    center2(1:M)=AA{j}.center;
    AA{j}.center=center2;
    cov2=eye(M+num_noise_dims);
    cov2(1:M,1:M)=AA{j}.cov;
    AA{j}.cov=cov2;
end;
[X,labels]=create_multimodal_nd(AA);
N=size(X,2);

figure; ms_view_clusters(X(1:2,:),labels);
title(sprintf('N = %d',size(X,2)));

% [inds_to_use,weights_thin,map]=sparsify(X,eps);
% X_thin=X(:,inds_to_use);
% labels_thin=labels(inds_to_use);
% figure; ms_view_clusters(X_thin(1:2,:));
% title(sprintf('N = %d',size(X_thin,2)));

X_thin=X;
weights_thin=[];
labels_thin=labels;

tic;
labels_iso2=isosplit2(X_thin);
time1=toc;
figure; ms_view_clusters(X_thin(1:2,:),labels_iso2);

isosplit_opts.isocut_threshold=1.5;
isosplit_opts.bin_width=eps;
isosplit_opts.weights=weights_thin;
fprintf('Using bin width: %g\n',isosplit_opts.bin_width);
tic;
labels_iso=isosplit3(X_thin,isosplit_opts);
time2=toc;
figure; ms_view_clusters(X_thin(1:2,:),labels_iso);
title(sprintf('N = %d',size(X_thin,2)));

figure; plot(1:10,1:10)
title(sprintf('%g, %g',time1,time2));

function [inds_to_use,weights,map]=sparsify(X,eps)
M=size(X,1);
N=size(X,2);
if (N==0)
    mask=[];
    return;
end;
mask=ones(1,N)*(-1);
while 1
    inds_undetermined=find(mask==-1);
    if (length(inds_undetermined)==0)
        break;
    end
    if (length(inds_undetermined)==1)
        mask(inds_undetermined)=1;
        break;
    end;
    nearest_inds=knnsearch(X(:,inds_undetermined)',X(:,inds_undetermined)','K',2)';
    nearest_inds=nearest_inds(2,:);
    nearest_inds=inds_undetermined(nearest_inds);
    dists=sqrt(sum((X(:,inds_undetermined)-X(:,nearest_inds)).^2,1));
    isolated_inds=inds_undetermined(find(dists>=eps));
    if (length(isolated_inds)>0)
        mask(isolated_inds)=1;
    else
        pct=0.1;
        inds_to_exclude=randsample(length(inds_undetermined),min(length(inds_undetermined),ceil(length(inds_undetermined)*pct)));
        inds_to_exclude=inds_undetermined(inds_to_exclude);
        mask(inds_to_exclude)=0;   
    end;
end;
inds_to_use=find(mask);
nearest=knnsearch(X(:,inds_to_use)',X','K',2)';
nearest=nearest(2,:);
weights=1+accumarray(nearest',1,[length(inds_to_use),1])';
map=nearest;

function mask=sparsify2(X,eps)
M=size(X,1);
N=size(X,2);
if (N<20)
    % If we have a small number of points, do the brute force sparsify
    mask=sparsify_dense(X,eps);
    return;
end;

% find the dimensions of largest range
mins=min(X,[],2);
maxs=max(X,[],2);
diams=maxs-mins;
[~,ddd]=max(diams);

if (diams(ddd)<=eps)
    mask=sparsify_small_box(X,eps);
    return;
end;

% split into two pieces with some overlap
% do the left first, then the right
divider=(mins(ddd)+maxs(ddd))/2;
inds_left=find(X(ddd,:)<divider+eps/2);
inds_right=find(X(ddd,:)>=divider-eps/2);
mask_left=sparsify(X(:,inds_left),eps);
mask1=zeros(1,N);
mask1(inds_right)=1;
mask1(inds_left(find(mask_left==0)))=0;
new_inds_right=find(mask1);
mask_right=sparsify(X(:,new_inds_right),eps);
mask=ones(1,N);
mask(inds_left(find(mask_left==0)))=0;
mask(new_inds_right(find(mask_right==0)))=0;

function mask=sparsify_small_box(X,eps)
mask=sparsify_dense(X,eps);

function mask=sparsify_dense(X,eps)
M=size(X,1);
N=size(X,2);
disp([M,N]);
mask=ones(1,N);
dists=zeros(N,N);
for m=1:M
    [A,B]=ndgrid(X(m,:),X(m,:));
    dists=dists+(A-B).^2;
end;
dists=sqrt(dists);
[iii,jjj]=ndgrid(1:N,1:N);
inds=find((dists<eps)&(iii<jjj));
for j=1:length(inds)
    ii=iii(j);
    jj=jjj(j);
    if ((mask(ii))&(mask(jj)))
        if (rand<0.5)
            mask(ii)=0;
        else
            mask(jj)=0;
        end;
    end;
end;

function [L,C]=local_kmeans_sorber(X,k)
%KMEANS Cluster multivariate data using the k-means++ algorithm.
%   [L,C] = kmeans(X,k) produces a 1-by-size(X,2) vector L with one class
%   label per column in X and a size(X,1)-by-k matrix C containing the
%   centers corresponding to each class.

%   Version: 2013-02-08
%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%
%   References:
%   [1] J. B. MacQueen, "Some Methods for Classification and Analysis of 
%       MultiVariate Observations", in Proc. of the fifth Berkeley
%       Symposium on Mathematical Statistics and Probability, L. M. L. Cam
%       and J. Neyman, eds., vol. 1, UC Press, 1967, pp. 281-297.
%   [2] D. Arthur and S. Vassilvitskii, "k-means++: The Advantages of
%       Careful Seeding", Technical Report 2006-13, Stanford InfoLab, 2006.

L = [];
L1 = 0;

while length(unique(L)) ~= k
    
    % The k-means++ initialization.
    C = X(:,1+round(rand*(size(X,2)-1)));
    L = ones(1,size(X,2));
    for i = 2:k
        D = X-C(:,L);
%        D = cumsum(sqrt(dot(D,D,1)));  % orig, seems to be dist (l=1)
        D = cumsum(dot(D,D,1));  % Arthur-Vassilvitskii use dist^2 (l=2)
        if D(end) == 0, C(:,i:k) = X(:,ones(1,k-i+1)); return; end
        C(:,i) = X(:,find(rand < D/D(end),1));
        [~,L] = max(bsxfun(@minus,2*real(C'*X),dot(C,C,1).'));
    end
    
    % The k-means algorithm.
    while any(L ~= L1)
        L1 = L;
        for i = 1:k, l = L==i; C(:,i) = sum(X(:,l),2)/sum(l); end
        [~,L] = max(bsxfun(@minus,2*real(C'*X),dot(C,C,1).'),[],1);
    end
    
end
