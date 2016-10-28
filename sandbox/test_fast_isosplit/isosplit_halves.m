function [labels,output]=isosplit_halves(X,opts)

if (nargin<1)
    test_isosplit_halves;
	return;
end;

if (nargin<2)
    opts=struct;
end;

if (isfield(opts,'weights'))
    weights=opts.weights;
else
    weights=[];
end;
if (length(weights)==0)
    weights=ones(1,size(X,2));
end;

if (~isfield(opts,'return_iterations'))
    opts.return_iterations=0;
end
if (opts.return_iterations)
    output.iterations={};
end;

if (~isfield(opts,'isocut_threshold')) opts.isocut_threshold=1.2; end;
if (~isfield(opts,'max_pts_in_initial_cluster')) opts.max_pts_in_initial_cluster=100; end;
if (~isfield(opts,'verbose')) opts.verbose=0; end;
if (~isfield(opts,'verbose3')) opts.verbose3=0; end;
if (~isfield(opts,'whiten_at_each_comparison')) opts.whiten_at_each_comparison=1; end;
if (~isfield(opts,'repeat_tolerance')) opts.repeat_tolerance=0.2; end;
if (~isfield(opts,'max_iterations')) opts.max_iterations=inf; end;

if numel(size(X))~=2, error('X must be a 2D array'); end
[M,N]=size(X);

output.num_iterations=0;

if N==0
    labels=[];
    return;
end;

if (N<=opts.max_pts_in_initial_cluster)
    labels=ones(1,N);
    return;
end;

labels_init=random_split(X,3);

inds_01=find(labels_init==1);
inds_02=find(labels_init==2);
inds_03=find(labels_init==3);

[labels_01,output_01]=isosplit_halves(X(:,inds_01),opts);
[labels_02,output_02]=isosplit_halves(X(:,inds_02),opts);
[labels_03,output_03]=isosplit_halves(X(:,inds_03),opts);

if (opts.return_iterations)
    for j=1:length(output_01.iterations)
        output.iterations{end+1}=output_01.iterations{j};
    end;
    for j=1:length(output_02.iterations)
        output.iterations{end+1}=output_02.iterations{j};
    end;
    for j=1:length(output_03.iterations)
        output.iterations{end+1}=output_03.iterations{j};
    end;
end;

K_01=max(labels_01);
K_02=max(labels_02);
K_03=max(labels_03);
labels_init(inds_01)=labels_01;
labels_init(inds_02)=K_01+labels_02;
labels_init(inds_03)=K_01+K_02+labels_03;

K_init=max(labels_init);

labels=labels_init;
active_labels=ones(1,K_init);

centers=zeros(M,K_init);
weighted_counts=zeros(1,K_init);
for k=1:K_init
    inds_k=find(labels==k);
    centers(:,k)=compute_cluster_center(X(:,inds_k),weights(inds_k));
    weighted_counts(k)=sum(weights(inds_k));
end;

attempted_comparisons.centers1=zeros(M,0);
attempted_comparisons.centers2=zeros(M,0);
attempted_comparisons.weighted_counts1=[];
attempted_comparisons.weighted_counts2=[];

while 1
    output.num_iterations=output.num_iterations+1;
    if (output.num_iterations>opts.max_iterations)
        break;
    end;
    [k1,k2]=find_next_comparison(active_labels,centers,weighted_counts,attempted_comparisons,opts.repeat_tolerance);
    if (k1<0) break; end;
    if (opts.return_iterations)
        iteration_info=struct;
        iteration_info.X=X;
        iteration_info.labels_before=labels;
        iteration_info.k1=k1;
        iteration_info.k2=k2;
    end;
    if (opts.verbose)
        fprintf('Iteration %d: Comparing %d (%g) with %d (%g)...',output.num_iterations,k1,weighted_counts(k1),k2,weighted_counts(k2));
    end;
    inds1=find(labels==k1);
    inds2=find(labels==k2);
    inds12=[inds1,inds2];
    attempted_comparisons.centers1=cat(2,attempted_comparisons.centers1,centers(:,k1));
    attempted_comparisons.centers2=cat(2,attempted_comparisons.centers2,centers(:,k2));
    attempted_comparisons.weighted_counts1(end+1)=sum(weights(inds1));
    attempted_comparisons.weighted_counts2(end+1)=sum(weights(inds2));
    attempted_comparisons.centers1=cat(2,attempted_comparisons.centers1,centers(:,k2));
    attempted_comparisons.centers2=cat(2,attempted_comparisons.centers2,centers(:,k1));
    attempted_comparisons.weighted_counts1(end+1)=sum(weights(inds2));
    attempted_comparisons.weighted_counts2(end+1)=sum(weights(inds1));
    [do_merge,labels0,info0]=test_redistribute(X(:,inds1),weights(inds1),X(:,inds2),weights(inds2),opts);
    if (opts.return_iterations)
        iteration_info.projection=info0.projection;
        iteration_info.projection_weights=info0.projection_weights;
        iteration_info.projection_cutpoint=info0.cutpoint;
        iteration_info.projection_labels=info0.labels;
        iteration_info.dip_score=info0.dip_score;
    end
    if (do_merge)||(max(labels0)==1)
        if (opts.verbose)
            fprintf('Merging %d (%g) %d (%g)\n',k1,weighted_counts(k1),k2,weighted_counts(k2));
        end;
        labels(find(labels==k2))=k1;
        centers(:,k1)=compute_cluster_center(X(:,inds12),weights(inds12));
        weighted_counts(k1)=sum(weights(inds12));
        weighted_counts(k2)=0;
        active_labels(k2)=0;
    else
        indsA=inds12(find(labels0==1));
        indsB=inds12(find(labels0==2));
        if (opts.verbose)
            fprintf('Redistributing (%d) (%d)\n',length(indsA),length(indsB));
        end;
        labels(indsA)=k1;
        labels(indsB)=k2;
        centers(:,k1)=compute_cluster_center(X(:,indsA),weights(indsA));
        centers(:,k2)=compute_cluster_center(X(:,indsB),weights(indsB));
        weighted_counts(k1)=sum(weights(indsA));
        weighted_counts(k2)=sum(weights(indsB));
    end;
    if (opts.return_iterations)
        iteration_info.labels=labels;
        output.iterations{end+1}=iteration_info;
    end
end;

labels_map=zeros(1,K_init);
kk=1;
for j=1:K_init
    if (active_labels(j))
        labels_map(j)=kk; kk=kk+1;
    end;
end;
labels=labels_map(labels);

function dists=all_dists_between_centers(centers)
[M,N]=size(centers);
dists=zeros(N,N);
for m=1:M
    [xxx,yyy]=ndgrid(centers(m,:),centers(m,:));
    dists=dists+(xxx-yyy).^2;
end;
dists=sqrt(dists);

function [k1,k2]=find_next_comparison(active_labels,centers,weighted_counts,attempted_comparisons,repeat_tolerance)
active_inds=find(active_labels);
centers_active=centers(:,active_inds);
weighted_counts_active=weighted_counts(active_inds);
Na=size(centers_active,2);
dists=all_dists_between_centers(centers_active);
for j=1:Na
    dists(j,j)=inf;
end
[~,ii]=sort(dists(:));
for j=1:length(ii)
    if (isinf(dists(ii(j))))
       k1=-1; k2=-1;
       return;
    end;
    [k1,k2]=ind2sub(size(dists),ii(j));
    if ((weighted_counts(active_inds(k1))>0)&&(weighted_counts(active_inds(k2))>0)) % just to make sure (this was actually happening! probably should track down why)
        if (~was_already_attempted(attempted_comparisons,centers_active(:,k1),centers_active(:,k2),weighted_counts_active(k1),weighted_counts_active(k2),repeat_tolerance))
            k1=active_inds(k1);
            k2=active_inds(k2);
            return;
        end;
    end;
end;
k1=-1; k2=-1;

function ret=was_already_attempted(attempted_comparisons,center1,center2,weighted_count1,weighted_count2,repeat_tolerance)
AC=attempted_comparisons;
tol=repeat_tolerance;
ii=find( ...
    (abs(AC.weighted_counts1-weighted_count1)<=tol*sqrt((AC.weighted_counts1+weighted_count1)/2)) & ...
    (abs(AC.weighted_counts2-weighted_count2)<=tol*sqrt((AC.weighted_counts2+weighted_count2)/2)) ...
    );

AC_centers1=AC.centers1(:,ii);
AC_centers2=AC.centers2(:,ii);

dists0=sqrt( sum((AC_centers1-AC_centers2).^2,1) );
dists1=sqrt( sum((AC_centers1-repmat(center1,1,length(ii))).^2,1) );
dists2=sqrt( sum((AC_centers2-repmat(center2,1,length(ii))).^2,1) );
aa=find(dists0>0);
dists0=dists0(aa); dists1=dists1(aa); dists2=dists2(aa);

fracs1=dists1./dists0;
fracs2=dists2./dists0;

jj=find( (fracs1<=tol*1/sqrt(weighted_count1)) & (fracs2<=tol*1/sqrt(weighted_count2)) );

ret=~isempty(jj);

function [X1b,X2b,V]=whiten_two_clusters(X1,weights1,X2,weights2)
M=size(X1,1);
N1=size(X1,2);
N2=size(X2,2);

% Important to subtract the two centroids before whitening!
centroid1=compute_cluster_centroid(X1,weights1);
centroid2=compute_cluster_centroid(X2,weights2);
Y1=X1-repmat(centroid1,1,N1);
Y2=X2-repmat(centroid2,1,N2);

% Combine the data
Y=cat(2,Y1,Y2);
N=N1+N2;

% Obtain the whitening matrix using svd
if (N>=M)
    % TODO: Incorporate weights into the whitening!!!
    [U,D,V] = svd(Y,'econ');
    D(D~=0)=1./D(D~=0); %% Problem here if one of the singular values is close to zero!!!!
    % Amd apply it to the original (non-mean subtracted) data
    X1b=sqrt(N-1)*U*D(1:M,1:M)*(U'*X1);
    X2b=sqrt(N-1)*U*D(1:M,1:M)*(U'*X2);
else
    %too few points to whiten
    X1b=X1;
    X2b=X2;
end;

% The best direction is now the one connecting the centroids.
centroid1b=compute_cluster_centroid(X1b,weights1);
centroid2b=compute_cluster_centroid(X2b,weights2);
V=centroid2b-centroid1b;

function [do_merge,labels,info]=test_redistribute(X1,weights1,X2,weights2,opts)
% There was some numerical instability in the whitening -- one of the singular values was close to zero (or was equal to zero) -- fix this
%if opts.whiten_at_each_comparison
%    [X1,X2,V]=whiten_two_clusters(X1,weights1,X2,weights2);
%else
    V=compute_cluster_center(X2,weights2)-compute_cluster_center(X1,weights1);
%end;

if (sum(V.^2)==0)
	warning('isosplit: vector V is null.');
else
    V=V/sqrt(sum(V.^2));
end;
XX=V'*cat(2,X1,X2); %Project onto the line connecting the centroids
WW=cat(2,weights1,weights2);
info.projection=XX;
info.projection_weights=WW;
N=length(XX);
if (N<=5) %avoid a crash - 2/22/2016 jfm
    do_merge=1;
    labels=ones(1,N);    
    info.cutpoint=0;
    info.labels=labels;
    info.dip_score=0;
    return;
end;
%XXs=sort(XX);
isocut4_opts.try_ranges=1;
[dip_score,cutpoint]=isocut4_matlab(XX,WW,isocut4_opts); %This is the core procedure -- split based on isotonic regression
if (dip_score>opts.isocut_threshold)
	%It was a statistically significant split -- so let's redistribute!
	ii1=find(XX<=cutpoint);
	ii2=find(XX>cutpoint);
    do_merge=0;
else
	ii1=1:N;
	ii2=[];
    do_merge=1;
end;
labels=zeros(1,N);
labels(ii1)=1;
labels(ii2)=2;
info.cutpoint=cutpoint;
info.labels=labels;
info.dip_score=dip_score;

if (opts.verbose3)
    %% TODO: this plot should involve weights
    figure;
    bins=linspace(min(XX),max(XX),100);
    if (length(ii2)>0)
        histogram(XX(find(XX<=cutpoint)),bins,'FaceColor','b','EdgeColor','b'); hold on;
        histogram(XX(find(XX>cutpoint)),bins,'FaceColor','r','EdgeColor','r');
    else
        histogram(XX,bins,'FaceColor','k','EdgeColor','k');
    end;
end;

function center=compute_cluster_center(X,weights)
center=weighted_geometric_median(X,weights);

function Y=compute_cluster_centroid(X,weights)
Y=sum(X.*repmat(weights,size(X,1),1),2)/sum(weights);

function [mm,changes]=weighted_geometric_median(X,multiplicities,num_iterations)
if nargin<3, num_iterations=10; end;
[M,N]=size(X);
if (N==1) 
    mm=X; changes=[0];
    return;
end;
weights=ones(1,N);
changes=[];
for it=1:num_iterations
    weights=weights/sum(weights.*multiplicities);
    mm=X*(weights.*multiplicities)';
    if (it>1)
        changes=[changes,sqrt(sum((mm-mm_old).^2))];
    end;
    mm_old=mm;
    diffs=X-repmat(mm,1,N);
    weights=sqrt(sum(diffs.^2,1));
    inds=find(weights~=0);
    weights(inds)=1./weights(inds);
end;

function test_isosplit_halves

close all;

rng(6);

N0=1e4;
%approx_num_parcels=600;
num_noise_dims=0;
p_opts.max_parcel_diameter=2;
bin_width=1;

A1.N=N0*2; A1.center=[0,0]; A1.cov=[1,0;0,1];
A2.N=N0/2; A2.center=[4.5,0]; A2.cov=[1,0;0,2.5];
%A3.N=N0/2; A3.center=[-5,5]; A3.cov=[2.5,1;1,2.5];
%A4.N=ceil(N0/16); A4.center=[0,-8]; A4.cov=[0.5,0;0,0.5];
A3.N=1000; A3.center=[-5,5]; A3.cov=[2.5,1;1,2.5];
A4.N=100; A4.center=[0,-8]; A4.cov=[0.5,0;0,0.5];

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
[samples,true_labels]=create_multimodal_nd(AA);
N=length(samples);

figure;
ms_view_clusters(samples,true_labels);
title(sprintf('N = %d',N));
drawnow;

tic;
[labels,info]=isosplit2(samples);
fprintf('isosplit2: num clusters = %d, num iterations = %d, elapsed = %g\n',max(labels),info.num_iterations,toc);
figure;
ms_view_clusters(samples,labels); title(sprintf('isosplit2: N = %d',N));
drawnow;

tic;
opts.return_iterations=1;
[labels,info]=isosplit_halves(samples,opts);
fprintf('isosplit halves: num clusters = %d, num iterations = %d, elapsed = %g\n',max(labels),info.num_iterations,toc);
figure;
ms_view_clusters(samples,labels); title(sprintf('isosplit halves: N = %d',N));
drawnow;

%show_iterations(samples,info);

%return;

[inds_thin,weights_thin]=thin4(samples,ceil(N/10));
samples_thin=samples(:,inds_thin);

tic;
opts.return_iterations=1;
opts.weights=weights_thin;
[labels_thin,info]=isosplit_halves(samples_thin,opts);
fprintf('thin isosplit halves: num clusters = %d, num iterations = %d, elapsed = %g\n',max(labels),info.num_iterations,toc);
figure;
ms_view_clusters(samples_thin,labels_thin); title(sprintf('thin isosplit halves: N=%d',length(samples_thin)));
drawnow;

%show_iterations(samples_thin,info);

tic;
[labels_thin,info]=isosplit_halves(samples_thin);
fprintf('thin isosplit halves - unweighted: num clusters = %d, num iterations = %d, elapsed = %g\n',max(labels),info.num_iterations,toc);
figure;
ms_view_clusters(samples_thin,labels_thin); title(sprintf('thin isosplit halves unweighted: N=%d',length(samples_thin)));
drawnow;

function show_iterations(X,info)
for j=1:length(info.iterations)
    tmp=info.iterations{j};
    f=figure; set(f,'position',[100,100,2000,600]);
    subplot(1,2,1);
        ms_view_clusters(tmp.X,tmp.labels_before);
        set(gca,'xtick',[]); set(gca,'ytick',[]); xlabel(''); ylabel('');
        title(sprintf('Iteration %d',j), 'FontSize', 20);
    subplot(1,2,2);
        view_clusters_1d_w(tmp.projection,tmp.projection_weights,tmp.projection_labels);
        if (max(tmp.projection_labels)>1)
            vline0(tmp.projection_cutpoint,'k-');
        end;
        title(sprintf('Compare %d,%d, dip score = %g',min(tmp.k1,tmp.k2),max(tmp.k1,tmp.k2),tmp.dip_score), 'FontSize', 20);
    %wait_for_key_press;
    %close(f);
end;

function view_clusters_1d(X,labels)
K=max(labels);
[~,bins]=hist(X,100);
for k=1:K
    inds=find(labels==k);
    if (length(inds)>0)
        vals=hist(X(inds),bins);
        if (k==1) col='k';
        else col='r';
        end;
        hh=bar(bins,vals); hold on;
        set(hh,'FaceColor',col,'EdgeColor',col,'BarWidth',1);
    end;
end;

set(gca,'xtick',[]); set(gca,'ytick',[]); xlabel(''); ylabel('');

function view_clusters_1d_w(X,W,labels)
K=max(labels);
bin_width=(max(X(:))-min(X(:)))/100;
bin_ints=round(X/bin_width);
i1=min(bin_ints);
i2=max(bin_ints);
bins=(i1:i2)*bin_width;
for k=1:K
    inds=find(labels==k);
    if (length(inds)>0)
        if (k==1) col='k';
        else col='r';
        end;
        ii=bin_ints(inds)-i1+1;
        vals=accumarray(ii',W(inds)',[length(bins),1]);
        hh=bar(bins,vals); hold on;
        set(hh,'FaceColor',col,'EdgeColor',col,'BarWidth',1);
    end;
end;

set(gca,'xtick',[]); set(gca,'ytick',[]); xlabel(''); ylabel('');

function vline0(x,linespec)
ylim0=get(gca,'ylim');
plot([x,x],[ylim0(1),ylim0(2)],linespec);

function [X,labels]=create_multimodal_nd(A)
M=length(A{1}.center);
X=zeros(M,0);
labels=zeros(1,0);
for j=1:length(A)
    A0=A{j};
    tmp=randn(M,A0.N);
    %b=[A0(4),A0(6);-A0(6),A0(5)];
    tmp=A0.cov*tmp;
    for m=1:M
        tmp(m,:)=A0.center(m)+tmp(m,:);
    end;
    X=cat(2,X,tmp);
    labels=cat(2,labels,j*ones(1,A0.N));
end

function labels=random_split(X,K)
[M,N]=size(X);
V=randn(1,M);
V=V/sqrt(V*V');
projection=V*X;
labels=zeros(1,N);
sorted=sort(projection);
cutoff_prev=-inf;
for k=1:K
    cutoff_k=sorted(ceil(k/K*N));
    inds_k=find((projection>cutoff_prev)&(projection<=cutoff_k));
    if (length(inds_k)==0)
        disp(k);
        disp(K);
        disp(N);
        disp(cutoff_k);
        error('Unexpected problem in random split');
    end;
    labels(inds_k)=k;
    cutoff_prev=cutoff_k;
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