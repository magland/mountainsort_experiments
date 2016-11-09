function [labels,info]=isosplit5(X,opts)

%% Default parameters and self test
if nargin<1, test_isosplit5; return; end;
if nargin<2, opts=struct; end;
if ~isfield(opts,'isocut_threshold') opts.isocut_threshold=1.5; end;
if ~isfield(opts,'min_cluster_size') opts.min_cluster_size=10; end;
if ~isfield(opts,'K_init') opts.K_init=30; end;
if ~isfield(opts,'refine_clusters'), opts.refine_clusters=true; end;
if ~isfield(opts,'max_iterations'), opts.max_iterations_per_pass=500; end;
if ~isfield(opts,'verbose') opts.verbose=0; end;
if ~isfield(opts,'verbose_pause_duration') opts.verbose_pause_duration=0.5; end;

%% Initialize the timers for diagnostic
timers.get_pairs_to_compare=0;
timers.compare_pairs=0;
timers.compute_centers=0;

[M,N]=size(X);

%% Compute the initial clusters
target_parcel_size=opts.min_cluster_size;
target_num_parcels=opts.K_init;
% !! important not to do a final reassign because then the shapes will not
% be conducive to isosplit iterations -- hexagons are not good for isosplit!
data.labels=parcelate2(X,target_parcel_size,target_num_parcels,struct('final_reassign',0));
Kmax=max(data.labels);

%debug
%labels=data.labels;
%return;

%% Compute the cluster centers
ttt=tic;
data.centers=compute_centers(X,data.labels);
timers.compute_centers=timers.compute_centers+toc(ttt);

%% Repeat while something has been merged in the pass
final_pass=false; % plus we do one final pass at the end
data.comparisons_made=zeros(Kmax,Kmax); % Keep a matrix of comparisons that have been made in this pass
while 1 % Passes
    something_merged=false; % Keep track of whether something has merged in this pass. If not, do a final pass.
    
    clusters_changed_vec=zeros(1,Kmax); % Keep track of the clusters that have changed in this pass so that we can update the comparisons_made matrix at the end
    iteration_number=0;
    while 1 % Iterations
        iteration_number=iteration_number+1;
        if (iteration_number>opts.max_iterations_per_pass)
            error('max iterations per pass exceeded');
            break;
        end;
        % The active labels are those that are still being used
        active_labels_vec=zeros(1,Kmax);
        active_labels_vec(data.labels)=1;
        active_labels=find(active_labels_vec);
        fprintf('num active labels = %d\n',length(active_labels));
        active_centers=data.centers(:,active_labels);

        % Find the pairs to compare on this iteration
        % These will be closest pairs of active clusters that have not yet
        % been compared in this pass
        ttt=tic;
        [inds1,inds2]=get_pairs_to_compare(active_centers,data.comparisons_made(active_labels,active_labels));
        timers.get_pairs_to_compare=timers.get_pairs_to_compare+toc(ttt);

        % If we didn't find any, break from this iteration
        if (length(inds1)==0)
            disp('Nothing else to compare.');
            break;
        end;
        old_labels=data.labels; % So we can determine the number of label changes for diagnostics

        % Actually compare the pairs -- in principle this operation could be parallelized
        ttt=tic;
        [data.labels,clusters_changed]=compare_pairs(X,data.labels,active_labels(inds1),active_labels(inds2),opts);
        clusters_changed_vec(clusters_changed)=1;
        timers.compare_pairs=timers.compare_pairs+toc(ttt);

        % Update which comparisons have been made
        for j=1:length(inds1)
            data.comparisons_made(active_labels(inds1(j)),active_labels(inds2(j)))=1;
            data.comparisons_made(active_labels(inds2(j)),active_labels(inds1(j)))=1;
        end;

        % Recompute the centers -- note: maybe this should only apply to those that changed?
        ttt=tic;
        data.centers=compute_centers(X,data.labels);
        timers.compute_centers=timers.compute_centers+toc(ttt);

        % For diagnostics, cound the number of changes
        total_num_label_changes=length(find(data.labels~=old_labels));
        fprintf('total num label changes = %d\n',total_num_label_changes);
    
        % Determine whether something has merged
        new_active_labels_vec=zeros(1,N);
        new_active_labels_vec(data.labels)=1;
        new_active_labels=find(new_active_labels_vec);
        if (length(new_active_labels)<length(active_labels))
            something_merged=1;
        end;

        % Show the clusters if we are in verbose mode
        if (opts.verbose)&&(length(clusters_changed)>0)
            if (length(new_active_labels)<50)
                labels_map=zeros(1,Kmax);
                active_labels_vec=zeros(1,Kmax);
                active_labels_vec(data.labels)=1;
                active_labels=find(active_labels_vec);
                for ii=1:length(active_labels)
                    labels_map(active_labels(ii))=ii;
                end;
                labels_mapped=labels_map(data.labels);
                figure; ms_view_clusters(X(1:2,:),labels_mapped);
                pause(opts.verbose_pause_duration);
            end;
        end;
        
        %break;
    end;
    % zero out the comparisons made matrix only for those that have changed
    clusters_changed=find(clusters_changed_vec);
    for j=1:length(clusters_changed)
        data.comparisons_made(clusters_changed(j),:)=0;
        data.comparisons_made(:,clusters_changed(j))=0;
    end;
    
    disp('---------------------------------------------');
    disp(clusters_changed);
    
    if (something_merged) final_pass=false; end;
    if (final_pass) break; end; % This was the final pass and nothing has merged
    if (~something_merged) final_pass=true; end; % If we are done, do one last pass for final redistributes
    
    %break;
end;

% This is the result
labels=data.labels;

% But we should remap the labels to occupy the first natural numbers
labels_map=zeros(1,Kmax);
active_labels_vec=zeros(1,Kmax);
active_labels_vec(labels)=1;
active_labels=find(active_labels_vec);
for ii=1:length(active_labels)
    labels_map(active_labels(ii))=ii;
end;
labels=labels_map(labels);

% Return the timers in the info
info.timers=timers;

% If the user wants to refine the clusters, then we repeat isosplit on each
% of the new clusters, recursively. Unless we only found one cluster.
if ((opts.refine_clusters)&&(max(labels)>1))
    opts2=opts;
    opts2.refine_clusters=true; % Maybe we should provide an option on whether to do recursive refinement
    K=max(labels);
    labels_split=zeros(1,N);
    for k=1:K
        inds_k=find(labels==k);
        X_k=X(:,inds_k);
        labels_k=isosplit5(X_k,opts2);
        labels_split(inds_k)=max(labels_split)+labels_k;
    end;
    labels=labels_split;
end;

function centers=compute_centers(X,labels)
[M,N]=size(X);
centers=zeros(M,N);
counts=accumarray(labels',1,[N,1])';
for m=1:M
    centers(m,:)=accumarray(labels',X(m,:)',[N,1])';
end;
centers(:,find(counts))=centers(:,find(counts))./repmat(counts(find(counts)),M,1);

function [new_labels,clusters_changed]=compare_pairs(X,labels,k1s,k2s,opts)
clusters_changed_vec=zeros(1,max(labels));
new_labels=labels;
for i1=1:length(k1s)
    k1=k1s(i1);
    k2=k2s(i1);
    inds1=find(labels==k1);
    inds2=find(labels==k2);
    if ((length(inds1)>0)&&(length(inds2)>0))
        if ((length(inds1)<opts.min_cluster_size)||(length(inds2)<opts.min_cluster_size))
            do_merge=1;
        else
            inds12=cat(2,inds1,inds2);
            L12_old=cat(2,ones(1,length(inds1)),2*ones(1,length(inds2)));
            [do_merge,L12,proj]=merge_test(X(:,inds1),X(:,inds2),opts);
        end;
        if (do_merge)
            new_labels(find(new_labels==k2))=k1;
            clusters_changed_vec(k1)=1;
            clusters_changed_vec(k2)=1;
        else
            %redistribute
            new_labels(inds12(find(L12==1)))=k1;
            new_labels(inds12(find(L12==2)))=k2;
            if (length(find(L12~=L12_old))>0)
                clusters_changed_vec(k1)=1;
                clusters_changed_vec(k2)=1;
            end;
        end;
    end;
end;
clusters_changed=find(clusters_changed_vec);

function [ret,new_labels,projection12]=merge_test(X1,X2,opts)
[~,N1]=size(X1); [~,N2]=size(X2);
if ((N1==0)||(N2==0))
    error('Error in merge test: N1 or N2 is zero');
end;
centroid1=mean(X1,2);
centroid2=mean(X2,2);
V=centroid2-centroid1;
V=V/sqrt(V'*V);
projection1=V'*X1;
projection2=V'*X2;
projection12=cat(2,projection1,projection2);
[dipscore,cutpoint]=isocut5(projection12,ones(size(projection12)));
ret=(dipscore<opts.isocut_threshold);
%cutpoint=isocut(projection12,opts.isocut_threshold);
%ret=(cutpoint~=0);
new_labels=ones(1,N1+N2);
new_labels(find(projection12>=cutpoint))=2;

function dists=make_dists_matrix(centers)
[M,N]=size(centers);
dists=zeros(N,N);
[aa,bb]=ndgrid(1:N,1:N);
for m=1:M
    dists=dists+reshape((centers(m,aa(:))-centers(m,bb(:))).^2,N,N);
end;
dists=sqrt(dists);

function [inds1,inds2]=get_pairs_to_compare(centers,comparisons_made)
[M,N]=size(centers);
inds1=[];
inds2=[];
dists=make_dists_matrix(centers);
dists(find(comparisons_made(:)))=inf;
for j=1:N
    dists(j,j)=inf;
end;
something_changed=1;
while (something_changed)
    something_changed=0;
    [~,best_inds]=min(dists,[],1);
    for j=1:N
        if (best_inds(j)>j)
            if (best_inds(best_inds(j))==j) % mutual
                if (dists(j,best_inds(j))<inf)
                    inds1(end+1)=j;
                    inds2(end+1)=best_inds(j);
                    dists(j,:)=inf;
                    dists(:,j)=inf;
                    dists(best_inds(j),:)=inf;
                    dists(:,best_inds(j))=inf;
                    something_changed=1;
                end;
            end;
        end;        
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,labels]=generate_dataset

N0=1e4;
num_noise_dims=0;

A1.N=N0*2; A1.center=[0,0]; A1.cov=[1,0;0,1];
A2.N=N0/2; A2.center=[6,0]; A2.cov=[1,0;0,2.5];
%A3.N=N0/2; A3.center=[-5,5]; A3.cov=[2.5,1;1,2.5];
%A4.N=ceil(N0/16); A4.center=[0,-8]; A4.cov=[0.5,0;0,0.5];
A3.N=1000; A3.center=[-5,5]; A3.cov=[2.5,1;1,2.5];
A4.N=500; A4.center=[0,-8]; A4.cov=[0.5,0;0,0.5];

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_isosplit5
rng(2);
close all;

[X,true_labels]=generate_dataset;
figure; ms_view_clusters(X(1:2,:),true_labels);
title('Truth');

ttt=tic;
[labels2,info]=isosplit5(X,struct('verbose',1,'refine_clusters',0));
fprintf('Time for isosplit5: %g\n',toc(ttt));
figure; ms_view_clusters(X(1:2,:),labels2);
title('isosplit5');
disp(info.timers);

ttt=tic;
labels_mex=isosplit5_mex(X);
fprintf('Time for isosplit5_mex: %g\n',toc(ttt));
figure; ms_view_clusters(X(1:2,:),labels_mex);
title('isosplit5 mex');

