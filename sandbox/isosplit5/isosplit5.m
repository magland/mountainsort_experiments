function [labels,info]=isosplit5(X,opts)

if nargin<1, test_isosplit5; return; end;
if nargin<2, opts=struct; end;
if ~isfield(opts,'isocut_threshold') opts.isocut_threshold=1.5; end;
if ~isfield(opts,'min_cluster_size') opts.min_cluster_size=10; end;
if ~isfield(opts,'verbose') opts.verbose=0; end;

timers.get_pairs_to_compare=0;
timers.compare_pairs=0;
timers.compute_centers=0;

[M,N]=size(X);
%data.labels=1:N; % each point in its own cluster
data.labels=initialize_labels(X);
Kmax=max(data.labels);

ttt=tic;
data.centers=compute_centers(X,data.labels);
timers.compute_centers=timers.compute_centers+toc(ttt);

max_iterations=500;
max_iterations_without_merges=5;

something_merged=1;
while something_merged
    something_merged=0;
    data.comparisons_made=zeros(Kmax,Kmax);
    iteration_number=1;
    num_iterations_without_merges=0;
    while 1
        iteration_number=iteration_number+1;
        if (iteration_number>max_iterations)
            error('max iterations exceeded');
            break;
        end;
        active_labels_vec=zeros(1,N);
        active_labels_vec(data.labels)=1;
        active_labels=find(active_labels_vec);
        fprintf('num active labels = %d\n',length(active_labels));
        %if (length(active_labels)<=50);
        %    figure; ms_view_clusters(X,data.labels);
        %end;
        active_centers=data.centers(:,active_labels);

        ttt=tic;
        [inds1,inds2]=get_pairs_to_compare(active_centers,data.comparisons_made(active_labels,active_labels));
        timers.get_pairs_to_compare=timers.get_pairs_to_compare+toc(ttt);

        if (length(inds1)==0)
            disp('Nothing else to compare.');
            break;
        end;
        old_labels=data.labels;

        ttt=tic;
        [data.labels,changes]=compare_pairs(X,data.labels,active_labels(inds1),active_labels(inds2),opts);
        timers.compare_pairs=timers.compare_pairs+toc(ttt);

        for j=1:length(inds1)
            data.comparisons_made(active_labels(inds1(j)),active_labels(inds2(j)))=1;
            data.comparisons_made(active_labels(inds2(j)),active_labels(inds1(j)))=1;
        end;

        ttt=tic;
        data.centers=compute_centers(X,data.labels);
        timers.compute_centers=timers.compute_centers+toc(ttt);

        total_num_label_changes=length(find(data.labels~=old_labels));
        fprintf('total num label changes = %d\n',total_num_label_changes);
    %     if (total_num_label_changes<10)
    %         break;
    %     end;
        new_active_labels_vec=zeros(1,N);
        new_active_labels_vec(data.labels)=1;
        new_active_labels=find(new_active_labels_vec);
        if (length(new_active_labels)<length(active_labels))
            something_merged=1;
        end;
%         if (length(new_active_labels)==length(active_labels))
%             num_iterations_without_merges=num_iterations_without_merges+1;
%             if (num_iterations_without_merges>max_iterations_without_merges)
%                 break;
%             end;
%         else
%             num_iterations_without_merges=0;
%         end;

        if (opts.verbose)&&(length(find(changes==1))>0)
            if (length(new_active_labels)<50)
                labels_map=zeros(1,N);
                active_labels_vec=zeros(1,N);
                active_labels_vec(data.labels)=1;
                active_labels=find(active_labels_vec);
                for ii=1:length(active_labels)
                    labels_map(active_labels(ii))=ii;
                end;
                labels_mapped=labels_map(data.labels);
                figure; ms_view_clusters(X(1:2,:),labels_mapped);
                pause(0.03);
            end;
        end;
    end;
end;

for pass=1:1
    active_labels_vec=zeros(1,N);
    active_labels_vec(data.labels)=1;
    active_labels=find(active_labels_vec);
    for i1=1:length(active_labels)
        for i2=i1+1:length(active_labels)
            k1=active_labels(i1);
            k2=active_labels(i2);

            ttt=tic;
            [new_labels,changes]=compare_pairs(X,data.labels,k1,k2,opts);
            timers.compare_pairs=timers.compare_pairs+toc(ttt);
            
            total_num_label_changes=length(find(data.labels~=new_labels));
            fprintf('total num label changes for %d/%d = %d (pass=%d,num_active=%d)\n',k1,k2,total_num_label_changes,pass,length(active_labels));

            data.labels=new_labels;

            if (opts.verbose)&&(length(find(changes==1))>0)
                labels_map=zeros(1,N);
                for ii=1:length(active_labels)
                    labels_map(active_labels(ii))=ii;
                end;
                labels_mapped=labels_map(data.labels);
                figure; ms_view_clusters(X(1:2,:),labels_mapped);
                title(sprintf('k1/k2 = %d/%d (pass %d)',k1,k2,pass));
                pause(0.03);
            end;
        end;
    end;
end;

labels=data.labels;

labels_map=zeros(1,N);
active_labels_vec=zeros(1,N);
active_labels_vec(labels)=1;
active_labels=find(active_labels_vec);
for ii=1:length(active_labels)
    labels_map(active_labels(ii))=ii;
end;
labels=labels_map(labels);

info.timers=timers;

function centers=compute_centers(X,labels)
[M,N]=size(X);
centers=zeros(M,N);
counts=accumarray(labels',1,[N,1])';
for m=1:M
    centers(m,:)=accumarray(labels',X(m,:)',[N,1])';
end;
centers(:,find(counts))=centers(:,find(counts))./repmat(counts(find(counts)),M,1);

function [new_labels,changes]=compare_pairs(X,labels,k1s,k2s,opts)
changes=zeros(1,length(k1s));
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
            changes(i1)=1;
        else
            %redistribute
            new_labels(inds12(find(L12==1)))=k1;
            new_labels(inds12(find(L12==2)))=k2;
            if (length(find(L12~=L12_old))>0)
                changes(i1)=1;
            end;
        end;
    end;
end;

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

function [inds1,inds2]=get_pairs_to_compare_old(centers)
[M,N]=size(centers);
inds1=[];
inds2=[];
while 1
    used_vec=zeros(1,N);
    used_vec([inds1,inds2])=1;
    unused_inds=find(used_vec==0);
    if (length(unused_inds)==0) break; end;
    [i1,i2]=get_pairs_to_compare_2(centers(:,unused_inds));
    if (length(i1)==0) break; end;
    inds1=[inds1,unused_inds(i1)];
    inds2=[inds2,unused_inds(i2)];
end;

function nearest_inds=find_nearest(pts)
[M,N]=size(pts);
nearest=knnsearch(pts',pts','K',2)';
nearest_inds=nearest(2,:);

function [inds1,inds2]=get_pairs_to_compare_2(centers)
[M,N]=size(centers);
if (N<=1)
    inds1=[]; inds2=[];
    return;
end;

nearest_inds=find_nearest(centers);

% the points that have a mutual nearest neighbor
identity=1:N;
aaa=find(nearest_inds); %beware of zeros! meaning that there was no associated nearest that was needed
mutual_inds=find(nearest_inds(nearest_inds(aaa))==identity(aaa));
mutual_inds=aaa(mutual_inds);

% Compare each of these with its mutual neighbor
inds1=mutual_inds;
inds2=nearest_inds(mutual_inds);

% But everything will be counted twice! so only use one of each
iii=find(inds1<inds2);
inds1=inds1(iii);
inds2=inds2(iii);

% Test! (not needed if code is correct)
bad_guys=find(ismember(inds1,inds2));
if (length(bad_guys)>0)
    error('Unexpected problem.');
end;

function labels=initialize_labels(X)
[M,N]=size(X);
labels1=initialize_labels_2(X);
K=max(labels1);
labels=zeros(1,N);
for k=1:K
    inds_k=find(labels1==k);
    labels_k=initialize_labels_2(X(:,inds_k));
    labels(inds_k)=max(labels)+labels_k;
end;

function labels=initialize_labels_2(X)
[M,N]=size(X);
K=30;
K=min(K,N);
centers=X(:,randsample(N,K));
for pass=1:1
    distsqrs=zeros(K,N);
    for m=1:M
        distsqrs=distsqrs+(repmat(X(m,:),K,1)-repmat(centers(m,:),N,1)').^2;
    end;
    [~,min_inds]=min(distsqrs,[],1);
    labels=min_inds;
    counts=accumarray(labels',1,[K,1])';
    sums=zeros(M,K);
    for m=1:M
        sums(m,:)=accumarray(labels',X(m,:)',[K,1])';
    end;
    centers=sums./repmat(counts,M,1);
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
ttt=tic;
[labels2,info]=isosplit5(X,struct('verbose',0));
fprintf('Time for isosplit5: %g\n',toc(ttt));
figure; ms_view_clusters(X(1:2,:),labels2);
title('isosplit5');
disp(info.timers);

ttt=tic;
labels_mex=isosplit5_mex(X);
fprintf('Time for isosplit5_mex: %g\n',toc(ttt));
figure; ms_view_clusters(X(1:2,:),labels_mex);
title('isosplit5 mex');


