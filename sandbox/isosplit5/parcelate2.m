function [labels,parcels]=parcelate2(X,target_parcel_size,target_num_parcels)
if nargin<1, test_parcelate2; return; end;

[M,N]=size(X);
labels=zeros(1,N);

parcels={};

P.indices=1:N;
P.centroid=mean(X(:,P.indices),2);
P.radius=compute_max_distance(P.centroid,X(:,P.indices));
labels(P.indices)=1;
parcels{end+1}=P;

split_factor=3; % split factor around 2.71 is in a sense ideal

target_radius=P.radius;
while (length(parcels)<target_num_parcels)
    if (get_max_parcel_size(parcels)<target_parcel_size)
        % nothing else will ever be split
        break;
    end;
    target_radius=target_radius*0.9;
    
    p_index=1;
    while (p_index<=length(parcels))
        inds=parcels{p_index}.indices;
        rad=parcels{p_index}.radius;
        sz=length(inds);
        if (sz>target_parcel_size)&&(rad>target_radius)
            iii=randsample(sz,split_factor);
            pts=X(:,inds(iii));
            dm=distance_matrix(pts,X(:,inds));
            [~,assignments]=min(dm,[],1);
            parcels{p_index}.indices=inds(find(assignments==1));
            parcels{p_index}.centroid=mean(X(:,parcels{p_index}.indices),2);
            parcels{p_index}.radius=compute_max_distance(parcels{p_index}.centroid,X(:,parcels{p_index}.indices));
            labels(parcels{p_index}.indices)=p_index;
            for jj=2:split_factor
                PP.indices=inds(find(assignments==jj));
                PP.centroid=mean(X(:,PP.indices),2);
                PP.radius=compute_max_distance(PP.centroid,X(:,PP.indices));
                parcels{end+1}=PP;
                labels(PP.indices)=length(parcels);
            end;
        else
            p_index=p_index+1;
        end;
    end;
end;

function ret=get_max_parcel_size(parcels)
tmp=zeros(1,length(parcels));
for j=1:length(parcels)
    tmp(j)=length(parcels{j}.indices);
end;
ret=max(tmp);

function ret=compute_max_distance(pt,X)
dm=distance_matrix(pt,X);
ret=max(dm);

function ret=distance_matrix(X1,X2)
[M,N1]=size(X1);
[M,N2]=size(X2);
ret=zeros(N1,N2);
[aa,bb]=ndgrid(1:N1,1:N2);
for m=1:M
    ret=ret+reshape((X1(m,aa(:))-X2(m,bb(:))).^2,N1,N2);
end;
ret=sqrt(ret);

function test_parcelate2

target_parcel_size=100;
target_num_parcels=500;
M=2;
N=1e6;
X=randn(M,N);
tic;
labels=parcelate2(X,target_parcel_size,target_num_parcels);
toc
figure; ms_view_clusters(X,labels);

