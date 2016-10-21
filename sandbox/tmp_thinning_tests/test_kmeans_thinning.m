function test_kmeans_thinning

close all;

rng(2);

N0=1e4;
approx_num_parcels=300;
num_noise_dims=8;

A1.N=N0*2; A1.center=[0,0]; A1.cov=[1,0;0,1];
A2.N=N0/2; A2.center=[4.5,0]; A2.cov=[1,0;0,2.5];
A3.N=N0/2; A3.center=[-5,5]; A3.cov=[2.5,1;1,2.5];
A4.N=ceil(N0/16); A4.center=[0,-8]; A4.cov=[0.5,0;0,0.5];
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

disp('Computing parcelation');
[parcel_labels,parcels]=compute_parcelation(X,approx_num_parcels,2);
disp('Done computing parcelation');

disp('Getting parcel diameters');
[diameters,pointwise_diameters]=get_parcel_diameters(parcels,N);
figure; hist(diameters,40);
title('Parcel diameters');
figure; hist(pointwise_diameters,80);
title('Pointwise parcel diameters');

disp('Getting parcel sizes');
[parcel_sizes,pointwise_sizes]=get_parcel_sizes(parcels,N);
figure; hist(parcel_sizes,40);
title('Parcel sizes');
figure; hist(pointwise_sizes,80);
title('Thinning factors');
disp(sum(1./pointwise_sizes)) % should be around the number of parcels, obviously

disp('Plotting data with parcel labels');
figure; ms_view_clusters(X(1:2,:),parcel_labels);
legend(gca,'off');

disp('Computing parcel representatives');
representative_inds=compute_parcel_representative_inds(X,parcels);
X_thin=X(:,representative_inds);
figure; ms_view_clusters(X_thin(1:2,:),labels(representative_inds));
legend(gca,'off');
title('After thinning');

isosplit_opts=struct;
isosplit_opts.isocut_threshold=1.5;

inds=find(diameters>0);
if (~isempty(inds))
    average_diameter=mean(diameters(inds));
    isosplit_opts.bin_width=average_diameter/6;
    fprintf('Using bin width: %g\n',isosplit_opts.bin_width);
end;

disp('unweighted isosplit3');
isosplit_labels=isosplit3(X_thin,isosplit_opts);
figure; ms_view_clusters(X_thin(1:2,:),isosplit_labels);
title('Unweighted iso-split clustering');

disp('weighted isosplit3');
isosplit_opts.weights=parcel_sizes;
isosplit_opts.diameters=diameters;
isosplit_opts.return_iterations=1;
[isosplit_labels,iso_info]=isosplit3(X_thin,isosplit_opts);
figure; ms_view_clusters(X_thin(1:2,:),isosplit_labels);
title('Weighted iso-split clustering');
drawnow;

show_iterations(X_thin,iso_info);


function ret=compute_parcel_representative_inds(X,parcels)
ret=zeros(1,length(parcels));
for j=1:length(parcels)
    ind0=compute_representative_ind(X(:,parcels{j}.inds));
    ret(j)=parcels{j}.inds(ind0);
end;

function ret=compute_representative_ind(X)
%ret=1;
%return;
N=size(X,2);
centroid=mean(X,2);
dists_to_centroid=sqrt(sum((X-repmat(centroid,1,N)).^2));
[~,closest_ind]=min(dists_to_centroid);
ret=closest_ind(1);

function [labels,parcels]=compute_parcelation(X,num_parcels,K)
N=size(X,2);
parcels={};
parcels{end+1}=create_parcel(X,1:N);
something_changed=1;
while ((length(parcels)<num_parcels)&&(length(parcels)<N)&&(something_changed)) %Better check needed to avoid infinite loop
    disp(length(parcels));
    something_changed=0;
    num_to_split=min(10,length(parcels));
    parcel_inds=select_parcels_to_split(parcels,num_to_split,N);
    for j=1:length(parcel_inds)
        ind0=parcel_inds(j);
        if (length(parcels{ind0}.inds)>1)
            new_parcels=split_parcel(X,parcels{ind0},min(K,length(parcels{ind0}.inds)));
            if (length(new_parcels)>1)
                parcels{ind0}=new_parcels{1};
                for k=2:length(new_parcels)
                    parcels{end+1}=new_parcels{k};
                end;
                something_changed=1;
            end;
        end;
    end;
end
labels=zeros(1,N);
for j=1:length(parcels)
    labels(parcels{j}.inds)=j;
end;

function ret=select_parcels_to_split(parcels,num,N)
diameters=get_parcel_diameters(parcels,N);
[~,inds]=sort(diameters,'descend');
ret=inds(1:num);

function [diameters,pointwise_diameters]=get_parcel_diameters(parcels,N)
diameters=zeros(1,length(parcels));
pointwise_diameters=zeros(1,N);
for j=1:length(parcels)
    diameters(j)=parcels{j}.diameter;
    pointwise_diameters(parcels{j}.inds)=parcels{j}.diameter;
end;

function [sizes,pointwise_sizes]=get_parcel_sizes(parcels,N)
sizes=zeros(1,length(parcels));
pointwise_sizes=zeros(1,N);
for j=1:length(parcels)
    sizes(j)=length(parcels{j}.inds);
    pointwise_sizes(parcels{j}.inds)=length(parcels{j}.inds);
end;

function ret=split_parcel(X,P,K)
ret={};
if (length(P.inds)<=1)
    error('Cannot split parcel with fewer than two points.');
end;

mins=min(X(:,P.inds),[],2);
maxs=max(X(:,P.inds),[],2);
diams=maxs-mins;
[~,iii]=max(diams); iii=iii(1);
cut=(maxs(iii)+mins(iii))/2;
inds1=find(X(iii,P.inds)<cut);
inds2=find(X(iii,P.inds)>=cut);
if (length(inds1)==0)||(length(inds2)==0)
    inds1=1:ceil(length(P.inds)/2);
    inds2=ceil(length(P.inds)/2)+1:length(P.inds);
end;
ret{1}=create_parcel(X,P.inds(inds1));
ret{2}=create_parcel(X,P.inds(inds2));

% Y=X(:,P.inds);
% labels=local_kmeans_sorber(Y,K);
% for k=1:K
%     inds_k=find(labels==k);
%     if (length(inds_k)==0)
%         warning(sprintf('Problem splitting parcel. N=%d, max_k=%d',size(Y,2),max(labels)));
%     end;
%     ret{end+1}=create_parcel(X,P.inds(inds_k));
% end

function show_iterations(X,info)
for j=1:length(info.iterations)
    tmp=info.iterations{j};
    f=figure; set(f,'position',[100,100,2000,600]);
    subplot(1,2,1);
        ms_view_clusters(X,tmp.labels_before);
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



function ret=create_parcel(X,inds)
if (length(inds)==0)
    error('Cannot create parcel. Length of inds is zero');
end;
Y=X(:,inds);
mins=min(Y,[],2);
maxs=max(Y,[],2);
ret.inds=inds;
% Important: need a better estimate of the size
% Probably diameter is best
ret.diameter=mean(maxs-mins); %???????????????????????


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
