function ret=kmdensity(X,opts)

if nargin<1, test_kmdensity; return; end;
if nargin<2, opts=struct; end;
if ~isfield(opts,'Ks'), opts.Ks=[10,10,10]; end;
if ~isfield(opts,'num_kmeans_iterations'), opts.num_kmeans_iterations=0; end;
if ~isfield(opts,'K_nn'), opts.K_nn=4; end;

[M,N]=size(X);
%local_kmeans(X,100,5);
labels=local_kmeans_multistep(X,opts.Ks,opts.num_kmeans_iterations);

ret=zeros(1,N);
KK=max(labels);
centroids=zeros(M,KK);
for k=1:max(labels)
    inds_k=find(labels==k);
    if (length(inds_k)>0)
        centroids(:,k)=mean(X(:,inds_k),2);
    end;
end;

for k=1:max(labels)
    inds_k=find(labels==k);
    if (length(inds_k)>=opts.K_nn)
        dists=sqrt(sum((X(:,inds_k)-repmat(centroids(:,k),1,length(inds_k))).^2,1));
        [~,sort_inds]=sort(dists);
        density=(1/dists(sort_inds(opts.K_nn)))^M;
        ret(inds_k)=density;
    end;
end;

function ret=estimate_density_at_point(X,p)


function labels=local_kmeans_multistep(X,Ks,max_iterations)
if (length(Ks)==1)
    labels=local_kmeans(X,Ks(1),max_iterations);
    return;
end;
[M,N]=size(X);
labels1=local_kmeans(X,Ks(1),max_iterations);
labels=zeros(1,N);
for k=1:Ks(1)
    inds_k=find(labels1==k);
    if (length(inds_k)>0)
        labels2=local_kmeans_multistep(X(:,inds_k),Ks(2:end),max_iterations);
        labels(inds_k)=max(labels)+labels2;
    end;
end

function labels=local_kmeans(X,K,max_iterations)
[M,N]=size(X);
K=min(K,N);
centers=X(:,randsample(N,K));
for pass=1:max_iterations
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
distsqrs=zeros(K,N);
for m=1:M
    distsqrs=distsqrs+(repmat(X(m,:),K,1)-repmat(centers(m,:),N,1)').^2;
end;
[~,min_inds]=min(distsqrs,[],1);
labels=min_inds;

function test_kmdensity

M=10;
N=1e5;

X=randn(M,N);

eps=sqrt(M)*1;

tic;
densities=kmdensity(X);
toc

figure; hist(densities,100);
