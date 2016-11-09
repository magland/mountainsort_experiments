function ret=approx_knn(X,K_nn)
[M,N]=size(X);
ret=zeros(1,N);
ttt=tic;
labels=split_into_equal_size_clusters(X,max(2000,K_nn*2),max(4000,K_nn*4));
fprintf('Time for split into equal size clusters: %g\n',toc(ttt));
K=max(labels);
for k=1:K
    inds_k=find(labels==k);
    X_k=X(:,inds_k);
    knn=knnsearch(X_k',X_k','K',K_nn+1)';
    knn=knn(K_nn+1,:);
    ret(inds_k)=inds_k(knn);
end;

function labels=split_into_equal_size_clusters(X,min_cluster_size,max_cluster_size)
[M,N]=size(X);

if (N<=max_cluster_size)
    labels=ones(1,N);
    return;
end;
K=5;
while (K*max_cluster_size>N)&&(K>2)
    K=K-1;
end;

labels0=local_kmeans(X,K);
while 1
    if (K<=1), break; end;
    counts=zeros(1,K);
    for k=1:K
        inds_k=find(labels0==k);
        counts(k)=length(inds_k);
    end;
    if (min(counts)>=min_cluster_size)
        break;
    end;
    %combine the smallest two
    [~,sort_inds]=sort(counts);
    k1=min(sort_inds(1),sort_inds(2));
    k2=max(sort_inds(2),sort_inds(1));
    labels0(find(labels0==k2))=k1;
    labels0(find(labels0>k2))=labels0(find(labels0>k2))-1;
    K=max(labels0);
end;

if (K==1)
    labels=ones(1,N);
    return;
end;

labels=zeros(1,N);
for k=1:K
    inds_k=find(labels0==k);
    if (length(inds_k)>max_cluster_size)
        labels_k=split_into_equal_size_clusters(X(:,inds_k),min_cluster_size,max_cluster_size);
        labels(inds_k)=max(labels)+labels_k;
    else
        labels(inds_k)=max(labels)+1;
    end;
end;

function labels=local_kmeans(X,K,max_iterations)
if nargin<3, max_iterations=2; end;
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

