function [inds,weights]=thin(X,desired_num)

KNN=100;
N=size(X,2);

nearest_inds=knnsearch(X',X','K',KNN+1);
nearest_inds=nearest_inds(:,2:end);
dists=zeros(KNN,N);
for k=1:KNN
    dists(k,:)=compute_dists(X,X(:,nearest_inds(:,k)));
end;
dists_sorted=sort(dists(:));
figure; plot(dists_sorted);

for aaa=0.05:0.05:0.95
    cutoff=dists_sorted(ceil(length(dists_sorted)*aaa));
    counts=sum((dists<=cutoff),1);
    probs=1./max(1,counts);
    expected_num=sum(probs);
    if (expected_num<desired_num)
        break;
    end;
end

rnd=rand(1,N);
inds=find(rnd<probs);
weights=1./probs(inds);

expected_frac=sum(probs)/N
actual_frac=length(inds)/N
desired_frac=desired_num/N

%figure; plot(dists_sorted);
%figure; hist(counts,0:KNN);

function ret=compute_dists(X,Y)
ret=sqrt(sum((X-Y).^2,1));