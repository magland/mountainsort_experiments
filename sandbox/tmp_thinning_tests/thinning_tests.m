function thinning_tests

close all;

N0=2000;
KNN=100;
desired_num=500;
num_noise_dims=10;

%[X,labels]=create_multimodal_2d({[N0,0,0,1,1,0],[N0/2,4.5,0,1,3,0]});
%opts.K_init=5;

A1.N=N0*2; A1.center=[0,0]; A1.cov=[1,0;0,1];
A2.N=N0/2; A2.center=[4.5,0]; A2.cov=[1,0;0,3];
A3.N=N0/2; A3.center=[-5,5]; A3.cov=[3,1;1,3];
A4.N=N0/2; A4.center=[0,-8]; A4.cov=[0.5,0;0,0.5];
AA={A1,A2,A3,A4};
for j=1:length(AA)
    M=length(AA{j}.center);
    center2=rand(1,M+num_noise_dims)*10;
    center2(1:M)=AA{j}.center;
    AA{j}.center=center2;
    cov2=eye(M+num_noise_dims);
    cov2(1:M,1:M)=AA{j}.cov;
    AA{j}.cov=cov2;
end;
[X,labels]=create_multimodal_nd(AA);
N=size(X,2);

figure; ms_view_clusters(X(1:2,:),labels);

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
inds_to_use=find(rnd<probs);
X2=X(:,inds_to_use);
labels2=labels(inds_to_use);
figure; ms_view_clusters(X2(1:2,:),labels2);
figure; ms_view_clusters(X2([1,3],:),labels2);
figure; ms_view_clusters(X2(1:3,:),labels2);

inds_to_not_use=find(rnd>=probs);
labels2=labels;
labels2(inds_to_not_use)=0;
figure; ms_view_clusters(X(1:2,:),labels2);

pause(0.5);
expected_frac=sum(probs)/N
actual_frac=size(X2,2)/N
desired_frac=desired_num/N

%figure; plot(dists_sorted);
%figure; hist(counts,0:KNN);

function ret=compute_dists(X,Y)
ret=sqrt(sum((X-Y).^2,1));


