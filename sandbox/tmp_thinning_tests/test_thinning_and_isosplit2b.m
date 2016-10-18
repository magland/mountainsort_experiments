function test_thinning_and_isosplit2b

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

[inds_thin,weights_thin]=thin(X,1000);
X_thin=X(:,inds_thin);
labels_thin=labels(inds_thin);
figure; ms_view_clusters(X_thin(1:2,:),labels_thin);

labels2=isosplit2(X_thin);
figure; ms_view_clusters(X_thin(1:2,:),labels2);
