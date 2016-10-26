function test_thinning2

close all;

rng(6);

N0=1e5;
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
[X,labels]=create_multimodal_nd(AA);
N=size(X,2);

figure; ms_view_clusters(X(1:2,:),labels);

split_into_sections(X);
