function test_isosplit_processor

close all;

N0=3000;
%KNN=100;
num_noise_dims=0;

A1.N=N0*4; A1.center=[0,0]; A1.cov=[1,0;0,1];
A2.N=N0/2; A2.center=[5.5,0]; A2.cov=[1,0;0,3];
A3.N=N0/2; A3.center=[-6,5]; A3.cov=[3,1;1,3];
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

desired_num=N/20;

%fig0; ms_view_clusters(X(1:2,:),labels);

tic;
%labels2=isosplit2(X);
%toc
%fig0; ms_view_clusters(X(1:2,:),labels2);

tic;
%labels3=isosplit2_wrapper(X);
%toc
%fig0; ms_view_clusters(X(1:2,:),labels3);

[inds_thin,W_thin]=thin2(X,desired_num);
X_thin=X(:,inds_thin);
tic;
labels5=isosplit2_wrapper(X_thin);
toc
fig0; ms_view_clusters(X_thin(1:2,:),labels5);
%figure; hist(W_thin,200);
title(sprintf('%g',size(X_thin,2)/size(X,2)));

%figure; plot3(X_thin(1,:),X_thin(2,:),W_thin,'r.');


% W=ones(1,N);
% %W=ones(1,N)+rand(1,N)*5;
% tic;
% labels4=isosplit2_w_wrapper(X,W);
% toc
% fig0; ms_view_clusters(X(1:2,:),labels4);

% if density <10 thinning=1
% if density < 20 thinning = 2

function labels=isosplit2_wrapper(X)
X_path=make_temp_mda_path;
labels_path=make_temp_mda_path;
writemda32(X,X_path);
cmd=sprintf('%s isosplit2 --data=%s --labels=%s',mscmd_exe,X_path,labels_path);
disp(cmd);
system(cmd);
labels=readmda(labels_path);
delete(X_path);
delete(labels_path);

function labels=isosplit2_w_wrapper(X,W)
X_path=make_temp_mda_path;
W_path=make_temp_mda_path;
labels_path=make_temp_mda_path;
writemda32(X,X_path);
writemda32(W,W_path);
cmd=sprintf('%s isosplit2_w --data=%s --weights=%s --labels=%s',mscmd_exe,X_path,W_path,labels_path);
disp(cmd);
system(cmd);
labels=readmda(labels_path);
delete(X_path);
delete(W_path);
delete(labels_path);

function fname=make_temp_mda_path
dir = [tempdir,'/mountainlab/tmp_short_term'];
if ~exist(dir,'dir'), mkdir(dir); end     % note can handle creation of parents
fname = [dir,'/',num2str(randi(1e15)),'.mda'];  % random filename


function fig0
figure;
