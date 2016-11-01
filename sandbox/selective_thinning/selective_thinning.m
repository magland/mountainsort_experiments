function selective_thinning
%test1;
%test2;
%test2b;
%test3;
%test4;
test5;

function test5
close all;
M=10;
N0=1e4;
K=6;
cc=0.3;
ds_factor=50;

X1=randn(M,ceil(N0/4));
X2=randn(M,ceil(N0/4));
X2(1,:)=X2(1,:)+4.5;
X3=randn(M,N0*12);
X3(1,:)=X3(1,:)+10.5;
X4=randn(M,N0*24);
X4(1,:)=X4(1,:)+15;
X=cat(2,X1,X2,X3,X4);
N=size(X,2);

hisosplit_opts.verbose=0;
labels1=hisosplit(X,hisosplit_opts);
figure; ms_view_clusters(X(1:2,:),labels1);
figure; hist(X(1,:),100);

knn=knnsearch(X',X','K',K+1)';
knn=knn(K+1,:);
Y=X(:,knn);
estlogdensities=-log2(sqrt(sum((X-Y).^2,1)))*M;
%[~,sort_inds]=sort(estlogdensities);
newestlogdensities=(estlogdensities-mean(estlogdensities))*cc;
probs=2.^(newestlogdensities-estlogdensities);
probs=probs/(sum(probs))*N/ds_factor;
[mean(probs),max(probs)]
%return;
inds_to_use=find(rand(size(probs))<probs);
X_thin=X(:,inds_to_use);


labels_thin=hisosplit(X_thin);
figure; ms_view_clusters(X_thin(1:2,:),labels_thin);
figure; hist(X_thin(1,:),100);

fprintf('ds factor = %g\n',size(X,2)/size(X_thin,2));

knn_thin=knnsearch(X_thin',X_thin','K',K+1)';
knn_thin=knn_thin(K+1,:);
Y_thin=X_thin(:,knn_thin);
estlogdensities_thin=log2(sqrt(sum((X_thin-Y_thin).^2,1)))*M;

%figure; hist(estlogdensities,100);
%figure; hist(estlogdensities_thin,100);

function test4

close all;

M=8;
N=4*1e3;
K=4;

X=rand(M,N);

knn=knnsearch(X',X','K',K+1)';
knn=knn(K+1,:);
Y=X(:,knn);
logdists=log2(sqrt(sum((X-Y).^2,1)))*M;
estlogdensity=logdists-1.25;

inds=find((max(X,[],1))<=1);
X=X(:,inds);
N=size(X,2);
   
mean(logdists)
sqrt(var(logdists))
log2(1/N)-mean(logdists)
disp([mean(estlogdensity),log2(1/N)])
figure; hist(logdists,100);


function test3
N0=1250;
N1=ceil(N0*0.85);

X=zeros(1,0);
X=cat(2,X,rand(1,N0));
X=cat(2,X,rand(1,N1)+1);
X=cat(2,X,rand(1,N0)+2);
[ds,cp]=isocut4_mat(X);

figure; hist(X,100);
hold on;
vline(cp);
title(sprintf('dip score = %g',ds));

function test2
close all;
M=4;
N0=1e3;
K=6;
densities=[1,2,4,8,16,32,64,128];
num_columns=4;
%desired_density=N0*32;
desired_frac=0.5;

X=zeros(M,0);
iii=[1];
for j=1:length(densities)
    density=densities(j);
    X0=rand(M,ceil(N0*density));
    column0=mod(j-1,num_columns);
    row0=floor((j-1)/num_columns);
    X0(1,:)=X0(1,:)+column0*2;
    X0(2,:)=X0(2,:)-row0*2;
    X=cat(2,X,X0);
    iii(end+1)=size(X,2);
end;
N=size(X,2);
desired_N=N*desired_frac;
figure; ms_view_clusters(X(1:2,:));

knn=knnsearch(X',X','K',K+1)';
knn=knn(K+1,:);
identity=1:N;
Y=X(:,knn);
dists=(sum((X-Y).^2,1)).^(1/1.5);
est_densities=dists.^(-M);
[~,sort_inds]=sort(est_densities);
probs=sort_inds/N
probs(find(knn>identity))=0;
probs=min(1,probs);
avg_probs=mean(probs);
probs=probs/avg_probs*desired_N/N;
rnd=rand(size(dists));
to_use=find(rnd<probs);
not_to_use=find(rnd>=probs);
X_thin=X(:,to_use);
figure; ms_view_clusters(X_thin(1:2,:));
%X_thin_complement=X(:,not_to_use);
%figure; ms_view_clusters(X_thin_complement(1:2,:));
frac=size(X_thin,2)/size(X,2)

pcts=[];
for j=1:length(densities)
    pcts(end+1)=length(find((iii(j)<=to_use)&(to_use<=iii(j+1))))/(iii(j+1)-iii(j)+1);
end;
figure; plot(densities,densities,'k',densities,pcts.*densities,'r');

function test2b
close all;
M=6;
N0=1e3;
K=6;
densities=[1,2,4,8,16,32,64,128];
num_columns=4;
%desired_density=N0*32;
desired_frac=0.2;

X1=randn(M,ceil(N0/4));
X2=randn(M,ceil(N0/4));
X2(1,:)=X2(1,:)+4.5;
X3=randn(M,N0*12);
X3(1,:)=X3(1,:)+10.5;
X4=randn(M,N0*24);
X4(1,:)=X4(1,:)+15;
X=cat(2,X1,X2,X3,X4);
N=size(X,2);
desired_N=N*desired_frac;

knn=knnsearch(X',X','K',K+1)';
knn=knn(K+1,:);
identity=1:N;
Y=X(:,knn);
dists=(sum((X-Y).^2,1)).^(1/2);
est_densities=dists.^(-M);

[~,sort_inds]=sort(est_densities);
probs=(1-(sort_inds/N)).^8;

%probs(find(knn>identity))=0;

%return;

avg_probs=mean(probs);
probs=probs/avg_probs*desired_N/N;
probs=min(1,probs);

rnd=rand(size(dists));
to_use=find(rnd<probs);
X_thin=X(:,to_use);

frac=size(X_thin,2)/size(X,2)

to_use2=find(rand(1,N)<frac);
X_thin2=X(:,to_use2);

%figure; ms_view_clusters(X(1:2,:));
%figure; ms_view_clusters(X_thin(1:2,:));
%figure; ms_view_clusters(X_thin2(1:2,:));

figure; hist(X(1,:),100);
figure; hist(X_thin(1,:),100);
figure; hist(X_thin2(1,:),100);

%labels1=hisosplit(X);
labels2=isosplit2(X_thin);
labels3=isosplit2(X_thin2);

%figure; ms_view_clusters(X(1:2,:),labels1);
figure; ms_view_clusters(X_thin(1:2,:),labels2);
figure; ms_view_clusters(X_thin2(1:2,:),labels3);

figure; hist(probs,100);


function test1

close all;

M=4;
N=1e4;
K=3;
densities=[1,2,4,8,16,32,64];
C=0.3471;

data={};
data_thinned={};
for j=1:length(densities)
    density=densities(j);
    data{end+1}=rand(M,N)/((density/N)^(1/M));
end;

means=zeros(1,length(densities));
stdevs=zeros(1,length(densities));
to_use=zeros(1,length(densities));
for j=1:length(densities)
    density=densities(j);
    X=data{j};
    knn=knnsearch(X',X','K',K+1)';
    knn=knn(K+1,:);
    Y=X(:,knn);
    dists=sqrt(sum((X-Y).^2,1));
    est_densities=dists.^(-M)/C;
    means(j)=mean(est_densities);
    stdevs(j)=sqrt(var(est_densities));
    
    probs=1./est_densities;
    probs=max(0,min(1,probs));
    to_use=find(rand(size(dists))<probs);
    %num_used(j)=length(to_use)/length(dists);
    num_used(j)=mean(probs);
end;

figure;
errorbar(densities,means,stdevs);
%set(gca,'xscale','log','yscale','log')

figure;
plot(densities,num_used);

means(1)

