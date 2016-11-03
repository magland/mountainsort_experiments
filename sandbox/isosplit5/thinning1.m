% ss is the approximate overall subsampling factor
% cc should be between 0 and 1
% cc=1 means subsample at a uniform factor ss
% cc=0 means subsample at a variable factor to equalize the density

function X_thin=thinning1(X,ss,cc,opts)

if nargin<1, test_thinning1; return; end;
if nargin<4, opts=struct; end;
if ~isfield(opts,'target_parcel_size'), opts.target_parcel_size=100; end;
if ~isfield(opts,'target_num_parcels'), opts.target_num_parcels=1000; end;


[M,N]=size(X);

estdensities=density_estimate1(X,opts);
inds0=find(estdensities==0);
inds1=find(estdensities>0);
estlogdensities=zeros(1,N);
estlogdensities(inds1)=log2(estdensities(inds1));

estlogdensities(inds1)=estlogdensities(inds1)-mean(estlogdensities(inds1));
targetlogdensities(inds1)=estlogdensities(inds1)*cc;
probs=2.^(targetlogdensities-estlogdensities);
for passes=1:10
    probs=min(1,probs/(sum(probs))*N/ss);
end;
probs(inds0)=1; %always use those with estimated density of zero
inds_to_use=find(rand(size(probs))<probs);
X_thin=X(:,inds_to_use);

function test_thinning1

M=2;
N=1e6;
ss=5e2;
cc=0.5;

X=randn(M,N);
X=cat(2,X,randn(M,N/100)*0.1+3);
X_thin=thinning1(X,ss,cc);

figure; ms_view_clusters(X_thin);
