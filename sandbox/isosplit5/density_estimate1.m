function densities=density_estimate1(X)
if nargin<1, test_density_estimate1; return; end;

[M,N]=size(X);
densities=zeros(1,N);

target_parcel_size=100;
target_num_parcels=1000;

%[~,parcels]=parcelate1(X,100);
[~,parcels]=parcelate2(X,target_parcel_size,target_num_parcels);
length(parcels)
for k=1:length(parcels)
    inds_k=parcels{k}.indices;
    pts=X(:,inds_k);
    centroid=mean(pts,2);
    dists=sqrt(sum((pts-repmat(centroid,1,length(inds_k))).^2,1));
    [~,sort_inds]=sort(dists);
    if (length(sort_inds)>=10)
        dval=1/dists(sort_inds(10));
    else
        dval=0;
    end;
    densities(inds_k)=dval;
end;

function test_density_estimate1

M=2;
N=1e6;
X=randn(M,N);
X=cat(2,X,randn(M,N/100)*0.1+3);
densities=density_estimate1(X);

figure; ms_view_clusters(X,1+(densities>mean(densities))+2*(densities<0.15*mean(densities)));

norms=sqrt(sum(X.^2,1));
figure; semilogy(norms,densities,'b.');
xlabel('norms');
ylabel('densities');

