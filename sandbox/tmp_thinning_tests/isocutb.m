function [ks_normalized,cutpoint]=isocutb(X,W)

if (nargin<1) test_isocutb; return; end;

if (nargin<2) W=[]; end;
if (length(W)==0) W=ones(size(X)); end;
[counts,bin_centers]=autobin(X,W);

counts_fit=jisotonic(counts,'updown');
ks=compute_ks(counts,counts_fit);
ks_normalized=ks*sqrt(length(X));
counts_resid=counts-counts_fit;
counts_resid_fit=jisotonic(counts_resid,'downup');
[~,ind]=min(counts_resid_fit);
ind=ind(1);

cutpoint=bin_centers(ind);

if 1
ks_normalized
cutpoint
figure;
subplot(1,2,1);
bar(bin_centers,counts);
subplot(1,2,2);
bar(bin_centers,counts_fit);
end;

function [counts,bin_centers]=autobin(X,W)
if (length(W)==0) W=ones(size(X)); end;
med=median(X);
Q1=median(X(X<=med));
Q3=median(X(X>=med));
N_IQR=length(X((Q1<=X)&(X<=Q3)));
IQR=Q3-Q1;
bin_width=IQR/N_IQR*2;
%bin_width=0.001;

minval=min(X(:));
maxval=max(X(:));
iminval=round(minval/bin_width);
imaxval=round(maxval/bin_width);
bin_centers=(iminval:imaxval)*bin_width;

inds=round(X/bin_width)-iminval+1;
accumarray(inds',W);

counts=hist(X,bin_centers);

function ks=compute_ks(counts1,counts2)
S1=cumsum(counts1);
S2=cumsum(counts2);
ks=max(abs(S1-S2))/sum(counts1);

function test_isocutb
close all;

N0=200;

num_trials=1000;

ks=[];
cutpoints=[];
for j=1:num_trials
    X1=randn(1,N0);
    X2=randn(1,N0);
    X=cat(2,X1,X2);
    [ks0,cutpoint0]=isocutb(X);
    ks(end+1)=ks0;
    cutpoints(end+1)=cutpoint0;
end;

figure; hist(ks,num_trials);
%figure; hist(cutpoints,num_trials);
mean(ks)
sqrt(var(ks))
mean(cutpoints)
