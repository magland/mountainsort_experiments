function [ks_normalized,cutpoint]=isocut_A1(X,W,opts)

ks_normalized=0;
cutpoint=0;

if (nargin<1) test_isocut_A1; return; end;

if (nargin<2) W=[]; end;
if (length(W)==0) W=ones(size(X)); end;

if (nargin<3) opts=struct; end;

[~,sort_inds]=sort(X);
X=X(sort_inds);
W=W(sort_inds);

x1=X(1);
x2=X(end);

bin_width=(x2-x1)/5;
while 1
    bin_ints=round(X/bin_width);
    i1=bin_ints(1);
    i2=bin_ints(end);
    bin_centers=(i1:i2)*bin_width;
    counts=accumarray(bin_ints'-i1+1,1,[length(bin_centers),1])';
    weighted_counts=accumarray(bin_ints'-i1+1,W',[length(bin_centers),1])';
    figure;
    
    counts_fit=jisotonic(counts,'updown');
    weighted_counts_fit=jisotonic(weighted_counts,'updown');
    
    ks=compute_ks(counts,counts_fit)*sqrt(sum(counts));
    %ks_w=compute_ks(weighted_counts,weighted_counts_fit)*sqrt(sum(weighted_counts));
    ks_w=compute_ks(weighted_counts,weighted_counts_fit)*sqrt(length(X));
    
    bimodal_fit=estimate_bimodal(weighted_counts);
    
    subplot(1,3,1); bar(bin_centers,counts);
    hold on; plot(bin_centers,counts_fit,'r','LineWidth',4);
    title(sprintf('ks = %g',ks));
    subplot(1,3,2); bar(bin_centers,weighted_counts);
    hold on; plot(bin_centers,weighted_counts_fit,'r','LineWidth',4);
    title(sprintf('ks_w = %g',ks_w));
    subplot(1,3,3); bar(bin_centers,weighted_counts);
    hold on; plot(bin_centers,bimodal_fit,'r','LineWidth',4);
    
    if (length(bin_centers)>2000)
        break;
    end;
    
    bin_width=bin_width/4;
end

function Y=estimate_bimodal(counts)
N=length(counts);
sum_counts=sum(counts);
best_index=1;
best_score=-inf;
Y=zeros(size(counts));
%fig=figure;
for j=1:1:length(counts)
    X1=counts(1:j);
    X2=counts(j:end);
    A1=jisotonic(X1,'updown');
    A2=jisotonic(X2,'updown');
    Y0=zeros(size(counts));
    Y0(1:j)=A1;
    Y0(j:end)=A2;
    score0=sum(counts.*Y0/sum_counts)-sum_counts*logsumexp(Y0/sum_counts);
    if (score0>best_score)
        best_score=score0;
        best_index=j;
        Y=Y0;
    end;
    %plot(1:length(counts),counts,'k',1:length(Y0),Y0,'b',1:length(Y0),Y,'r');
    %title(sprintf('%d/%d (%g)',j,length(counts),best_score));
    %pause(0.1);
end
%close(fig);

function Y=logsumexp(X)
val=max(X(:));
Y=val*log(sum(exp(X-val)));

function [ks,ks_ind]=compute_ks(counts1,counts2)
S1=cumsum(counts1);
S2=cumsum(counts2);
[ks,ks_ind]=max(abs(S1-S2)/sum(counts1));
ks_ind=ks_ind(1);

function investigate_k0_normalization
num_trials=500;
ks=[];
weight_vals=[1,2,1,2,1,2,1,2,16,16];
for j=1:num_trials
    N=1000;
    X=randi([1,2],1,N);
    weights=weight_vals(randi(length(weight_vals),size(X)));
    counts=accumarray(X',weights,[2,1])';
    expected_counts=[sum(weights)*0.5,sum(weights)*0.5];
    ks0=compute_ks(counts,expected_counts);
    %ks0=ks0*sqrt(sum(weights));
    ks0=ks0*sqrt(length(X));
    ks(end+1)=ks0;
    
end;
figure; hist(ks,100);
mean(ks)


function test_isocut_A1
close all;

%investigate_k0_normalization;
%return;

rng(1);

N0=2000;

num_trials=1;

ks=[];
cutpoints=[];
for j=1:num_trials
    X1=randn(1,N0);
    X2=randn(1,N0)+2;
    X3=randn(1,N0)+6;
    X=cat(2,X1,X2,X3);
    
    opts=struct;
    [ks0,cutpoint0]=isocut_A1(X,[],opts);
    
    [inds_thin,weights_thin]=thin(X,300);
    X_thin=X(:,inds_thin);
    [ks0,cutpoint0]=isocut_A1(X_thin,weights_thin,opts);
end;

