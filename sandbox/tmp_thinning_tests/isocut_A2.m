function isocut_A2(X,W,opts)

if (nargin<1) test_isocut_A2; return; end;

if (nargin<2) W=[]; end;
if (length(W)==0) W=ones(size(X)); end;

if (nargin<3) opts=struct; end;

[~,sort_inds]=sort(X);
X=X(sort_inds);
W=W(sort_inds);

x1=X(1);
x2=X(end);

bin_width=(x2-x1)/1000;
while 1
    bin_ints=round(X/bin_width);
    i1=bin_ints(1);
    i2=bin_ints(end);
    bin_centers=(i1:i2)*bin_width;
    counts=accumarray(bin_ints'-i1+1,1,[length(bin_centers),1])';
    weighted_counts=accumarray(bin_ints'-i1+1,W',[length(bin_centers),1])';
    
    nbins=length(bin_centers);
    updown_scores=[];
    downup_scores=[];
    j1s=[];
    j2s=[];
    for j1=1:3:nbins
        for j2=j1+1:3:nbins
            wc=weighted_counts(j1:j2);
            updown=jisotonic(wc,'updown');
            downup=jisotonic(wc,'downup');
            updown_ks=compute_ks(wc,updown)*sqrt(sum(counts(j1:j2)));
            downup_ks=compute_ks(wc,downup)*sqrt(sum(counts(j1:j2)));
            updown_scores(end+1)=updown_ks;
            downup_scores(end+1)=downup_ks;
            j1s(end+1)=j1;
            j2s(end+1)=j2;
        end;
    end;
    
    figure; plot(updown_scores,downup_scores,'b.');
    xlabel('updown ks');
    ylabel('downup ks');
    
    figure; plot(j1s,updown_scores-downup_scores,'b.');
    [~,ind]=max(updown_scores-downup_scores);
    disp([j1s(ind),j2s(ind)]);
    wc=weighted_counts(j1s(ind):j2s(ind));
    downup=jisotonic(wc,'downup');
    figure; plot(bin_centers,weighted_counts,'k',bin_centers(j1s(ind):j2s(ind)),downup,'r');
    
    
    %if (length(bin_centers)>2000)
    %    break;
    %end;
    break;
    
    bin_width=bin_width/4;
end



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


function test_isocut_A2
close all;

%investigate_k0_normalization;
%return;

rng(1);

N0=200;

num_trials=1;

ks=[];
cutpoints=[];
for j=1:num_trials
    X1=randn(1,N0);
    X2=randn(1,N0)+2;
    X3=randn(1,N0)*2+9;
    X=cat(2,X1,X2,X3);
    
    opts=struct;
    isocut_A2(X,[],opts);
    
    [inds_thin,weights_thin]=thin(X,300);
    X_thin=X(:,inds_thin);
    isocut_A2(X_thin,weights_thin,opts);
end;

