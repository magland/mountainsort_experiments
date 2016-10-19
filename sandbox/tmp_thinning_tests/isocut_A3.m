function isocut_A3(X,W,opts)

if (nargin<1) test_isocut_A3; return; end;

if (nargin<2) W=[]; end;
if (length(W)==0) W=ones(size(X)); end;

if (nargin<3) opts=struct; end;

[~,sort_inds]=sort(X);
X=X(sort_inds);
W=W(sort_inds);

x1=X(1);
x2=X(end);

bin_width=(x2-x1)/45;

bin_ints=round(X/bin_width);
i1=bin_ints(1);
i2=bin_ints(end);
bin_centers=(i1:i2)*bin_width;
counts=accumarray(bin_ints'-i1+1,1,[length(bin_centers),1])';
weighted_counts=accumarray(bin_ints'-i1+1,W',[length(bin_centers),1])';

tic;
[score,i0,i1,i2]=compute_dipscore(counts,weighted_counts);
toc

if (opts.verbose)
    figure;
    bar(bin_centers,weighted_counts,'k');
    hold on;
    vline(bin_centers(i1),'b');
    vline(bin_centers(i2),'b');
    vline(bin_centers(i0),'r');
    title(sprintf('dip score = %g (num.pts = %d)\n',score,sum(counts)));
end


function [dipscore,i0,i1,i2]=compute_dipscore(counts,weighted_counts)

nbins=length(counts);
updown_ks=[];
downup_ks=[];
i0s=[];
i1s=[];
i2s=[];
for j1=1:nbins
    for j2=j1+1:nbins
        wc=weighted_counts(j1:j2);
        updown=jisotonic(wc,'updown');
        downup=jisotonic(wc,'downup');
        updown_ks(end+1)=compute_ks(wc,updown)*sqrt(sum(counts(j1:j2)));
        downup_ks(end+1)=compute_ks(wc,downup)*sqrt(sum(counts(j1:j2)));
        [~,minind]=min(downup); minind=minind(1);
        i0s(end+1)=minind+j1-1;
        i1s(end+1)=j1;
        i2s(end+1)=j2;
    end;
end;
    
dipscores=updown_ks-downup_ks;
[~,ind]=max(dipscores);
dipscore=dipscores(ind);
i1=i1s(ind);
i2=i2s(ind);
i0=i0s(ind);

function Y=logsumexp(X)
val=max(X(:));
Y=val*log(sum(exp(X-val)));

function [ks,ks_ind]=compute_ks(counts1,counts2)
S1=cumsum(counts1);
S2=cumsum(counts2);
[ks,ks_ind]=max(abs(S1-S2)/sum(counts1));
ks_ind=ks_ind(1);

function test_isocut_A3
close all;

rng(3);

N0=200;

num_trials=1;

ks=[];
cutpoints=[];
for j=1:num_trials
    X1=randn(2,N0);
    X2=randn(2,N0); X2(1,:)=X2(1,:)+2;
    X3=randn(2,N0); X3(1,:)=X3(1,:)+7;
    X=cat(2,X1,X2,X3);
    
    opts=struct;
    opts.verbose=1;
    
    isocut_A3(X(1,:),[],opts);
    
    [inds_thin,weights_thin]=thin(X,200);
    X_thin=X(:,inds_thin);
    isocut_A3(X_thin(1,:),weights_thin,opts);
end;

