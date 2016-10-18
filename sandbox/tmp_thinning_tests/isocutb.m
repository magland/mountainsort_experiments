function [ks_normalized,cutpoint]=isocutb(X,W,opts)

if (nargin<1) test_isocutb; return; end;

if (nargin<2) W=[]; end;
if (length(W)==0) W=ones(size(X)); end;

if (nargin<3) opts=struct; end;
if (~isfield(opts,'verbose')) opts.verbose=0; end;
if (~isfield(opts,'ks_threshold')) opts.ks_threshold=1.5; end;
if (~isfield(opts,'try_ends')) opts.try_ends=0; end;

[~,sort_inds]=sort(X);
X=X(sort_inds);
W=W(sort_inds);

if (opts.try_ends)
    N=length(X);
    N0=8;
    opts0=opts;
    opts0.try_ends=0;
    opts0.verbose=0;
    while (N0<N/2)
        for pass=1:2
            if (pass==1)
                X0=X(1:N0);
                W0=W(1:N0);
            else
                X0=X(end-N0+1:end);
                W0=W(end-N0+1:end);
            end;
            [ks0,cutpoint0]=isocutb(X0,W0,opts);
            if (ks0>=opts.ks_threshold)
                ks_normalized=ks0;
                cutpoint=cutpoint0;
                return;
            end
        end;
        N0=N0*2;
    end
end

[counts,bin_centers]=autobin(X,W);

counts_fit=jisotonic(counts,'updown');
[ks,ks_ind]=compute_ks(counts,counts_fit);
%ks_normalized=ks*sqrt(length(X)); % How to normalize in the case of weighted data?
%ks_normalized=ks*sqrt(sum(W.^(1/3)));
ks_normalized=ks*sqrt(sum(W));

[~,cut_index]=jisotonic_updownupdown(counts);
cutpoint=bin_centers(cut_index);

% counts_resid=counts-counts_fit;
% counts_resid_fit=jisotonic(counts_resid,'updown');
% [~,peak1]=max(counts_fit); peak1=peak1(1);
% [~,peak2]=max(counts_resid_fit); peak2=peak2(1);
% peak_lower=min(peak1,peak2);
% peak_upper=max(peak1,peak2);
% tmp=jisotonic(counts(peak_lower:peak_upper),'downup');
% [~,cut_index]=min(tmp); cut_index=cut_index(1);
% cut_index=peak_lower-1+cut_index;
% cutpoint=bin_centers(cut_index);

%counts_resid=counts-counts_fit;
%counts_resid_fit=jisotonic(counts_resid,'downup');
%[~,ind]=min(counts_resid_fit);
%ind=ind(1);
%cutpoint=bin_centers(ind);

if opts.verbose
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

bin_width=IQR/10;
%bin_width=IQR/N_IQR;
%bin_width=0.001;

minval=min(X(:));
maxval=max(X(:));
iminval=round(minval/bin_width);
imaxval=round(maxval/bin_width);
bin_centers=(iminval:imaxval)*bin_width;

inds=round(X/bin_width)-iminval+1;
counts=accumarray(inds',W,[length(bin_centers),1])';

%counts=hist(X,bin_centers);

function [ks,ks_ind]=compute_ks(counts1,counts2)
S1=cumsum(counts1);
S2=cumsum(counts2);
[ks,ks_ind]=max(abs(S1-S2)/sum(counts1));
ks_ind=ks_ind(1);


function test_isocutb
close all;

N0=300;

num_trials=10;

ks=[];
cutpoints=[];
for j=1:num_trials
    X1=randn(1,N0*7);
    X2=randn(1,N0/2)+5;
    X3=randn(1,N0*7)+10;
    X=cat(2,X1,X2,X3);
    if (j==1) opts.verbose=1;
    else opts.verbose=0;
    end;
    [ks0,cutpoint0]=isocutb(X,[],opts);
    ks(end+1)=ks0;
    cutpoints(end+1)=cutpoint0;
end;

figure; hist(ks,num_trials);
figure; hist(cutpoints,num_trials);
mean(ks)
sqrt(var(ks))
mean(cutpoints)
