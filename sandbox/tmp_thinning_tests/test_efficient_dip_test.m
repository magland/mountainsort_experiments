function test_efficient_dip_test

%test1;
test2

function test2
N0=800;
desired_num_crit_pts=10;
interval=10;
%samples=cat(2,randn(1,N0),randn(1,N0*2)+3,randn(1,N0)+6.2,randn(1,2000)*0.1+1);
samples=cat(2,randn(1,N0),randn(1,N0*2)+2,randn(1,N0*2)+5);

tic;
samples=sort(samples);
num_bins=determine_num_bins(samples,desired_num_crit_pts,interval)
[counts,bins]=hist(samples,num_bins);
N=length(counts);
X_piecewise_fit=zeros(1,N)*inf;
%X_piecewise_fit_control=zeros(1,N)*inf;

inds=find_local_maxima(counts,interval);
%inds=[inds,1,length(counts)];
inds=sort(inds);
length(inds)
%inds=[min(find(bins>0)),min(find(bins>4))];

ks_scores=[];
results={};
for j1=1:length(inds)
    for j2=j1+1:length(inds)
        i1=inds(j1)+1;
        i2=inds(j2)-1;
        fit=jisotonic(counts(i1:i2),'downup');
        %X_piecewise_fit(i1:i2)=fit;
        fit_control=jisotonic(counts(i1:i2),'updown');
        
        %X_piecewise_fit_control(i1:i2)=fit;
        
        ks=compute_ks(fit,counts(i1:i2))*sqrt(sum(counts(i1:i2)));
        ks_control=compute_ks(fit_control,counts(i1:i2))*sqrt(sum(counts(i1:i2)));
        
        result.range=i1:i2;
        result.ks=ks;
        result.ks_control=ks_control;
        result.fit=fit;
        result.fit_control=fit_control;
        ks_scores(end+1)=ks_control-ks;
        results{end+1}=result;
    end;
end;

[~,best_ind]=max(ks_scores); best_ind=best_ind(1);
result=results{best_ind};
N
toc

figA=figure;
bar(bins,counts,'g');
hold on;
plot(bins(result.range),result.fit_control,'b','LineWidth',1);
plot(bins(result.range),result.fit,'r','LineWidth',2);
title(sprintf('ks score = %g',ks_scores(best_ind)));

function num_bins=determine_num_bins(samples,desired_num_crit_pts,interval)
for num=10:10:length(samples)
    counts=hist(samples,num);
    inds=find_local_maxima(counts,interval);
    if (length(inds)>=desired_num_crit_pts)
        num_bins=num;
        return;
    end;
end
num_bins=length(samples);


function [ks,ks_ind]=compute_ks(counts1,counts2)
S1=cumsum(counts1);
S2=cumsum(counts2);
[ks,ks_ind]=max(abs(S1-S2)/sum(counts1));
ks_ind=ks_ind(1);



function inds=find_local_maxima(X,interval,threshold)
if (nargin<3) threshold=-inf; end;
N=length(X);
use_it=zeros(1,N);
best_ind=1;
best_val=X(1);
candidates=find(X>=threshold);
for tt=candidates
    if (best_ind<tt-interval)
        [~,best_ind]=max(X(tt-interval:tt-1));
        best_ind=best_ind+tt-interval-1;
        best_val=X(best_ind);
    end;
    if (X(tt)>=best_val)
        use_it(tt)=1;
        use_it(best_ind)=0;
        best_ind=tt;
        best_val=X(tt);
    end;
end;

inds=find(use_it==1);


function test1

N=100;
mid=50;
X=rand(1,N);
X1_fit=jisotonic(X(1:mid),'downup');
X2_fit=jisotonic(X(mid+1:end),'downup');
X0_fit=cat(2,X1_fit,X2_fit);
X_fit=jisotonic(X,'downup');

figure; plot(1:N,X,'k');
hold on;
plot(1:N,X0_fit,'b','LineWidth',4);
plot(1:N,X_fit,'g');
