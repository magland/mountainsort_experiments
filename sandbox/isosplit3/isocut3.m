function [dip_score,cutpoint,info]=isocut3(samples,weights,diameters,opts)

if nargin<1, test_isocut3; return; end;
if nargin<2, weights=[]; end;
if (length(weights)==0) weights=ones(size(samples)); end;
if nargin<3, diameters=[]; end;
if (length(diameters)==0) diameters=zeros(size(samples)); end;
if nargin<4, opts=struct; end;

tic;

%desired_num_crit_pts=10;
%separation_interval=10;
%num_bins=determine_optimal_num_bins(samples,weights,diameters,desired_num_crit_pts,separation_interval);
%[weighted_counts,bin_centers]=compute_hist(samples,weights,diameters,num_bins);
%num_bins=length(bin_centers);
%critical_inds=find_local_maxima(weighted_counts,separation_interval);
%critical_inds=sort(critical_inds);
%if (length(critical_inds)==0)
%    error('No critical indices found!');
%end;

%bin_width=compute_optimal_bin_width(samples,weights,diameters);
if (isfield(opts,'bin_width'))
    bin_width=opts.bin_width;
else
    minval=min(samples-diameters/2);
    maxval=max(samples+diameters/2);
    if (maxval==minval)
        bin_width=1;
    else
        bin_width=(maxval-minval)/100;
    end
end;
[weighted_counts,bin_centers]=compute_hist(samples,weights,diameters,bin_width);
num_bins=length(bin_centers);
fprintf('isocut: Using %d bins, bin_width=%g\n',num_bins,bin_width);
interval=4;
critical_inds=get_critical_inds(weighted_counts,interval);

dip_scores=[];
trials={};
for j1=1:length(critical_inds)
    for j2=j1+1:length(critical_inds)
        i1=critical_inds(j1);
        i2=critical_inds(j2);
        if (i2>i1)
            fit=jisotonic(weighted_counts(i1:i2),'downup');
            fit_control=jisotonic(weighted_counts(i1:i2),'updown');

            ks=compute_ks2(fit,weighted_counts(i1:i2));
            ks_control=compute_ks2(fit_control,weighted_counts(i1:i2));

            trial.indices=i1:i2;
            trial.ks=ks;
            trial.ks_control=ks_control;
            trial.fit=fit;
            trial.fit_control=fit_control;
            trial.dip_range=[bin_centers(i1),bin_centers(i2)];
            [~,minind]=min(fit); minind=minind(1);
            minind=minind+i1-1;
            trial.cut_index=minind;
            trial.cutpoint=bin_centers(minind);
            dip_scores(end+1)=ks_control-ks;
            trials{end+1}=trial;
        end;
    end;
end;

if (~isempty(dip_scores))
    [~,best_ind]=max(dip_scores); best_ind=best_ind(1);
    best_trial=trials{best_ind};
    dip_score=dip_scores(best_ind);
    cutpoint=best_trial.cutpoint;
else
    best_ind=0;
    best_trial.indices=[];
    best_trial.ks=0;
    best_trial.ks_control=0;
    best_trial.fit=[];
    best_trial.fit_control=[];
    best_trial.dip_range=[0,0];
    dip_score=0;
    cutpoint=0;
end;

if (nargout>=3)
    info.dip_score=dip_score;
    info.cutpoint=cutpoint;
    info.dip_range=best_trial.dip_range;
    info.weighted_counts=weighted_counts;
    info.bin_centers=bin_centers;
    info.critical_inds=critical_inds;
    info.best_trial=best_trial;
    info.elapsed_time=toc;
end;

function ret=get_critical_inds(weighted_counts,interval)
ret=[1];
j1=2;
while (j1<length(weighted_counts))
    j2=min(j1+interval-1,length(weighted_counts)-1);
    [~,ii]=max(weighted_counts(j1:j2)); ii=ii(1);
    ret(end+1)=j1+ii-1;
    j1=j1+interval;
end;
ret=[ret,length(weighted_counts)];

function num_bins=determine_optimal_num_bins(samples,weights,diameters,desired_num_crit_pts,separation_interval)
for num=10:10:length(samples)
    [weighted_counts]=compute_hist(samples,weights,diameters,num);
    inds=find_local_maxima(weighted_counts,separation_interval);
    if (length(inds)>=desired_num_crit_pts)
        num_bins=num;
        return;
    end;
end
num_bins=length(samples);

function bin_width=compute_optimal_bin_width(samples,weights,diameters)
% Not sure about this function!
minval=min(samples-diameters/2);
maxval=max(samples+diameters/2);
if (maxval==minval)
    bin_width=1;
    return;
end
bin_width=(maxval-minval)/100;
return;
golden_bin_width=(maxval-minval)/1000;
[golden_hist,golden_bin_centers]=compute_hist(samples,weights,diameters,golden_bin_width);
factor=1;
done=0;
ks_threshold=0.5; %Not sure if this is a good criteria!!
while 1
    trial_factor=ceil(factor*1.5);
    trial_hist=downsample_hist(golden_hist,trial_factor);
    ks0=compute_ks2(golden_hist,trial_hist);
    if ((ks0>=ks_threshold)||(golden_bin_width*factor>(maxval-minval)/3))
        bin_width=golden_bin_width*factor;
        factor
        return;
    end;
    factor=trial_factor;
end;

function Y=downsample_hist(X,factor)
N=length(X);
N2=ceil(N/factor)*factor;
X2=zeros(1,N2);
X2(1:N)=X;
Y2=reshape(repmat(mean(reshape(X2,factor,N2/factor),1),factor,1),1,N2);
Y=Y2(1:N); %Does not handle the last few bins quite right, if N is not divisible by factor


function [weighted_counts,bin_centers]=compute_hist(samples,weights,diameters,bin_width)

bin_ints_lower=round((samples-diameters/2)/bin_width);
bin_ints_upper=round((samples+diameters/2)/bin_width);
bin_spans=bin_ints_upper-bin_ints_lower+1;

min_bin_int=min(bin_ints_lower);
max_bin_int=max(bin_ints_upper);
bin_centers=(min_bin_int:max_bin_int)*bin_width; % question: is it bad to snap to a grid here?
num_bins=length(bin_centers);

bin_inds_lower=bin_ints_lower-min_bin_int+1;
bin_inds_upper=bin_ints_upper-min_bin_int+1;

weighted_counts=zeros(1,num_bins);
for offset=0:(max(bin_spans)-1)
    to_use=find(bin_inds_lower+offset<=bin_inds_upper);
    weighted_counts=weighted_counts+accumarray(bin_inds_lower(to_use)'+offset,weights(to_use)'./bin_spans(to_use)',[num_bins,1])';
end;

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
        use_it(best_ind)=0; %important to do this first in case best_ind=tt=1 (fixed 10/20/16)
        use_it(tt)=1;
        best_ind=tt;
        best_val=X(tt);
    end;
end;

inds=find(use_it==1);

function [ks,ks_ind]=compute_ks2(weighted_counts1,weighted_counts2)
S1=cumsum(weighted_counts1)/sum(weighted_counts1);
S2=cumsum(weighted_counts2)/sum(weighted_counts2);
%[ks,ks_ind]=max(abs(S1-S2).*sqrt(2*normalization_count));
[ks,ks_ind]=max(abs(S1-S2));
ks=ks*sqrt((sum(weighted_counts1)+sum(weighted_counts2))/2);
ks_ind=ks_ind(1);

function y=round2(x)
y=round(x,2,'significant');

function test_isocut3

close all;

N0=100;
num_trials=100;

dip_scores=[];
cutpoints=[];
timers=[];
info1=struct;
for trial=1:num_trials
    %samples=cat(2,randn(1,N0*1.5)*2,randn(1,N0*2)+3,randn(1,N0*2)*2+6);
    samples=cat(2,randn(1,N0*1.5)-1,randn(1,N0*2)+3,randn(1,N0*2)+6);
    %samples=cat(2,randn(1,N0*1.5),randn(1,10)*0.1+0);
    %samples=cat(2,randn(1,300),randn(1,300)+4);
    [dip_score,cutpoint,info]=isocut3(samples,[]);
    result=info.best_trial;
    dip_scores(end+1)=dip_score;
    cutpoints(end+1)=cutpoint;
    timers(end+1)=info.elapsed_time;

    if trial==1
        info1=info;
        figA=figure;
        bar(info.bin_centers,info.weighted_counts,'g');
        hold on;
        plot(info.bin_centers(result.indices),result.fit_control,'b','LineWidth',1);
        plot(info.bin_centers(result.indices),result.fit,'r','LineWidth',2);
        title(sprintf('dip score = %g, N = %d, #bins = %d',round2(dip_score),length(samples),length(info.weighted_counts)));
        drawnow;
    end;
end;

figB=figure;
set(figB,'position',[100,100,1400,600]);
subplot(1,2,1);
hist(dip_scores,length(dip_scores)/3);
xlim0=xlim; xlim([0,xlim0(2)]);
title(sprintf('dip scores: %g +/- %g',round2(mean(dip_scores)),round2(sqrt(var(dip_scores)))));
subplot(1,2,2);
bar(info1.bin_centers,info1.weighted_counts,'FaceColor','g','EdgeColor','g');
hold on;
for j=1:length(cutpoints)
    vline(cutpoints(j),'k');
end
bar(info1.bin_centers,info1.weighted_counts,'FaceColor','g','EdgeColor','g');
title(sprintf('cutpoints: %g +/- %g',round2(mean(cutpoints)),round2(sqrt(var(cutpoints)))));

fprintf('Average time per run: %g sec\n',mean(timers));
