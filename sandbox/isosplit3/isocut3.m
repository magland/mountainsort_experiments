function [dip_score,cutpoint,info]=isocut3(samples,weights,opts)

if nargin<1, test_isocut3; return; end;
if nargin<2, weights=[]; end;
if (length(weights)==0) weights=ones(size(samples)); end;
if nargin<3, opts=struct; end;

tic;

desired_num_crit_pts=10;
separation_interval=10;

num_bins=determine_optimal_num_bins(samples,weights,desired_num_crit_pts,separation_interval);

[counts,weighted_counts,bin_centers]=compute_hist(samples,weights,num_bins);
num_bins=length(counts);

critical_inds=find_local_maxima(weighted_counts,separation_interval);
critical_inds=sort([critical_inds,0,num_bins+1]);

dip_scores=[];
trials={};
for j1=1:length(critical_inds)
    for j2=j1+1:length(critical_inds)
        i1=critical_inds(j1)+1;
        i2=critical_inds(j2)-1;
        if (i2>i1)
            fit=jisotonic(weighted_counts(i1:i2),'downup');
            fit_control=jisotonic(weighted_counts(i1:i2),'updown');

            ks=compute_ks2(counts(i1:i2),fit,weighted_counts(i1:i2));
            ks_control=compute_ks2(counts(i1:i2),fit_control,weighted_counts(i1:i2));

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

[~,best_ind]=max(dip_scores); best_ind=best_ind(1);
best_trial=trials{best_ind};
dip_score=dip_scores(best_ind);
cutpoint=best_trial.cutpoint;

if (nargout>=3)
    info.dip_score=dip_score;
    info.cutpoint=cutpoint;
    info.dip_range=best_trial.dip_range;
    info.counts=counts;
    info.weighted_counts=weighted_counts;
    info.bin_centers=bin_centers;
    info.critical_inds=critical_inds;
    info.best_trial=best_trial;
    info.elapsed_time=toc;
end;

function num_bins=determine_optimal_num_bins(samples,weights,desired_num_crit_pts,separation_interval)
for num=10:10:length(samples)
    [counts,weighted_counts]=compute_hist(samples,weights,num);
    inds=find_local_maxima(weighted_counts,separation_interval);
    if (length(inds)>=desired_num_crit_pts)
        num_bins=num;
        return;
    end;
end
num_bins=length(samples);

function [counts,weighted_counts,bin_centers]=compute_hist(samples,weights,num_bins)
minval=min(samples);
maxval=max(samples);
bin_width=(maxval-minval)/(num_bins-2);
bin_indices=round((samples-minval)/bin_width)+1;
counts=accumarray(bin_indices',1,[num_bins,1])';
weighted_counts=accumarray(bin_indices',weights',[num_bins,1])';
bin_centers=minval+bin_width*(0:num_bins-1);

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

function [ks,ks_ind]=compute_ks2(counts_for_normalization,weighted_counts1,weighted_counts2)
normalization_count_forward=cumsum(counts_for_normalization);
normalization_count_reverse=cumsum(counts_for_normalization(end:-1:1));
normalization_count=min(normalization_count_forward,normalization_count_reverse(end:-1:1));
S1=cumsum(weighted_counts1)/sum(weighted_counts1);
S2=cumsum(weighted_counts2)/sum(weighted_counts2);
[ks,ks_ind]=max(abs(S1-S2).*sqrt(2*normalization_count));
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
        title(sprintf('dip score = %g, N = %d, #bins = %d',round2(dip_score),length(samples),length(info.counts)));
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
