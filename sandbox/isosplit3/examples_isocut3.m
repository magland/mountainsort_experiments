function examples_isocut3

example1

function example1

N0=20;

opts.show_updown=0;

fA=figure;
total_isocut_time=0;
total_compare_time=0;
for j=1:16
    pops=[N0,N0,N0];
    centers=[-3,0,1+0.5*j];
    spreads=[1,1,1];
    samples=create_sample(pops,centers,spreads);

    tic;
    %[lefts,rights,counts]=split_into_bins(samples);
    %split_into_bins(samples);
    total_compare_time=total_compare_time+toc;
    tic;
    [dip_score,cutpoint,info]=isocut3(samples);
    total_isocut_time=total_isocut_time+toc;

    figure(fA);
    subplot(4,4,j);
    plot_result(samples,info,opts);
end;

total_compare_time
total_isocut_time

% function split_into_bins(samples)
% X=sort(samples);
% tmp=sum(X.*(1:length(X)));
% figure; plot(X,(1:length(X)).*tmp);
% max(X-(1:length(X)).*tmp);


function plot_result(samples,info,opts)
dip_score=info.dip_score;
cutpoint=info.cutpoint;
result=info.best_trial;
bar(info.bin_centers,info.weighted_counts,'g');
hold on;
if (opts.show_updown)
    plot(info.bin_centers(result.indices),result.fit_control,'b','LineWidth',1);
end;
plot(info.bin_centers(result.indices),result.fit,'r','LineWidth',2);
vline(info.cutpoint);
title(sprintf('dip score = %g, N = %d, #bins = %d',round2(dip_score),length(samples),length(info.weighted_counts)));
drawnow;

function y=round2(x)
y=round(x,2,'significant');

function samples=create_sample(pops,centers,spreads)
samples=zeros(1,0);
for j=1:length(pops)
    X0=randn(1,pops(j))*spreads(j)+centers(j);
    samples=cat(2,samples,X0);
end;