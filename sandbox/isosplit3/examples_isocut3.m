function examples_isocut3

example1

function example1

N0=100;

opts.show_updown=0;

fA=figure;
for j=1:9
    pops=[N0,N0,N0];
    centers=[-3,0,1+0.5*j];
    spreads=[1,1,1];
    samples=create_sample(pops,centers,spreads);

    [dip_score,cutpoint,info]=isocut3(samples);

    figure(fA);
    subplot(3,3,j);
    plot_result(samples,info,opts);
end;

function plot_result(samples,info,opts)
dip_score=info.dip_score;
cutpoint=info.cutpoint;
result=info.best_trial;
bar(info.bin_centers,info.counts,'g');
hold on;
if (opts.show_updown)
    plot(info.bin_centers(result.indices),result.fit_control,'b','LineWidth',1);
end;
plot(info.bin_centers(result.indices),result.fit,'r','LineWidth',2);
vline(info.cutpoint);
title(sprintf('dip score = %g, N = %d, #bins = %d',round2(dip_score),length(samples),length(info.counts)));
drawnow;

function y=round2(x)
y=round(x,2,'significant');

function samples=create_sample(pops,centers,spreads)
samples=zeros(1,0);
for j=1:length(pops)
    X0=randn(1,pops(j))*spreads(j)+centers(j);
    samples=cat(2,samples,X0);
end;