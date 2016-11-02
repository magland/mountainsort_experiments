function [dip_score,cutpoint,info]=isocut5(samples,sample_weights,opts)
if nargin<1, test_isocut5; return; end;
if nargin<2, sample_weights=[]; end;
if (nargin<3), opts=struct; end;

if (~isfield(opts,'already_sorted')) opts.already_sorted=0; end;
if (~isfield(opts,'try_ranges')) opts.try_ranges=0; end;
if (~isfield(opts,'range')) opts.range=[]; end;
if (~isfield(opts,'return_info')) opts.return_info=[]; end;

dip_score=0;
cutpoint=0;
info=struct;

[~,N]=size(samples);
if (N==0)
    error('Error in isocut5: N is zero.');
end;
num_bins_factor=1;
num_bins=ceil(sqrt(N/2)*num_bins_factor);

if (length(sample_weights)==0), sample_weights=ones(1,N); end;

if (opts.try_ranges)
    if (opts.already_sorted)
        X=samples;
    else
        [X,sort_inds]=sort(samples);
        sample_weights=sample_weights(sort_inds);
    end;
    opts2=opts;
    opts2.try_ranges=0;
    opts2.already_sorted=1;
    ranges_to_try={};
    N2=ceil(N/2);
    len=10;
    while ((len<N2)&&(len<N-N2))
        ranges_to_try{end+1}=[X(1),X(len)];
        ranges_to_try{end+1}=[X(N-len+1),X(N)];
        len=len*2;
    end;
    ranges_to_try{end+1}=[-inf,inf];
    best_dip_score=-inf;
    for j=1:length(ranges_to_try)
        opts2.range=ranges_to_try{j};
        if (nargout>=3)
            [ds,cp,info2]=isocut5(X,sample_weights,opts2);
        else
            [ds,cp]=isocut5(X,sample_weights,opts2);
            info2=struct;
        end;
        if (ds>best_dip_score)
            best_dip_score=ds;
            dip_score=ds;
            cutpoint=cp;
            info=info2;
        end;
    end;
    return;
end;

if (length(opts.range)==2)
    inds_in_range=find((opts.range(1)<=samples)&(samples<=opts.range(2)));
    samples_in_range=samples(inds_in_range);
    sample_weights_in_range=sample_weights(inds_in_range);
    opts2=opts;
    opts2.range=[];
    [dip_score,cutpoint,info]=isocut5(samples_in_range,sample_weights_in_range,opts2);
    if (nargout>=3)
        [counts,info.hist_bins]=hist_with_weights(samples,sample_weights,num_bins*3);
        info.hist_densities=counts/(info.hist_bins(2)-info.hist_bins(1));
    end;
    return;
end

if (opts.already_sorted)
    X=samples;
else
    [X,sort_inds]=sort(samples);
    sample_weights=sample_weights(sort_inds);
end;

num_bins_1=ceil(num_bins/2);
num_bins_2=num_bins-num_bins_1;
intervals=[1:num_bins_1,num_bins_2:-1:1];
alpha=(N-1)/sum(intervals);
intervals=intervals*alpha;
inds=floor([1,1+cumsum(intervals)]);
N_sub=length(inds);

cumsum_sample_weights=cumsum(sample_weights);

X_sub=X(inds);
spacings=X_sub(2:end)-X_sub(1:end-1);
%multiplicities=inds(2:end)-inds(1:end-1);
multiplicities=cumsum_sample_weights(inds(2:end))-cumsum_sample_weights(inds(1:end-1));
densities=multiplicities./spacings;

densities_unimodal_fit=jisotonic(densities,'updown',multiplicities);
densities_resid=densities-densities_unimodal_fit;
densities_resid_fit=jisotonic(densities_resid,'downup',multiplicities./densities);
[~,cutpoint_ind]=min(densities_resid_fit);
cutpoint=(X_sub(cutpoint_ind)+X_sub(cutpoint_ind+1))/2;

ks=compute_ks4(multiplicities,densities_unimodal_fit./densities.*multiplicities);
dip_score=ks;

if (opts.return_info)
    info.lefts=X_sub(1:end-1);
    info.rights=X_sub(2:end);
    info.centers=(info.lefts+info.rights)/2;
    info.densities=densities;
    info.densities_unimodal=densities_unimodal_fit;
    info.densities_bimodal=densities_resid_fit+densities_unimodal_fit;
    info.plot_xx=zeros(1,(N_sub-1)*2);
    info.plot_xx(1:2:end)=info.lefts;
    info.plot_xx(2:2:end)=info.rights;
    info.plot_densities=zeros(1,(N_sub-1)*2);
    info.plot_densities(1:2:end)=info.densities;
    info.plot_densities(2:2:end)=info.densities;
    info.plot_densities_unimodal=zeros(1,(N_sub-1)*2);
    info.plot_densities_unimodal(1:2:end)=info.densities_unimodal;
    info.plot_densities_unimodal(2:2:end)=info.densities_unimodal;
    info.plot_densities_bimodal=zeros(1,(N_sub-1)*2);
    info.plot_densities_bimodal(1:2:end)=info.densities_bimodal;
    info.plot_densities_bimodal(2:2:end)=info.densities_bimodal;
    [counts,info.hist_bins]=hist_with_weights(samples,sample_weights,num_bins*3);
    info.hist_densities=counts/(info.hist_bins(2)-info.hist_bins(1));
else
    info=struct;
end;

function [ks,ks_ind]=compute_ks4(counts1,counts2)
S1=cumsum(counts1)/sum(counts1);
S2=cumsum(counts2)/sum(counts2);
[ks,ks_ind]=max(abs(S1-S2));
ks=ks*sqrt((sum(counts1)+sum(counts2))/2);
ks_ind=ks_ind(1);

function [counts,bins]=hist_with_weights(X,weights,num_bins)
bin_width=(max(X(:))-min(X(:)))/num_bins;
bin_ints=round(X/bin_width);
i1=min(bin_ints);
i2=max(bin_ints);
bins=(i1:i2)*bin_width;
ii=bin_ints-i1+1;
counts=accumarray(ii',weights',[length(bins),1]);

function test_isocut5
close all;

opts.try_ranges=0;
num_trials=100;
cutpoints=[];
dip_scores=[];
run_times=[];

for trial=1:num_trials;

N0=1e4;
%X=[randn(1,N0)*1.2,randn(1,N0*2)*1.2+3.5,randn(1,N0/10)*1.2+12];
%X=[randn(1,2*N0)];

X=[randn(1,2*N0),randn(1,N0)+4];
sample_weights=[1*ones(1,2*N0),1*ones(1,N0)];

%X=[randn(1,N0/20)*0.4,randn(1,N0/10)+6,rand(1,N0*3*0)*12+6];
opts.return_info=(trial==1);

%sample_weights=[];
%sample_weights=abs(X);

tA=tic;
[dip_score,cutpoint,info]=isocut5(X,sample_weights,opts);
%cutpoint=isocut(X,1);
%dip_score=cutpoint;
run_time=toc(tA);

%tA=tic;
%[dip_score,cutpoint]=isocut4(X);
%run_time=toc(tA);

cutpoints(end+1)=cutpoint;
dip_scores(end+1)=dip_score;
run_times(end+1)=run_time;

if (trial==1)
    X_first=X;
end;

if (trial==1)    
    figure; hold on;
    bar(info.hist_bins,info.hist_densities,'FaceColor',[0.8,0.8,0.8],'EdgeColor',[0.8,0.8,0.8]);
    plot(info.plot_xx,info.plot_densities,'k','LineWidth',2);
    plot(info.plot_xx,info.plot_densities_unimodal,'g','LineWidth',2);
    plot(info.plot_xx,info.plot_densities_bimodal,'b','LineWidth',2);
    vline(cutpoint);
    title(sprintf('dip score = %g, N=%g',dip_score,length(X)));
    drawnow;
end;

end;

figure;
subplot(1,3,1);
hist(dip_scores,num_trials);
xlim([min(0,min(dip_scores)),max(dip_scores)]);
title(sprintf('dip scores - avg = %g',mean(dip_scores)));
subplot(1,3,2);
hist(cutpoints,num_trials);
xlim([min(X_first),max(X_first)]);
title(sprintf('cutpoints - avg = %g',mean(cutpoints)));
subplot(1,3,3);
hist(run_times,num_trials);
title(sprintf('run times - avg = %g seconds',round(mean(run_times)*100000)/100000));
set(gcf,'position',[50,50,1000,500]);
