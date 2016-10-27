function canhist(samples,tmp)
if nargin<1, test_canhist; return; end;

tic;
X=sort(samples);
N=length(X);

num_bins=ceil(sqrt(N/2)*tmp);
num_bins_1=ceil(num_bins/2);
num_bins_2=num_bins-num_bins_1;
intervals=[1:num_bins_1,num_bins_2:-1:1];
alpha=(N-1)/sum(intervals);
intervals=intervals*alpha;
inds=floor([1,1+cumsum(intervals)]);

N_sub=length(inds);

X_sub=X(inds);
spacings=X_sub(2:end)-X_sub(1:end-1);
multiplicities=inds(2:end)-inds(1:end-1);
lefts=X_sub(1:end-1);
rights=X_sub(2:end);
densities=multiplicities./spacings;

densities_fit=jisotonic(densities,'updown',multiplicities);
densities_resid=densities-densities_fit;
densities_resid_fit=jisotonic(densities_resid,'downup',multiplicities);
densities_fit2=densities_resid_fit+densities_fit;

ks=compute_ks2(multiplicities,densities_fit./densities.*multiplicities);
ks3=compute_ks3(multiplicities,densities_fit./densities.*multiplicities);

xx=zeros(1,(N_sub-1)*2);
xx(1:2:end)=lefts;
xx(2:2:end)=rights;
yy=zeros(1,(N_sub-1)*2);
yy(1:2:end)=densities;
yy(2:2:end)=densities;
yy_fit=zeros(1,(N_sub-1)*2);
yy_fit(1:2:end)=densities_fit;
yy_fit(2:2:end)=densities_fit;
yy_fit2=zeros(1,(N_sub-1)*2);
yy_fit2(1:2:end)=densities_fit2;
yy_fit2(2:2:end)=densities_fit2;

[counts,bins]=hist(X,num_bins*3);
bin_width=bins(2)-bins(1);
figure; bar(bins,counts/bin_width,'EdgeColor',[0.8,0.8,0.8],'FaceColor',[0.8,0.8,0.8]);
hold on;
plot(xx,yy,'k','LineWidth',3);
plot(xx,yy_fit,'g','LineWidth',2);
plot(xx,yy_fit2,'b','LineWidth',2);

title(sprintf('ks = %g, ks3 = %g',ks,ks3));

% best_j=1;
% best_k=1;
% best_dip_score=-inf;
% best_densities_downup=[];
% best_densities_updown=[];
% for j=1:4:length(densities)
%     for k=j+1:4:length(densities)
%         densities0=densities(j:k);
%         multiplicities0=multiplicities(j:k);
%         densities_downup=jisotonic(densities0,'downup',multiplicities0);
%         densities_updown=jisotonic(densities0,'updown',multiplicities0);
%         ks_downup=compute_ks2(densities_downup./densities0.*multiplicities0,multiplicities0);
%         ks_updown=compute_ks2(densities_updown./densities0.*multiplicities0,multiplicities0);
%         %dip_score=(ks_updown-ks_downup)*sqrt(sum(multiplicities0));
%         dip_score=(ks_updown-ks_downup);
%         if (dip_score>best_dip_score)
%             best_j=j;
%             best_k=k;
%             best_dip_score=dip_score;
%             best_densities_downup=densities_downup;
%             best_densities_updown=densities_updown;
%         end;
%     end;
% end;



% xx0=zeros(1,2*length(best_j:best_k));
% yy0=zeros(1,2*length(best_j:best_k));
% xx0(1:2:end)=lefts(best_j:best_k);
% xx0(2:2:end)=rights(best_j:best_k);
% yy0(1:2:end)=best_densities_downup;
% yy0(2:2:end)=best_densities_downup;
% plot(xx0,yy0,'r','LineWidth',3);
% yy0(1:2:end)=best_densities_updown;
% yy0(2:2:end)=best_densities_updown;
% plot(xx0,yy0,'g','LineWidth',3);
% 
% title(sprintf('best dip score = %g' ,best_dip_score));

function [ks,ks_ind]=compute_ks2(weighted_counts1,weighted_counts2)
S1=cumsum(weighted_counts1)/sum(weighted_counts1);
S2=cumsum(weighted_counts2)/sum(weighted_counts2);
%[ks,ks_ind]=max(abs(S1-S2).*sqrt(2*normalization_count));
[ks,ks_ind]=max(abs(S1-S2));
ks=ks*sqrt((sum(weighted_counts1)+sum(weighted_counts2))/2);
ks_ind=ks_ind(1);

function [ks,ks_ind]=compute_ks3(weighted_counts1,weighted_counts2)
S1=cumsum(weighted_counts1)/sum(weighted_counts1);
S2=cumsum(weighted_counts2)/sum(weighted_counts2);

N=length(weighted_counts1);
N2=ceil(N/2);

W=(weighted_counts1+weighted_counts2)/2;
normalization_counts=zeros(1,N);
normalization_counts(1:N2)=cumsum(W(1:N2));
normalization_counts(N:-1:N2+1)=cumsum(W(N:-1:N2+1));

normalization=sqrt(sum(W))*sqrt( max(normalization_counts)./normalization_counts );

[ks,ks_ind]=max(abs(S1-S2).*normalization);
ks_ind=ks_ind(1);

function test_canhist

close all;

N0=4000;
%X=[randn(1,N0)*1.2,randn(1,N0*2)*1.2+3.5,randn(1,N0)*1.2+9];
%X=[randn(1,2*N0)];
X=[randn(1,N0/20)*0.4,randn(1,N0/10)+6,rand(1,N0*3)*12+6];
canhist(X,0.5);
canhist(X,1);
canhist(X,2);



