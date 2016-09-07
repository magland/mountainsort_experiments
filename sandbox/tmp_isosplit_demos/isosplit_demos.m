function isosplit_demos

close all; 

%demo1;
demo2;

end

function demo2
N0=1000;

[X,labels]=create_multimodal_2d({[N0,0,0,1,1],[N0/2,4.5,0,1,3]});
opts.K_init=4;

%[X,labels]=create_multimodal_2d({[N0,0,0,1,1],[N0/2,4.5,0,1,3],[N0/2,8,8,1,1]});
%opts.K_init=4;

opts.return_iterations=1;

[labels2,info]=isosplit2(X,opts);
for j=1:length(info.iterations)
    tmp=info.iterations{j};
    f=figure; set(f,'position',[100,100,1000,600]);
    subplot(1,2,1);
        ms_view_clusters(X,tmp.labels_before);
    subplot(1,2,2);
        view_clusters_1d(tmp.projection,100);
        title(sprintf('Compare %d,%d',min(tmp.k1,tmp.k2),max(tmp.k1,tmp.k2)));
    wait_for_key_press;
    close(f);
end;
figure; ms_view_clusters(X,labels);
figure; ms_view_clusters(X,labels2);
end

function demo1

N0=5000;

X1=create_multimodal({[N0,0,1],[N0/2,2.5,1]});
N1=length(X1);
X1_fit=fit_unimodal(X1);

X2=create_multimodal({[N0,0,1],[N0/2,4,1]});
N2=length(X2);
X2_fit=fit_unimodal(X2);

figure; set(gcf,'position',[100,100,1200,600]);
subplot(2,3,1);
    hist2(X1,200);
    format_hist(gca,'Unimodal sample');
subplot(2,3,2);
    hist2(X1_fit,200);
    format_hist(gca,'Best unimodal fit');
subplot(2,3,3);
    plot(linspace(0,1,N1),X1,'b',linspace(0,1,N1),X1_fit,'r');
    legend('Sample','Best unimodal fit','Location','North');
    format_hist(gca,'Cumulative density function');
subplot(2,3,4);
    hist2(X2,200);
    format_hist(gca,'Bimodal sample');
subplot(2,3,5);
    hist2(X2_fit,200);
    format_hist(gca,'Best unimodal fit');
subplot(2,3,6);
    plot(linspace(0,1,N2),X2,'b',linspace(0,1,N2),X2_fit,'r');
    legend('Sample','Best unimodal fit','Location','North');
    format_hist(gca,'Cumulative density function');
end

function [X,labels]=create_multimodal_2d(A)
X=zeros(2,0);
labels=zeros(1,0);
for j=1:length(A)
    A0=A{j};
    tmp=randn(2,A0(1));
    b=[A0(4),0;0,A0(5)];
    tmp=b*tmp;
    tmp(1,:)=A0(2)+tmp(1,:);
    tmp(2,:)=A0(3)+tmp(2,:);
    X=cat(2,X,tmp);
    labels=cat(2,labels,j*ones(1,A0(1)));
end
end

function hist2(X,nbins)
[counts,bins]=hist(X,nbins);
bar(bins,counts,'k');
end

function format_hist(h,title0)
set(h,'xtick',[]);
set(h,'ytick',[]);
title(h,title0);
end

function X_fit=fit_unimodal(X)
spacing=diff(X);
spacing_fit=jisotonic(spacing,'downup');
X_fit=[X(1),X(1)+cumsum(spacing_fit)];
end

function X=create_multimodal(A)

X=zeros(1,0);
for j=1:length(A)
    A0=A{j};
    tmp=A0(3)*randn(1,A0(1))+A0(2);
    X=cat(2,X,tmp);
end;
X=sort(X);

end

function wait_for_key_press
while (waitforbuttonpress==0)
end;
end