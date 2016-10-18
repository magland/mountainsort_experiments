function test_thinning_and_isosplit2b

close all;

rng(4);

N0=2000;
KNN=100;
desired_num=2000;
num_noise_dims=0;

%[X,labels]=create_multimodal_2d({[N0,0,0,1,1,0],[N0/2,4.5,0,1,3,0]});
%opts.K_init=5;

A1.N=N0*2; A1.center=[0,0]; A1.cov=[1,0;0,1];
A2.N=N0/2; A2.center=[5.5,0]; A2.cov=[1,0;0,3];
A3.N=N0/2; A3.center=[-5,5]; A3.cov=[3,1;1,3];
A4.N=N0/2; A4.center=[0,-8]; A4.cov=[0.5,0;0,0.5];
AA={A1,A2,A3,A4};
for j=1:length(AA)
    M=length(AA{j}.center);
    center2=rand(1,M+num_noise_dims)*10;
    center2(1:M)=AA{j}.center;
    AA{j}.center=center2;
    cov2=eye(M+num_noise_dims);
    cov2(1:M,1:M)=AA{j}.cov;
    AA{j}.cov=cov2;
end;
[X,labels]=create_multimodal_nd(AA);
N=size(X,2);

figure; ms_view_clusters(X(1:2,:),labels);

opts.return_iterations=1;
%[labels2,info]=isosplit2b(X,[],opts);
%figure; ms_view_clusters(X(1:2,:),labels2);
%show_iterations(X,info);

%[labels2,info]=isosplit2(X,opts);
%figure; ms_view_clusters(X(1:2,:),labels2);

[inds_thin,weights_thin]=thin(X,desired_num);
X_thin=X(:,inds_thin);
labels_thin=labels(inds_thin);
%figure; ms_view_clusters(X_thin(1:2,:),labels_thin);

opts.max_iterations=50;
[labels2,info]=isosplit2b(X_thin,weights_thin,opts);
figure; ms_view_clusters(X_thin(1:2,:),labels2);

show_iterations(X_thin,info);

function show_iterations(X,info)
for j=1:length(info.iterations)
    tmp=info.iterations{j};
    f=figure; set(f,'position',[100,100,2000,600]);
    subplot(1,2,1);
        ms_view_clusters(X,tmp.labels_before);
        set(gca,'xtick',[]); set(gca,'ytick',[]); xlabel(''); ylabel('');
        title(sprintf('Iteration %d',j), 'FontSize', 20);
    subplot(1,2,2);
        view_clusters_1d_w(tmp.projection,tmp.projection_weights,tmp.projection_labels);
        if (max(tmp.projection_labels)>1)
            vline0(tmp.projection_cutpoint,'k-');
        end;
        title(sprintf('Compare %d,%d',min(tmp.k1,tmp.k2),max(tmp.k1,tmp.k2)), 'FontSize', 20);
    %wait_for_key_press;
    %close(f);
end;

function view_clusters_1d(X,labels)
K=max(labels);
[~,bins]=hist(X,100);
for k=1:K
    inds=find(labels==k);
    if (length(inds)>0)
        vals=hist(X(inds),bins);
        if (k==1) col='k';
        else col='r';
        end;
        hh=bar(bins,vals); hold on;
        set(hh,'FaceColor',col,'EdgeColor',col,'BarWidth',1);
    end;
end;

set(gca,'xtick',[]); set(gca,'ytick',[]); xlabel(''); ylabel('');

function view_clusters_1d_w(X,W,labels)
K=max(labels);
bin_width=(max(X(:))-min(X(:)))/100;
bin_ints=round(X/bin_width);
i1=min(bin_ints);
i2=max(bin_ints);
bins=(i1:i2)*bin_width;
for k=1:K
    inds=find(labels==k);
    if (length(inds)>0)
        if (k==1) col='k';
        else col='r';
        end;
        ii=bin_ints(inds)-i1+1;
        vals=accumarray(ii',W(inds)',[length(bins),1]);
        hh=bar(bins,vals); hold on;
        set(hh,'FaceColor',col,'EdgeColor',col,'BarWidth',1);
    end;
end;

set(gca,'xtick',[]); set(gca,'ytick',[]); xlabel(''); ylabel('');

function vline0(x,linespec)
ylim0=get(gca,'ylim');
plot([x,x],[ylim0(1),ylim0(2)],linespec);