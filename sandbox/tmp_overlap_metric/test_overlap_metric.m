function test_overlap_metric

%close all;

N=2000;
K=1;
eta=2;

scores1=[];
scores2=[];

%for cutoff=0:0.5:1
for eta=0.1:0.025:5
for cutoff=0
    for M=20
        X1=randn(M,N/2);
        X2=randn(M,N/2)*0.5;
        X1(1,:)=X1(1,:)-eta;
        X2(1,:)=X2(1,:)+eta;
        X=cat(2,X1,X2);
        labels=zeros(1,N);
        labels(find(X(1,:)<cutoff))=1;
        labels(find(X(1,:)>=cutoff))=2;
        labels_true=cat(2,ones(1,N/2),2*ones(1,N/2));
        %figure; ms_view_clusters(X,labels);

        knn_inds=knnsearch(X',X','K',K+1)';
        knn_inds=knn_inds(2:end,:);
        labels2=zeros(K,N);
        for k=1:K
            labels2(k,find(X(1,knn_inds(k,:))<cutoff))=1;
            labels2(k,find(X(1,knn_inds(k,:))>=cutoff))=2;
        end;

        labels3=mean(labels2,1);
        labels3(find(labels3<1.5))=1;
        labels3(find(labels3>=1.5))=2;
        if (M==3) figure; ms_view_clusters(X,labels3); end;

        frac=length(find(labels3~=labels))/N;
        isolation_metric(M)=1-frac;
    end;

    %figure; plot(1:M,isolation_metric);
    %xlabel('Number of dimensions');
    %title('Isolation metric');
    %ylim([0,1]);
    %if (cutoff==0)
    %    title({'Isolation metrics for evenly split Gaussians','in various dimensions',sprintf('cutoff = 0 sigma')});
    %else
    %    title({'Isolation metrics for unevenly split Gaussians','in various dimensions',sprintf('cutoff = %g sigma',cutoff)});
    %end;
    
    scores1(end+1)=length(find(labels==labels_true))/N;
    scores2(end+1)=isolation_metric(end);
    
end;
end;

figure; plot(scores1,scores2,'b.');
