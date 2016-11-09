function knn_timing_tests

close all;

Ms=[2,5,10,20];
Ns=[1e3,2e3,4e3,8e3,16e3,32e3];
Ks=[1,10,20];

results=[];

for iK=1:length(Ks)
    disp(iK);
    for iM=1:length(Ms)
        for iN=1:length(Ns)
            results(iM,iN,iK)=do_time(Ms(iM),Ns(iN),Ks(iK));
        end;
    end;
end;

for iN=1:length(Ns)
    figure;
    hold on;
    lll={};
    for iM=1:length(Ms)
        plot(Ks,squeeze(results(iM,iN,:)));
        xlabel('K');
        ylabel('time');
        lll{end+1}=sprintf('M = %d',Ms(iM));
    end;
    title(sprintf('N = %d',Ns(iN)));
    legend(lll);
end;

for iK=1:length(Ks)
    figure;
    hold on;
    lll={};
    for iM=1:length(Ms)
        plot(Ns,squeeze(results(iM,:,iK)));
        xlabel('N');
        ylabel('time');
        lll{end+1}=sprintf('M = %d',Ms(iM));
    end;
    title(sprintf('K = %d',Ks(iK)));
    legend(lll);
end;

disp(results);

function ret=do_time(M,N,K)
X=rand(M,N);
ttt=tic;
nearest=knnsearch(X',X','K',K+1);
ret=toc(ttt);
