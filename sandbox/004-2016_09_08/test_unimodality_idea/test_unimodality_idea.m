function test_unimodality_idea

M=4;
N=3000;
X1=randn(M,N);
X2=randn(M,N); X2(1,:)=X2(1,:)+4;
X=cat(2,X1,X2);
tmp=get_diff_hist(X,15000);
figure; hist(tmp,200);

end

function ret=get_diff_hist(X,num)
inds1=randsample(length(X),num,true);
inds2=randsample(length(X),num,true);
ret=X(:,inds1)-X(:,inds2);
ret=sqrt(sum(ret.^2,1));
end