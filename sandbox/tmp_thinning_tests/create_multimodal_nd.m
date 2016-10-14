function [X,labels]=create_multimodal_nd(A)
M=length(A{1}.center);
X=zeros(M,0);
labels=zeros(1,0);
for j=1:length(A)
    A0=A{j};
    tmp=randn(M,A0.N);
    %b=[A0(4),A0(6);-A0(6),A0(5)];
    tmp=A0.cov*tmp;
    for m=1:M
        tmp(m,:)=A0.center(m)+tmp(m,:);
    end;
    X=cat(2,X,tmp);
    labels=cat(2,labels,j*ones(1,A0.N));
end