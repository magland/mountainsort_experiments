function ret=very_approx_knn(X)

[M,N]=size(X);
if (N==1)
    ret=1;
    return;
end;
if (N==2)
    ret=[2,1];
    return;
end;
V=X(:,end)-X(:,1);
V=V/sqrt(V'*V);
proj=V'*X;
[~,sort_inds]=sort(proj);
if (N<=10)
    ret=sort_inds([2:N,N-1]);
    return;
end;

i1=ceil(1*N/4);
i2=ceil(2*N/4);
i3=ceil(3*N/4);

ret1=very_approx_knn(X(:,sort_inds(1:i1)));
ret2=very_approx_knn(X(:,sort_inds(i1+1:i2)));
ret3=very_approx_knn(X(:,sort_inds(i2+1:i3)));
ret4=very_approx_knn(X(:,sort_inds(i3+1:N)));

ret=zeros(1,N);
ret(sort_inds(1:i1))=sort_inds(ret1);
ret(sort_inds(i1+1:i2))=sort_inds(ret2);
ret(sort_inds(i2+1:i3))=sort_inds(ret3);
ret(sort_inds(i3+1:N))=sort_inds(ret4);
