function [Y,cut_index]=jisotonic_updownupdown(X,weights)
if (nargin<2) weights=[]; end;
if (length(weights)==0) weights=ones(size(X)); end;
N=length(X);
best_index=1;
best_SSE=inf;
Y=zeros(size(X));
figure;
for j=1:5:length(X)
    X1=X(1:j); W1=weights(1:j);
    X2=X(j:end); W2=weights(j:end);
    A1=jisotonic(X1,'updown',W1);
    A2=jisotonic(X2,'updown',W2);
    Y0=zeros(size(X));
    Y0(1:j)=A1;
    Y0(j:end)=A2;
    SSE=sum(weights.*(Y0-X).^2);
    if (SSE<=best_SSE)
        best_SSE=SSE;
        best_index=j;
        Y=Y0;
    end;
    plot(1:length(X),X,'k',1:length(Y0),Y0,'b',1:length(Y0),Y,'r');
    title(sprintf('%d/%d',j,length(X)));
    pause(0.1);
end
cut_index=best_index;
