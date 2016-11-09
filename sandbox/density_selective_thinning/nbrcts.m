function ret=nbrcts(X,Y,eps,opts)
if nargin<1 test_nbrcts; return; end;
if nargin<4 opts=struct; end;
if ~isfield(opts,'recursion_level') opts.recursion_level=1; end;

global num_calls;
num_calls=num_calls+1;

[M,NX]=size(X);
[M,NY]=size(Y);
disp([NX,NY]);

if ((NX<=200)||(NY<=200))
    ret=exhaustive_count_neighbors(X,Y,eps);
    return;
end;


if (M==0)
    ret=zeros(1,NY);
    return;
end;

ret=zeros(1,NY);

X1=X(1,:);
Y1=Y(1,:);

xmin=min(X1);
X_bin_ints=floor((X1-xmin)/eps);
Y_bin_ints=floor((Y1-xmin)/eps);

i1=min(X_bin_ints);
i2=max(X_bin_ints);

%fprintf('rl=%d, num_bins=%d\n',opts.recursion_level,i2-i1+1);

opts2=opts;
opts2.recursion_level=opts.recursion_level+1;

for ii=i1:i2
    if (M==10) disp([ii,i1,i2]); end;
    X_inds=find(X_bin_ints==ii);
    if (length(X_inds)>0)
        if (M==1)
            Y_inds=find(Y_bin_ints==ii);
            if (length(Y_inds)>0)
                ret(Y_inds)=ret(Y_inds)+length(X_inds);
            end;
            Y_inds_left=find(Y_bin_ints==ii-1);
            if (length(Y_inds_left)>0)
                ret(Y_inds_left)=ret(Y_inds_left)+nbrcts2(X(1,X_inds),Y(1,Y_inds_left),eps,opts);
            end;
            Y_inds_right=find(Y_bin_ints==ii+1);
            if (length(Y_inds_right)>0)
                ret(Y_inds_right)=ret(Y_inds_right)+nbrcts2(X(1,X_inds),Y(1,Y_inds_right),eps,opts);
            end;
        else
            Y_inds=find((ii-1<=Y_bin_ints)&(Y_bin_ints<=ii+1));
            if (length(Y_inds)>0)
                ret(Y_inds)=ret(Y_inds)+nbrcts(X(2:end,X_inds),Y(2:end,Y_inds),eps,opts2);
            end;
        end;
    end;
end;

function ret=nbrcts2(x,y,eps,opts)
Nx=length(x);
Ny=length(y);

%ret=ones(1,Ny)*Nx;
%return;

xs=sort(x);
ret=zeros(1,Ny);

aa=20;
for j=1:aa:Nx
    j1=j;
    j2=min(j+aa-1,Nx);
    x1=xs(j1);
    x2=xs(j2);
    y_inds=find((x1-eps<=y)&(y<=x2+eps));
    if (length(y_inds)>0)
        ret(y_inds)=ret(y_inds)+nbrcts3(xs(j1:j2),y(y_inds),eps,opts);
    end;
end;
    
function ret=nbrcts3(x,y,eps,opts)
[xx,yy]=ndgrid(x,y);
diffs=abs(xx-yy);
ret=sum(diffs<=eps,1);

function ret=exhaustive_count_neighbors(X,Y,eps)
[M,NX]=size(X);
[M,NY]=size(Y);
maxdists=zeros(NX,NY);
for m=1:M
    maxdists=max(maxdists,abs(repmat(X(m,:),NY,1)'-repmat(Y(m,:),NX,1)));
end;
tmp=(maxdists<=eps);
ret=sum(tmp,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_nbrcts

global num_calls;
num_calls=0;

close all;

M=10;
N=60000;

X=randn(M,N);

eps=sqrt(M)*0.6;

tic;
counts=nbrcts(X,X,eps);
toc

figure; hist(counts,100);
figure; ms_view_clusters(X(1:2,:));
figure; ms_view_clusters(X(1:2,find(counts>200)));

disp(num_calls);
