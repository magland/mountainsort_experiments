% X - M x NX
% Y - M x NY
% eps - a small number
% ret - 1 x NY
% ret(j) is the approx number of points in X that are within eps euclidean
% distance of Y(:,j)

function ret=count_neighbors(X,Y,eps,opts)
if nargin<1, test_count_neighbors; return; end;
if nargin<4, opts=struct; end;
if ~isfield(opts,'max_num_trials') opts.max_num_trials=5; end;
if ~isfield(opts,'recursion_level') opts.recursion_level=1; end;

[M,NX]=size(X);
[M,NY]=size(Y);

if (opts.recursion_level>5)
    ret=ones(1,NY)*NX;
    return;
end;

if ((NX<=10)||(NY<=10))
    ret=exhaustive_count_neighbors(X,Y,eps);
    return;
end;

max_num_trials=opts.max_num_trials;

opts2=opts;
opts2.recursion_level=opts.recursion_level+1;

trial_number=0;
while 1
    % project onto a random direction
    V=randn(M,1);
    V=V/sqrt(V'*V);
    X_projection=V'*X;
    Y_projection=V'*Y;
    
    % look to see which points are even candidates for counting
    Y_candidate_inds=find((Y_projection>=min(X_projection)-eps)&(Y_projection<=max(X_projection)+eps));
    if (length(Y_candidate_inds)==0)
        % There are no candidates
        ret=zeros(1,NY);
        return;
    end;
    
    % check to see if we have no eliminations all together
    if ((abs(min(Y_projection)-max(X_projection))<=eps)&&(abs(max(Y_projection)-min(X_projection))<=eps))
        % No eliminations here, we can try again up to max_num_trials times
        trial_number=trial_number+1;
        if (trial_number>=max_num_trials)
            % We failed to find a direction for which any of the Y's
            % are further than eps away from one any of the X's
            % so, let's count them all!
            ret=ones(1,NY)*NX;
            if (NX>1) fprintf('%d, rl=%d\n',NX,opts.recursion_level); end;
            return;
        end;
    else
        %Split into bins of size eps
        X_minval=min(X_projection);
        X_bin_indices=floor((X_projection-X_minval)/eps);
        Y_bin_indices=floor((Y_projection-X_minval)/eps);
        ii1=min(X_bin_indices);
        ii2=max(X_bin_indices);
        num_bins=ii2-ii1+1;
        %bin_lefts=X_minval+(ii1:ii2)*eps;

        if (num_bins==1)
            if (NX==1)
                % There is only one X-point, and based on the above checks,
                % this situation should be impossible.
                % Because being a candidate is the same thing as being in range
                % of every X point
                error('Unexpected problem');
            end;
            jj1=1:ceil(NX/2);
            jj2=ceil(NX/2)+1:NX;
            ret1=count_neighbors(X(:,jj1),Y,eps,opts2);
            ret2=count_neighbors(X(:,jj2),Y,eps,opts2);
            ret=ret1+ret2;
            return;
        else
            ret=zeros(1,NY);
            for bb=ii1:ii2
                inds_X=find(X_bin_indices==bb);
                inds_Y=find((Y_bin_indices>=bb-1)&(Y_bin_indices<=bb+1));
                if ((length(inds_X)>0)&&(length(inds_Y)>0))
                    ret(inds_Y)=ret(inds_Y)+count_neighbors(X(:,inds_X),Y(:,inds_Y),eps,opts2);
                end;
            end;
            return;
        end;
    end;
end;

function ret=exhaustive_count_neighbors(X,Y,eps)
[M,NX]=size(X);
[M,NY]=size(Y);
distsqrs=zeros(NX,NY);
for m=1:M
    distsqrs=distsqrs+(repmat(X(m,:),NY,1)'-repmat(Y(m,:),NX,1)).^2;
end;
tmp=(distsqrs<=eps^2);
ret=sum(tmp,1);

function test_count_neighbors

close all;


M=20;
N=6000;

X=randn(M,N);

%[xx,yy]=ndgrid(1:3,1:3);
%X=zeros(2,length(xx(:)));
%X(1,:)=xx(:);
%X(2,:)=yy(:);

eps=3;

tic;
counts=count_neighbors(X,X,eps,struct('max_num_trials',4));
toc

%counts
%exhaustive_count_neighbors(X,X,1)

tic;
counts_check=count_neighbors(X,X,eps,struct('max_num_trials',24));
toc
errs=(counts-counts_check)./(counts+counts_check);
figure; hist(errs,100);

figure; hist(counts,100);

