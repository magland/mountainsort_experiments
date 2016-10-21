function test_distance_estimate

M=40;
K=2; % dimension of random subspace
num_trials=10000;

random_rotations=randn(M,M,num_trials);
for j=1:num_trials
    random_rotations(:,:,j)=orth(random_rotations(:,:,j));
end;
%random_subspaces=random_subspaces./repmat(sqrt(sum(random_subspaces.^2,1)),M,1,1);

dists_from_zero_in_subspace=squeeze(sqrt(sum(random_rotations(1,1:K,:).^2,2)));
%dists_from_zero_in_subspace=squeeze(max(abs(random_subspaces(1,:,:)),[],2));

figure;
hist(dists_from_zero_in_subspace,100);
xlim0=xlim; xlim([0,xlim0(2)]);

figure;
hist(1./dists_from_zero_in_subspace,100);
