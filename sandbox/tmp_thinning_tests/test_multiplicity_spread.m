function test_multiplicity_spread

M=18;
num_points=1e7;
epsilon=1;

points=(rand(M,num_points)*2-1)*epsilon;
dists_from_zero=sqrt(sum(points.^2,1));

nearby_points=points(:,dists_from_zero<=epsilon);
vals=nearby_points(1,:);
figure; hist(vals,100);

tt=linspace(-1,1,1000);
figure; plot(tt,sqrt(1-tt.^2).^(M-1));