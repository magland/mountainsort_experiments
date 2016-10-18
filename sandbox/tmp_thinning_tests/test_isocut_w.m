function test_isocut_w

close all;

N0=1000;

X1=randn(1,N0);
X2=randn(1,N0)+5;
X=cat(2,X1,X2);
X=sort(X);
figure; hist(X,200);

[X_thin_inds,weights]=thin(X,500);
X_thin=X(1,X_thin_inds);

figure; hist(X_thin,100);
title('X thin');
figure; plot(weights);
title('weights');

cutpoint=isocut_w_dev(X_thin,weights,1.5)

function cutpoint=isocut_w_dev(X,weights,thresh)

cutpoint=0;

[X,inds]=sort(X);
weights=weights(inds);

spacings_weights=(weights(1:end-1)+weights(2:end))/2;

spacings=X(2:end)-X(1:end-1);
spacings=spacings./spacings_weights;
figure; plot(spacings);
title('spacings');

spacings_fit=jisotonic(spacings,'downup');
figure; plot(spacings_fit);
title('spacings fit');

spacings_resid=spacings./spacings_fit;
figure; plot(spacings_resid);
title('spacings resid');

spacings_resid_fit=jisotonic(spacings_resid,'updown');
figure; plot(spacings_resid_fit);
title('spacings resid fit');