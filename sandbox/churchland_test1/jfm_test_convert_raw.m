function jfm_test_convert_raw

addpath('/home/magland/Dropbox (Simons Foundation)/misc/nlxtomatlab/releaseDec2015');
addpath('/home/magland/Dropbox (Simons Foundation)/misc/nlxtomatlab/releaseDec2015/binaries');

csc_fname='/home/magland/Dropbox (Simons Foundation)/Shared neural data/CSC26.ncs';
[ts,X]=getRawCSCData(csc_fname,1,inf,4);
%X=reshape(X,512,length(X)/512);

writemda16i(X','/home/magland/prvdata/churchland/CSC26.ncs.mda');