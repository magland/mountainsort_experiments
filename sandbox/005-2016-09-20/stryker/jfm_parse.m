function jfm_parse

raw_path='/home/magland/prvdata/stryker/2015_12_10';

mkdir([raw_path,'/raw_mda']);

G=load('emap64D.mat');
csvwrite('geom.csv',G.emap(2:end-1,:));

files=dir([raw_path,'/data003']);
arrays=cell(0,1);

%for j=3:length(files)
for j=3:(3+10-1)
    X=read_Intan_RHD2000_file([raw_path,'/data003/',files(j).name]);
    arrays{end+1}=X.amplifier_data;
end;

X=cat_timeseries(arrays);

writemda32(X,[raw_path,'/raw_mda/concat_first10.mda']);

end

function Y=cat_timeseries(X)
M=size(X{1},1);
N=0;
for j=1:length(X)
    N=N+size(X{j},2);
end;
Y=zeros(M,N);
ii=1;
for j=1:length(X)
    Y(:,ii:ii+size(X{j},2)-1)=X{j};
    ii=ii+size(X{j},2);
end;

end