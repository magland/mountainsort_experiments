function jfm_parse

raw_path='/home/magland/prvdata/stryker/2015_12_10';
samplerate=20000;
duration=120; %seconds
dataset_name='first_2min';
num_timepoints=duration*samplerate;
channels=1:8;

mkdir([raw_path,'/raw_mda']);

G=load('emap64D.mat');
geom=G.emap(2:end-1,:);

files=dir([raw_path,'/data003']);
arrays=cell(0,1);

num_timepoints_read=0;
for j=3:length(files)
    if (num_timepoints_read<num_timepoints)
        X=read_Intan_RHD2000_file([raw_path,'/data003/',files(j).name]);
        X=X.amplifier_data;
        X=X(channels,:);
        num_timepoints_read=num_timepoints_read+size(X,2);
        arrays{end+1}=X;
    end;
end;

X=cat_timeseries(arrays);
if (size(X,2)>num_timepoints)
    X=X(:,1:num_timepoints);
end;

%writemda32(X,[raw_path,'/raw_mda/concat_first10.mda']);
writemda16ui(X,[raw_path,'/raw_mda/',dataset_name,'.mda']);
mkdir(['datasets/',dataset_name]);
csvwrite(['datasets/',dataset_name,'/geom.csv'],geom(channels,:));
system(sprintf('prv-create %s %s',[raw_path,'/raw_mda/',dataset_name,'.mda'],['datasets/',dataset_name,'/raw.mda.prv']));
write_text_file(['datasets/',dataset_name,'/params.json'],'{"samplerate":20000}');

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

function write_text_file(fname,txt)
F=fopen(fname,'w');
fprintf(F,'%s',txt);
fclose(F);
end
