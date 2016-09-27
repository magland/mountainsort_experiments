function jfm_parse

mfile_path=fileparts(mfilename('fullpath'))

raw_path = '/home/magland/prvdata/EJ/'

jar_file = strcat(raw_path, 'Vision/Vision.jar')
javaaddpath(jar_file)

add_path = '2005-04-26-0/data002/'
full_path = strcat(raw_path, add_path)
rawFile = edu.ucsc.neurobiology.vision.io.RawDataFile(full_path)

% rawFile = edu.ucsc.neurobiology.vision.io.RawDataFile( \
%   '/usr/local/google/home/bhaishahster/Downloads/data002/');

% One can check dataset details (duration, number of electrodes, etc):
% the following :
% rawFile.getHeader()
% the geometry of the electrodes is in 'Electrode positions' folder.

% gets data from sample 0 to 20000 (first second) across all electrodes
time_len  = 60*20; % in seconds
samplingRate = 20000;
%channels=[131,132,138,139,140,147,148];
channels=101:120;

raw=zeros(length(channels),samplingRate*time_len);
for s=1:time_len
    fprintf('Parsing raw data: %d / %d seconds\n',s,time_len);
    data_len = samplingRate * 1; % 1 = one second
    data = rawFile.getData(data_len*(s-1)+1, data_len); % start time, length of datasets
    raw(:,data_len*(s-1)+1:data_len*s)=data(:,channels)';
end

rawFile.close();

outpath=[raw_path,'/raw_mda/2005-04-26-0/data002'];
outpath2=[mfile_path,'/datasets/2005-04-26-0/data002'];
mkdir(outpath);
mkdir(outpath2);
writemda(raw,[outpath,'/raw.mda'],'uint16');
system(sprintf('prv-create %s/raw.mda %s/raw.mda.prv',outpath,outpath2));
write_text_file([outpath2,'/params.json'],'{"samplerate":20000}');

geom=csvread('512coords.txt');
geom=geom(channels,:);
csvwrite([outpath2,'/geom.csv'],geom);


% Stimulus TTLs (most times every 100 frames) on java electrode
% position 0 (1 in matlab) other electrodes are raw voltage traces
%figure;
%plot(data(:,1:3));

end

function write_text_file(fname,txt)

fid = fopen(fname,'wt');
fprintf(fid, '%s',txt);
fclose(fid);

end
