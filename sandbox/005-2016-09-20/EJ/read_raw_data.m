
mfile_path=fileparts(mfilename('fullpath'))

% Add path to Vision.jar (It can be found on Spike sorting datasets
% folder titled Vision)
main_path = '/disk1/prvdata/EJ/'
jar_file = strcat(main_path, 'Vision/Vision.jar')
% javaaddpath('/usr/local/google/home/bhaishahster/GITs/Java/vision7/Vision.jar');
javaaddpath(jar_file)

% Open raw data file give path to datasets folder (like data002/ - it will
% choose the corresponding binary automatically)
add_path = '2005-04-26-0/data002/'
full_path = strcat(main_path, add_path)
rawFile = edu.ucsc.neurobiology.vision.io.RawDataFile(full_path)

% rawFile = edu.ucsc.neurobiology.vision.io.RawDataFile( \
%   '/usr/local/google/home/bhaishahster/Downloads/data002/');

% One can check dataset details (duration, number of electrodes, etc):
% the following :
% rawFile.getHeader()
% the geometry of the electrodes is in 'Electrode positions' folder.

% gets data from sample 0 to 20000 (first second) across all electrodes
time_len  = 2; % in seconds
samplingRate = 20000;
data_len = samplingRate * time_len;
data = rawFile.getData(0, data_len); % start time, length of datasets

% Stimulus TTLs (most times every 100 frames) on java electrode
% position 0 (1 in matlab) other electrodes are raw voltage traces
figure;
plot(data(:,1:3));

rawFile.close();
