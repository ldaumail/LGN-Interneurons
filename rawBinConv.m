%This script enables one to save neurophysological data in the requested
%raw binary int16 format for spike sorting with Kilosort 3.
%Loic Daumail -12/5/2021
%{
%list of modulated units and channels 
allUnits =  {'x160616_I_p01_uclust22_cinterocdrft_stab_fft_sig'}
    {'x160616_I_p01_uclust4_cinterocdrft_stab_fft_sig' }
    {'x160623_I_p04_uclust86_cinterocdrft_stab_fft_sig'}
    {'x160623_I_p05_uclust94_cinterocdrft_stab_fft_sig'}
    {'x160628_I_p01_uclust1_cinterocdrft_stab_fft_sig' }
    {'x160628_I_p01_uclust47_cinterocdrft_stab_fft_sig'}
    {'x180827_I_p02_uclust8_cinterocdrft_stab_fft_sig' }
    {'x181212_B_p01_uclust7_cinterocdrft_stab_fft_sig' }
    {'x190119_B_p01_uclust2_cinterocdrft_stab_fft_sig' }
    {'x190124_B_p01_uclust4_cinterocdrft_stab_fft_sig' }
    {'x190205_B_p01_uclust14_cinterocdrft_stab_fft_sig'};

allChannels = {'eD16'}
    {'eD22'}
    {'eD07'}
    {'eD07'}
    {'eD12'}
    {'eD21'}
    {'eD12'}
    {'eC14'}
    {'eC18'}
    {'eC23'}
    {'eC13'}
%}

%modulxFilenames =  {'x160616_I','x160623_I','x160628_I', ...
%    'x180827_I','x181212_B', 'x190119_B' ,'x190124_B' ,'x190205_B'}; %only
%    include sessions with binocularly modulated units
allfilename = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\all_cont_lgn_data_pooled_09292021';
selected_data = load(strcat(allfilename, '.mat'));
selected_data =selected_data.selected_data_pooled;
xFilenames = fieldnames(selected_data);
modulxFilenames = xFilenames; %all sessions, quicker to substitute with same variable name
for i = 1:length(modulxFilenames)
selectDate = modulxFilenames{i};
session = selectDate(2:9);
%session = '180827_I';
dataDir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\ns6_all', '\', sprintf('%s',session), '\');
ns6Filenames = dir(strcat(dataDir, '*cinterocdrft*.ns6'));

ns6_filename = ns6Filenames(1).name;
clear NS_header banks neural
% Read in NS Header
NS_Header    = openNSx(strcat(dataDir,ns6_filename),'noread');
if strcmp(session(end), 'I')
    el = 'eD';
elseif strcmp(session(end), 'B')
    el  = 'eC';
end

% get basic info about recorded data
neural       = strcmp({NS_Header.ElectrodesInfo.ConnectorBank},el(2)); % logicals where contact bank name matches electrode of interest
N.neural     = sum(neural); % number of neural channels
NeuralLabels = {NS_Header.ElectrodesInfo(neural).Label}; %get labels
Fs           = NS_Header.MetaTags.SamplingFreq; % get sampling frequency
%nyq          = Fs/2;
%r2           = Fs/30000;

% open data for all channels.
clear NS DAT
c = find(neural == 1);
nchan(i) = length(c);

electrode = sprintf('c:%u:%u',c(1),c(end));
NS        = openNSx(strcat(dataDir,ns6_filename),electrode,'read','uV');
DAT       = NS.Data; NS.Data = [];  % this is the whole signal on all channels, 30 kHz!

ns6Chunk = DAT;
clear DAT

% sort data from top of electrode to bottom. (electrodes aren't in
% the right order)
% get indices


idx = zeros(1,length(NeuralLabels));
for j = 1:length(NeuralLabels)
    Str  = cell2mat(NeuralLabels(j));
    Key   = el;
    Str(strfind(Str, '%02d')) = [];
    Index = strfind(Str, Key);
    idx(1, j) = sscanf(Str(Index(1) + length(Key):end), '%g', 1);
end
Srt_ns6Chunk = nan(length(NeuralLabels),length(ns6Chunk(1,:)));
sortedLabels = cell(1, length(NeuralLabels));
%apply sorting
Srt_ns6Chunk(idx,:) = ns6Chunk(:,:);
sortedLabels(idx) = NeuralLabels;

%% filtering step for some sessions
Fb=[800 3000];
filtDat = nan(size(Srt_ns6Chunk));
for c =1:size(Srt_ns6Chunk,1)
    filtDat(c,:) = bandpass(Srt_ns6Chunk(c,:),Fb,Fs);
end
%%
%datI = int16(Srt_ns6Chunk);
datI = int16(filtDat);

% This dataset is 24 channels by 30000 time samples.
saveDir =strcat( 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\kilosort_binary\all_sessions','\', sprintf('%s',session));
mkdir(saveDir);
fid = fopen(strcat(saveDir,'\',sprintf('%s_orig2_filt.bin',ns6_filename(1:end-4))), 'w'); 
fwrite(fid, datI, 'int16');
fclose(fid);
%}
end

%tSrt_ns6Chunk = bsxfun(@minus,Srt_ns6Chunk' , mean(Srt_ns6Chunk,1)');
%save(strcat(dataDir,sprintf('%s.dat',ns6_filename(1:end-4))), 'reshDat','-v7.3');
%{
fid = fopen(strcat(dataDir,sprintf('%s.bin',ns6_filename(1:end-4))));
oneSlice = fread(fid, [24 60000], '*int16'); % It will be a 2D array after this.
fclose(fid);

figure();
plot(Srt_ns6Chunk(1,1:5000))
%% Load data with Blackrock function
modulxFilenames = {'x160628_I'};
chan = {'eD13'};
selectDate = modulxFilenames{1};
session = selectDate(2:end);
dataDir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\ns6_bino_modul', '\', sprintf('%s',session), '\');
ns6Filenames = dir(strcat(dataDir, '*cinterocdrft*.ns6'));
ns6_filename = ns6Filenames(1).name;
data = openNSxHL(strcat(dataDir,ns6_filename));

reshdat =reshape(data(1:24*floor(size(data,1)/24)), [24 floor(size(data,1)/24)]);

datafile = strcat(dataDir,sprintf('%s_ns6_mat3.bin',ns6_filename(1:end-4)));
FIDw = fopen(datafile, 'w+', 'ieee-le');
fwrite(FIDw, reshdat, 'int16');
fclose(FIDw);


datafile = strcat(dataDir,sprintf('%s.bin',ns6_filename(1:end-4)));
fid = fopen(datafile, 'w'); 
fwrite(fid, reshdat, 'int16');
fclose(fid);


%% load data with hard code
modulxFilenames = {'x160628_I'};
chan = {'eD13'};
selectDate = modulxFilenames{1};
session = selectDate(2:end);
dataDir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\ns6_bino_modul', '\', sprintf('%s',session), '\');
ns6Filenames = dir(strcat(dataDir, '*cinterocdrft*.ns6'));
ns6_filename = ns6Filenames(1).name;

fname =strcat(dataDir,ns6_filename);
[path,fname, fext] = fileparts(fname);
FID          = fopen([path '\' fname '.ns6'], 'r', 'ieee-le');
FileTypeID   = fread(FID, [1,8]   , '*char');
dataHeaderBytes = 9;
BasicHeader   = fread(FID, 306, '*uint8');
HeaderBytes   = double(typecast(BasicHeader(3:6), 'uint32')) + dataHeaderBytes;
fseek(FID, HeaderBytes, 'bof');
simpleData = fread(FID, inf, '*int16');

figure();
plot(simpleData(1:50000))
%}
