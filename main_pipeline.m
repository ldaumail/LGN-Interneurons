% This script attempts to process single unit data from monocular and
% binocular stimulation at multiple contrast levels to detect binocularly modulated single
% units. The goal is then to look at the corresponding multi-unit activity
% surrounding those modulated channels, as modulated single units may be in
% fact modulated by interneurons (inhibitory neurons) that could be
% detected in the MUA
% Written by Loic Daumail, last updated on 09/30/2021


%% Pre-start: use the list of single units file names that were selected in the adaptation analysis with
% high contrast
selectUnitsFilenames =load('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\s_potentials_analysis\analysis\single_units_ns6_metadata.mat');
filenames = selectUnitsFilenames.STIMFileName;


unitsDir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\';
unitsDataDir = [unitsDir 'refined_dataset']; 
unitsData= load(unitsDataDir);


cellClass = {'K','M','P','K','K','K','M','P','P','','M','M','','','M','','','P','','M','','M','M','','P','M','','P', ...
'P','','','K','P','M','M','M','P','','P','K','P','P','','P','P','M','','P','M','P','M','P','','P','M','M','P','','M','M','P','M', ...
'','','M','M','M','P','M','M','M','M','P','P'};
cellClass([1,46,55]) = [];

directory = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\loic_code_042021\lgn_interneuron_suppression';
addpath(genpath(directory))

%% contrast levels encountered in the dominant eye

 allContLevels =[];
 for i =1:71
     if ~isempty(filenames{i})
         contLevels = unique(unitsData.new_data(i).channel_data.contrast);
         allContLevels = unique([allContLevels; contLevels]);
     end
 end
 
 
%% I. Get selected trials data

selected_data = suaTrialSelection(unitsData, filenames);%get 

allfilename = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\all_cont_lgn_data_09292021';
save(strcat(allfilename, '.mat'), 'selected_data');

%% II. Test for binocular modulation

%relevant comparisons:
%(]0:0.1] NDE && ]0:0.1] DE) vs (0 NDE && ]0.01]DE)
%(]0.1:0.3] NDE && ]0:0.1] DE) vs (0 NDE && ]0.01]DE)
%(]0.3:0.5] NDE && ]0:0.1] DE) vs (0 NDE && ]0.01]DE)
%(]0.7:1] NDE && ]0:0.1] DE) vs (0 NDE && ]0.01]DE)

%Repeat for ]0.1:0.3] DE, ]0.3:0.5] DE,]0.7:1] DE ===> total of 16 tests


lowContNDE = [0, 10, 30, 70];
lowContDE = [0, 10, 30, 70];

%binoconds:
binoConds = combvec(lowContDE,lowContNDE);
monoConds = binoConds(1,:);

sigResults = binocularModulationPowerTest(selected_data, binoConds, monoConds);
allfilename = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\binmod_lgn_power_test_09292021';
save(strcat(allfilename, '.mat'), 'sigResults');

%% I'. & II'.Try with more pooling of contrast levels to save trials
selected_data_pooled = suaTrialSelectionContPool(unitsData, filenames);%get 
allfilename = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\all_cont_lgn_data_pooled_09292021';
save(strcat(allfilename, '.mat'), 'selected_data_pooled');

lowContNDE = [50];
lowContDE = [0];

%binoconds:
binoConds = combvec(lowContDE,lowContNDE);
monoConds = binoConds(1,:);

sigPooledResults = binocularModulationPowerTest(selected_data_pooled, binoConds, monoConds);
allfilename = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\binmod_lgn_pooled_power_test_09292021';
save(strcat(allfilename, '.mat'), 'sigPooledResults');
%% Plot mean response of modulated data 

%% III. Get .ns6 files and .nev files for the MUA that corresponds to the modulated SUA

%get modulated single units indices
allfilename = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\binmod_lgn_pooled_power_test_09292021';
sigPooledResults = load(strcat(allfilename, '.mat'));
sigBinoModul = sigPooledResults.sigPooledResults;
sigBinoModul(isnan(sigBinoModul)) = 0;
sigBinoModulIdx = find(sigBinoModul);

%get modulated single units filenames and channel number
allfilename = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\all_cont_lgn_data_pooled_09292021';
selected_data = load(strcat(allfilename, '.mat'));
selected_data =selected_data.selected_data_pooled;
xfilenames = fieldnames(selected_data);
modulxFilenames = xfilenames(sigBinoModulIdx);

channel = cell(length(xfilenames),1);
for i = 1:length(xfilenames)
    channel{i} = selected_data.(xfilenames{i}).chan; 
end
modulChannel = channel(sigBinoModulIdx);

%Load MUA and select corresponding MUA
%{
%MUA files were obtained with "mua_adaptation_analysis.m" 
muaChannelsdir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\multi_units\selected_channels\';
muaFilenames = dir(strcat(muaChannelsdir, '*_DE*')); %only include the DE files, as contrast range differs between eyes stimulated. If we use the CRFs to label channels, we will need a wide range of contrasts ==> "DE' files

%create data structure to select trial responses
unitsData = struct();

for i = 1: length(muaFilenames)
    filename = muaFilenames(i).name;
    MUA_datafile = load(strcat(muaChannelsdir, filename));
    unitsData.new_data(i).channel_data = MUA_datafile;
end
%}

%transfer necessary .ns6 files from TEBA onto main disk
for i = 2:length(modulxFilenames)
    
    selectxFile = modulxFilenames{i};
    selectDate = selectxFile(2:9);
    %ns6Dir = strcat('E:\LGN_data_TEBA\rig021\', selectDate);
    %ns6Dir = strcat('E:\LGN_data_TEBA\rig021_2\', selectDate);
    %ns6Dir = strcat('E:\LGN_data_TEBA\rig022_2\', selectDate);
     ns6Dir = strcat('E:\LGN_data_TEBA\rig022\', selectDate);
    if exist(ns6Dir, 'dir')
        mkdir(strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\ns6_bino_modul\', selectDate))
        sourcedir = ns6Dir;
        destdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\ns6_bino_modul\',selectDate);
    end
    [status, msg] = copyfile( strcat(sourcedir,'\*cinterocdrft*.ns6'),  destdir);
end


%% Process .ns6 files to get the MUA and the right trial indices

%.ns6 file directory
selectxFile = modulxFilenames{1};
selectDate = selectxFile(2:9);
dataDir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\ns6_bino_modul\',selectDate,'\');

ns6Filenames = dir(strcat(dataDir, '*cinterocdrft*'));
for fn = 1:size(ns6Filenames.name)
    ns6_filename = ns6Filenames(fn).name;
    clear ext NS_header banks neural
    % Read in NS Header
    ext          = 'ns6';
    NS_Header    = openNSx(strcat(dataDir,ns6_filename),'noread');
    el  = 'eD';
    % get basic info about recorded data
    neural       = strcmp({NS_Header.ElectrodesInfo.ConnectorBank},el(2)); % logicals where contact bank name matches electrode of interest
    N.neural     = sum(neural); % number of neural channels
    NeuralLabels = {NS_Header.ElectrodesInfo(neural).Label}; %get labels
    Fs           = NS_Header.MetaTags.SamplingFreq; % get sampling frequency
    nyq          = Fs/2;
    r            = Fs/1000;
    %r2           = Fs/30000;

    % counters
    clear nct
    nct = 0;

    tic
    % process data electrode by electrode
    for e = 1:length(neural)
        if neural(e) == 1    % why? because neural is a vector of logicals, 1 = contacts we want
            nct = nct+1;

            % open data for this channel.
            clear NS DAT
            electrode = sprintf('c:%u',e);
            NS        = openNSx(strcat(dataDir,ns6_filename),electrode,'read');%'uV');
            DAT       = NS.Data; NS.Data = [];  % this is the whole signal on one channel, 30 kHz!


            % preallocate data matrices
            if nct == 1
                N.samples = length(DAT);
                %LFP       = zeros(ceil(N.samples/r),N.neural); % preallocating for downsampled data
                MUA      = zeros(ceil(N.samples/r),N.neural);
            end
        % extract the MUA:
            clear hpc hWn bwb bwa hpMUA
            hpc       = 750;  %high pass cutoff
            hWn       = hpc/nyq;
            [bwb,bwa] = butter(4,hWn,'high');
            hpMUA     = filtfilt(bwb,bwa,double(DAT)); %high pass filter

            % low pass at 5000 Hz and rectify 
            clear lpc lWn bwb bwa 
            lpc       = 5000;  % cutoff
            lWn       = lpc/nyq;
            [bwb,bwa] = butter(4,lWn,'low');
            hpMUA     = abs(filtfilt(bwb,bwa,hpMUA)); %low pass filter &rectify

            % low pass filter at x Hz. 
            clear lpc lWn bwb bwa lpMUA
            lpc       = 200; %low pass cutoff
            lWn       = lpc/nyq;
            [bwb,bwa] = butter(4,lWn,'low'); 
            lpMUA     = filtfilt(bwb,bwa,hpMUA);  %low pass filter to smooth


            % decimate both LFP and analog MUA (aMUA) to get 1kHz samp freq
            MUA(:,nct) = decimate(lpMUA,r); 
            
            % high pass filter
             %{
            clear hpc hWn bwb bwa hpMUA
            hpc       = 250; %high pass cutoff
            hWn       = hpc/nyq;
            [bwb,bwa] = butter(4,hWn,'high');
            hpsLFP      = filtfilt(bwb,bwa,DAT);  %low pass filter
            %}
            % decimate analog MUA (aMUA) to get 1kHz samp freq
            % MUA(:,nct) = decimate(lpMUA,r);

            % decimate sLFP to get 20kHz samp freq
            %LFP(:,nct) = decimate(hpsLFP,r2);

            clear DAT

        end

    end
    toc
    % Warning! THESE DATA ARE IN THE SAME ORDER AS THE BR PINS, NOT THE ORDER OF THE PROBE
    %As the .ns6 data was retrieved using openNSx() and not getLFP(), we need
    %to sort the channels ourselves. With getLFP(), sorting is done
    %automatically
    % sort data from top of electrode to bottom.

    % get indices
    idx = zeros(1,length(NeuralLabels));
    for i = 1:length(NeuralLabels)

        Str  = cell2mat(NeuralLabels(i));
        Key   = 'eD';
        Str(strfind(Str, '%02d')) = [];

        Index = strfind(Str, Key);
        idx(1, i) = sscanf(Str(Index(1) + length(Key):end), '%g', 1);

    end
    Srt_MUA = nan(length(MUA(:,1)), length(NeuralLabels));
   %apply sorting
    Srt_MUA(:,idx) = MUA(:,:);
    sortedLabels = NeuralLabels(idx);



filtdir ='C:\Users\daumail\Documents\LGN_data\single_units\s_potentials_analysis\data\filt_data';

% save(strcat(trigdir,'\',ns6_filename, '_rectified_1000hz_aMUA.mat'), 'STIM_aMUA');
RelChan = Srt_sLFP(:,11); %store data of the relevant channel
save(strcat(filtdir,'\',ns6_filename, '_250hzhighpass_30khz_sLFP.mat'), 'RelChan','-v7.3');


