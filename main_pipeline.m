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

%transfer necessary .ns6, .nev, .text (.gCINTEROCDRFTGrating_di) and .bhv files from TEBA onto main disk
PATHS = {'E:\LGN_data_TEBA\rig021\', 'E:\LGN_data_TEBA\rig021_2\', 'E:\LGN_data_TEBA\rig022\', 'E:\LGN_data_TEBA\rig022_2\'};
for p =1:length(PATHS)
    for i = 1:length(modulxFilenames)
        selectxFile = modulxFilenames{i};
        selectDate = selectxFile(2:9);
        ns6Dir = strcat(PATHS{p}, selectDate);
        if exist(ns6Dir, 'dir')
            %mkdir(strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\ns6_bino_modul\', selectDate))
            sourcedir = ns6Dir;
            destdir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\ns6_bino_modul\',selectDate);
            %[status, msg] = copyfile( strcat(sourcedir,'\*cinterocdrft*.ns6'),  destdir);
            %[status, msg] = copyfile( strcat(sourcedir,'\*cinterocdrft*.nev'),  destdir);
            %[status, msg] = copyfile( strcat(sourcedir,'\*.gCINTEROCDRFTGrating_di'),  destdir);
            [status, msg] = copyfile( strcat(sourcedir,'\*cinterocdrft*.bhv'),  destdir);
        end
    end
end

%% IV Process .ns6 files to get the MUA and the right trial indices

% (1) Get MUA from .ns6 files 
%.ns6 file directory
directory = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\loic_code_042021\lgn_interneuron_suppression';
addpath(genpath(directory))

selFilenames = modulxFilenames([4;9:11]); %also need to do #8 but there is a bug
rectMUA = rectifiedMUA(selFilenames, modulChannel);

%mkdir(strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\mua_bino_modul\', selectDate))    
%filtdir =strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\mua_bino_modul\', selectDate);
%save(strcat(filtdir,'_cinterocdrft_750hzhighpass_1khz_MUA.mat'), 'stackedMUA','-v7.3');

%%
%npmkdir    = 'C:\Users\maier\Documents\MATLAB\NPMK-master\'; 
%nbanadir   = 'C:\Users\maier\Documents\bootcamp-selected\nbanalysis\'; 
%addpath(genpath(dataDir))
%addpath(genpath(npmkdir))
%addpath(genpath(nbanadir))
penetrations = cell(length(modulxFilenames),1);
for i =1:length(modulxFilenames)
    fname = modulxFilenames{i};
    penetrations{i} = fname(2:13);
end
%for i = 1: length(unique(penetrations)) %number of sessions
selectFile = penetrations{1}; %try with one file, but use 'i' index otherwise !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
selectDate = selectFile(1:8);
dataDir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\ns6_bino_modul\',selectDate,'\');

textFilenames = dir(strcat(dataDir, '*.gCINTEROCDRFTGrating_di'));
%for n = 1:length(textFilenames.name)
BRdatafile = textFilenames(1).name(1:24); %try with first file use 'n' index !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
filename   = [dataDir BRdatafile]; 


%% V. LOAD STIMULUS CONDITIONS with text file 

patterns   = {'rforidrft','rfsfdrft','posdisparitydrft','disparitydrft','cinterocdrft','coneinterocdrft','conedrft', ...
                'colorflicker','bwflicker','rfori','rfsize','cinteroc','color','rfsf','mcosinteroc','dotmapping'}; 
clear p match
for p = 1:length(patterns)

   pattern      = patterns{p}; 
  
   if any(strfind(BRdatafile,pattern))
       startlog = strfind(BRdatafile,pattern); 
       if ~isequal(BRdatafile(startlog:end-3),pattern),continue
       else
       match    = patterns{p}; 
       end
   end
   
end

if isequal(match,'dotmapping')
ext  = '.gDotsXY_di';
else
ext  = ['.g' upper(match) 'Grating_di']; 
end

if contains(ext,'DRFT')
      grating     = readgDRFTGrating([filename ext]); % from nbanalysis (or even MLAnalysisOnline--might be out of date)
elseif contains(ext,'Dots')
      grating     = readgDotsXY([filename ext]);
else
      grating     = readgGrating([filename ext]);
end


%% VI. LOAD EVENT TIMES/CODES TO TRIGGER DATA
%The NEV files contain the digital event codes that are sent from MonkeyLogic to the 
%BlackRock machine. They are sent for different stimulus/behavioral events, 
%including when the stimuli appear on the screen. 
nevFilenames = dir(strcat(dataDir, '*cinterocdrft*.nev'));  
nevfilename = strcat(dataDir, nevFilenames(1).name); %use 'n' index !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
NEV             = openNEV(nevfilename ,'noread','overwrite');

%The codes themselves can found in this field after subtracting 128 (unknown quirk). 
%Good codes to know are 23, 25, 27, 29, 31 (Stimulus Onsets—23 
%is first stimulus on in a trial, 25 is second stimulus on in a trial, etc). 
%Stimulus offsets are the evens in that range—24, 26, 28, 30, 32. 
%The starts of trials are signaled by 999 and the ends by 18 18 18. 
%There are others in there as well. 

EventCodes      = NEV.Data.SerialDigitalIO.UnparsedData - 128;        % we don't know whywe have to subtract 128 but we do
%The time the code was sent can be found in this field. 
%This is in samples (at the 30kHz rate). The Times on the next line are in seconds. 
%We multiply by 1000 to get into ms. 
 
EventSamples    = NEV.Data.SerialDigitalIO.TimeStamp;                 % in samples 
EventTimes      = floor(NEV.Data.SerialDigitalIO.TimeStampSec.*1000); % convert to ms 

%Checkout the parsEventCodesML function. It takes those 999s and 18 18 18s and 
%breaks up the codes/times into separate trials.
%Rfori files typically contain 5 stimulus presentations per trial.

[pEvC, pEvT]    = parsEventCodesML(EventCodes,EventSamples);          % sorts codes, samps or times into trials


% So far all of these data are from EVERY trial, including trials where
% animal breaks fixation. Lets get rid of those and make a new structure
% with the grating info and the stimulus onsets 

STIM            = sortStimandTimeData(grating,pEvC,pEvT,'stim'); % this is in nbanalysis. definitely double check it before you use it. 

%% VII. Trigger MUA to behavioral event codes
muaDir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\mua_bino_modul\';
stacked_MUA = load(strcat(muaDir,sprintf('%s_cinterocdrft_750hzhighpass_1khz_MUA.mat',selectFile)));
muaFiles = fieldnames(stacked_MUA.stackedMUA);
pre   = -600;
post  = 1300; 
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
STIM.aMUA = trigData(stacked_MUA.stackedMUA.(muaFiles{1}).MUA,floor(STIM.onsets./30),-pre,post); % use 'n' index !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


%% VIII. re-trigger MUA to photodiode = use photoRetriggerSTIM..?
TP = repmat([pre post], length(STIM.trial),1);
STIM.photoMUA = photoReTrigger(TP, filename,STIM.ypos, STIM);







