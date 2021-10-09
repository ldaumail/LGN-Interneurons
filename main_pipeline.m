% This script attempts to process single unit data from monocular and
% binocular stimulation at multiple contrast levels to detect binocularly modulated single
% units. The goal is then to look at the corresponding multi-unit activity
% surrounding those modulated channels, as modulated single units may be in
% fact modulated by interneurons (inhibitory neurons) that could be
% detected in the MUA
% Written by Loic Daumail, last updated on 09/30/2021


%% Pre-start: use the list of single units file names that were selected in the adaptation analysis with
% high contrast
unitsDir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\inverted_power_channels\good_single_units_data_4bumps_more\new_peak_alignment_anal\';
unitsDataDir = [unitsDir 'refined_dataset']; 
unitsData= load(unitsDataDir);

selectUnitsFilenames =load(strcat(unitsDir,'single_units_ns6_metadata.mat'));
filenames = selectUnitsFilenames.STIMFileName;

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
codeDir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\loic_code_042021\lgn_interneuron_suppression';
addpath(genpath(codeDir))

selFilenames = modulxFilenames(8:11); %also need to do #8 but there is a bug
selChan = modulChannel(8:11);
rectMUA = rectifiedMUA(selFilenames, selChan);

%mkdir(strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\mua_bino_modul\', selectDate))    
%filtdir =strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\mua_bino_modul\', selectDate);
%save(strcat(filtdir,'_cinterocdrft_750hzhighpass_1khz_MUA.mat'), 'stackedMUA','-v7.3');

%% Trigger data to event codes and retrigger to photodiode

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
relpene = unique(penetrations);
dataDir =  'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\';
[all_photoMUA] = photoTrigMUA(relpene(9), dataDir);


%% Plot multiunit corresponding to each single unit

%%(1)load single units, take data used to perform ROC analysis of power at 4Hz
allfilename = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\all_cont_lgn_data_pooled_09292021';
single_units = load(strcat(allfilename, '.mat'));

%%(2)get single unit trials indices
tridx = struct();
for i = 1:length(sigBinoModulIdx)
    tridx.(modulxFilenames{i}) = single_units.selected_data_pooled.(modulxFilenames{i}).trials;
end
%%(3)select MUA trials
dataDir =  'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\';
photoTrigMuaDir = strcat(dataDir, 'photoReTrig_mua\');
muaFilenames = dir(strcat(photoTrigMuaDir, '*cinterocdrft*'));
mua = struct();
for i = 1:length(modulxFilenames) %modulated sua filenames
    xFilename = modulxFilenames{i};
    selectDate = xFilename(2:9);
    for j = 1:length({muaFilenames.name}) %potentially corresponding mua files
        if strfind(muaFilenames(j).name,selectDate)
            MUAfile = load(strcat(photoTrigMuaDir, muaFilenames(j).name));
            segNames = fieldnames(MUAfile.photoDiodeMUA);
            all_MUA = [];
            all_MUA_trials = [];
            all_MUA_photo_on = [];
            for n = 1:length(segNames)
                all_MUA = cat(3,all_MUA, MUAfile.photoDiodeMUA.(segNames{n}).aMUA);
                all_MUA_trials =cat(1,all_MUA_trials, MUAfile.photoDiodeMUA.(segNames{n}).trial);
                all_MUA_photo_on = {all_MUA_trials; MUAfile.photoDiodeMUA.(segNames{n}).photo_on};
            end
            stimCond = fieldnames(tridx.(xFilename));
            for bin =1:2
                mua.(xFilename).mua.(stimCond{bin}) = all_MUA(:,:,tridx.(xFilename).(stimCond{bin}) );
            end
        end
    end
end
    
    
%%(4)plot MUA
xFilename = modulxFilenames{4};
chan = modulChannel{4};
mean_MUA_mono =nanmean(mua.(xFilename).mua.(stimCond{2}),3);
mean_MUA_bino = nanmean(mua.(xFilename).mua.(stimCond{1}),3);

dmeans = mean_MUA_mono - mean_MUA_bino;
figure();
imagesc(dmeans')
 colorbar;
title(sprintf('MUA %s channel modulated %s',xFilename, chan), 'interpreter', 'none')

 saveDir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\analysis\';
 saveas(gcf,strcat(saveDir,'\',sprintf('mua_%s_%s_contacts.png', xFilename, chan)));

 
 fint_data = filterCSD(dmeans')';

h = figure();
 xabs = [-600:1300];
 yvec = [1:24];
 imagesc(xabs,yvec, fint_data')
title(sprintf('MUA %s channel modulated %s',xFilename, chan), 'interpreter', 'none')
 colorbar;
 saveDir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\analysis\';
 saveas(gcf,strcat(saveDir,'\',sprintf('mua_%s_%s.png', xFilename, chan)));
 %climit = max(abs(get(gca,'CLim'))*.8);
 %{
 if any(strfind(inputname(1), 'CSD'))
 set(gca,'CLim',[-800 800],'Box','off','TickDir','out')
 else
 end
 %}
 hold on;
plot([0 0], ylim,'k')
plot([1100 1100], ylim,'k')
clrbar = colorbar;

