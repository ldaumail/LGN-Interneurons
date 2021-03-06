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

%% III. Store .ns6 files and .nev files in specific folders for the MUA that corresponds to the modulated SUA

%get modulated single units indices
allfilename = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\binmod_lgn_pooled_power_test_09292021';
sigPooledResults = load(strcat(allfilename, '.mat'));
sigBinoModul = sigPooledResults.sigPooledResults;
sigBinoModul(isnan(sigBinoModul)) = 0;
sigBinoModulIdx = find(sigBinoModul);
facOrSup =sigBinoModul( sigBinoModul ~= 0);

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

%% IV. Process .ns6 files to store the MUA and the right trial indices

% (1) Get MUA from .ns6 files 
%.ns6 file directory
codeDir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\loic_code_042021\lgn_interneuron_suppression';
addpath(genpath(codeDir))

selFilenames = modulxFilenames(1); %also need to do #8 but there is a bug
selChan = modulChannel(1);
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
    penetrations{i} = fname(2:9);
end
relpene = unique(penetrations);
dataDir =  'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\';
[all_photoMUA] = photoTrigMUA(relpene(4:end), dataDir);


%% Plot multiunit (mono-bino) array corresponding to each single unit

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
        clear MUAfile
        if strfind(muaFilenames(j).name,selectDate)
            MUAfile = load(strcat(photoTrigMuaDir, muaFilenames(j).name));
            segNames = fieldnames(MUAfile.photoDiodeMUA);
            all_MUA = [];
            all_MUA_trials = [];
            all_MUA_photo_on = [];
            all_MUA_onsets = [];
            for n = 1:length(segNames)
                all_MUA = cat(3,all_MUA, MUAfile.photoDiodeMUA.(segNames{n}).aMUA);
                all_MUA_trials =cat(1,all_MUA_trials, MUAfile.photoDiodeMUA.(segNames{n}).trial);
                all_MUA_photo_on = horzcat(all_MUA_photo_on, MUAfile.photoDiodeMUA.(segNames{n}).photo_on);
                all_MUA_onsets = cat(1,all_MUA_onsets, MUAfile.photoDiodeMUA.(segNames{n}).onsets);
            end
            photo_on = nan(1,length(all_MUA_photo_on)); 
            for l = 1:length(all_MUA_photo_on) %only get first value of photo_on for each trial
                temp = all_MUA_photo_on{l};
                photo_on(l) = temp(1);
            end
            stimCond = fieldnames(tridx.(xFilename));
            for bin =1:2
                suaTrials = tridx.(xFilename).(stimCond{bin});

                suaToMuaTrials = suaTrials(suaTrials < size(all_MUA_trials,1)); %number of trials may be lower for mua, as there may be missing chunks of data
                mua.(xFilename).mua.(stimCond{bin}) = all_MUA(:,:,suaToMuaTrials);
                mua.(xFilename).photo_on.(stimCond{bin}) = photo_on(suaToMuaTrials);
                mua.(xFilename).onsets.(stimCond{bin}) = all_MUA_onsets(suaToMuaTrials);
                mua.(xFilename).time_diff.(stimCond{bin}) =  floor((photo_on(suaToMuaTrials)-all_MUA_onsets(suaToMuaTrials)')./30);

            end
        end
    end
end
 

%%(4) trigger mua to photodiode onset times !!! Currently not working
pre = -600;
post = 1300;

%%(5) trigger to photodiode onset times

for i = 1:length(modulxFilenames) %modulated sua filenames
    xFilename = modulxFilenames{i};
    selectDate = xFilename(2:9);
    stimCond = fieldnames(tridx.(xFilename));
            for bin =1:2

                tp =mua.(xFilename).time_diff.(stimCond{bin});
                origDat =mua.(xFilename).mua.(stimCond{bin});
                trigDat = nan(size(mua.(xFilename).mua.(stimCond{bin}),1)+max(mua.(xFilename).time_diff.(stimCond{bin})), size(mua.(xFilename).mua.(stimCond{bin}),2),size(mua.(xFilename).mua.(stimCond{bin}),3));
                for tr = 1:size(mua.(xFilename).mua.(stimCond{bin}),3)
                    trigDat(tp(tr):end,:,:) =origDat(:,:,tr);
                end
               % mua.(xFilename).trig_mua.(stimCond{bin}) = trigData(mua.(xFilename).mua.(stimCond{bin}),mua.(xFilename).time_diff.(stimCond{bin})+pre, pre, post);
               mua.(xFilename).trig_mua.(stimCond{bin}) = trigDat;
            end
end

%%(5)plot MUA
for i = 1:length(modulxFilenames)
    xFilename = modulxFilenames{i};
    chan = modulChannel{i};
    if facOrSup(i) == 1
        unitModul = 'facilitated';
    else
        unitModul = 'suppressed';
    end
    mean_MUA_mono =nanmean(mua.(xFilename).mua.(stimCond{2}),3); %!!! Not photo triggered
    mean_MUA_bino = nanmean(mua.(xFilename).mua.(stimCond{1}),3);%!!! Not photo triggered
    xabs = [-600:1300];
    yvec = [1:24];
    dmeans = mean_MUA_bino - mean_MUA_mono;
    figure();
    imagesc(xabs, yvec, dmeans')
    colorbar;
    set(gca, 'box', 'off');
    xlabel('time from stimulus onset(ms)');
    ylabel('contact #');
    title(sprintf('MUA diff (bino- mono) %s \n channel %s %s',xFilename, chan, unitModul), 'interpreter', 'none')
    
    saveDir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\analysis\plots\';
    saveas(gcf,strcat(saveDir,'\',sprintf('mua_%s_%s_contacts.png', xFilename, chan)));
    saveas(gcf,strcat(saveDir,'\',sprintf('mua_%s_%s_contacts.svg', xFilename, chan)));    
    
    fint_data = filterCSD(dmeans')';
    
    h = figure();
    
    imagesc(xabs,yvec, fint_data')

    title(sprintf('MUA diff (bino- mono) %s \n channel %s %s',xFilename, chan,unitModul), 'interpreter', 'none')
    colorbar;
    set(gca, 'box', 'off')
    xlabel('time from stimulus onset(ms)');
    ylabel('contact #');
    saveDir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\analysis\plots\';
    saveas(gcf,strcat(saveDir,'\',sprintf('mua_%s_%s.png', xFilename, chan)));
    saveas(gcf,strcat(saveDir,'\',sprintf('mua_%s_%s.svg', xFilename, chan)));
end
 
%% Spike-trigger MUA to SUA spikes
%Load spike data from BOSS
%Start with 1 session: 160616_I, chan 16
bossDataDir ='C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\ns6_bino_modul\160616_I\';
nevfilename ='160616_I_cinterocdrft_211111_001.nev';
bossData = loadNEV(strcat(bossDataDir,nevfilename), 'read');

%% 1) get spike timestamps from electrode 16
electrodesInfo = table(bossData.ElectrodesInfo);
el = 'eD22';
for i =1:size(electrodesInfo.Var1,2)
    if contains(erase(strjoin(string(bossData.ElectrodesInfo(i).ElectrodeLabel)),' '), el)
        electrodeID = electrodesInfo.Var1(1,i).ElectrodeID; 
        break
    end
end

unit1SpikeTimes = bossData.Data.Spikes.TimeStamp(bossData.Data.Spikes.Electrode == electrodeID &bossData.Data.Spikes.Unit ==1);
unit2SpikeTimes = bossData.Data.Spikes.TimeStamp(bossData.Data.Spikes.Electrode == electrodeID &bossData.Data.Spikes.Unit ==2);

%% Plot waveforms
unit1Waveform = 4.*double(bossData.Data.Spikes.Waveform(:,bossData.Data.Spikes.Electrode == electrodeID & bossData.Data.Spikes.Unit ==1));
unit2Waveform = 4.*double(bossData.Data.Spikes.Waveform(:,bossData.Data.Spikes.Electrode == electrodeID & bossData.Data.Spikes.Unit ==2));

meanUnit1Wf = mean(unit1Waveform,2);
meanUnit2Wf = mean(unit2Waveform,2);

u1_ci_low = meanUnit1Wf - std(unit1Waveform,[],2);
u1_ci_high = meanUnit1Wf + std(unit1Waveform,[],2);

u2_ci_low = meanUnit2Wf - std(unit2Waveform,[],2);
u2_ci_high = meanUnit2Wf + std(unit2Waveform,[],2);

x = linspace(-0.8,0.8,48);
figure();
subplot(2,1,1)
plot(x, meanUnit1Wf)
hold on
ciplot(u1_ci_low, u1_ci_high,x)
title('Unit 1')
set(gca, 'box','off')
subplot(2,1,2)
plot(x, meanUnit2Wf)
hold on
ciplot(u2_ci_low, u2_ci_high,x)
set(gca, 'box','off')
title('Unit 2')
xlabel('time (ms)')
ylabel('voltage (microVolt)')
plotdir= 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\analysis\plots\';
saveas(gcf,strcat(plotdir, sprintf('units_wf_%s_%s.png',nevfilename(1:8),el)));
saveas(gcf,strcat(plotdir, sprintf('units_wf_%s_%s.svg',nevfilename(1:8),el)));

%% 2) spike trigger MUA
%load corresponding MUA
muaDir =  'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\mua_bino_modul\';
muafilenames = dir(fullfile(muaDir,'*.mat'));
muafilename = muafilenames(2).name;
mua = load(strcat(muaDir, muafilename));

%1) spike trigger to unit 1 (looks more like a relay cell spike)
muaChunk = fieldnames(mua.stackedMUA);
relMua = mua.stackedMUA.(muaChunk{1}).MUA(:,:); %take surrounding channels (-3 +3)
sptrigMUA = nan(1651, length(unit1SpikeTimes(2:end-2)),size(relMua,2));
for c = 1:size(relMua,2)
    for j =2:length(unit1SpikeTimes(2:end-2))
        sptrigMUA(1:1651,j,c) = relMua(double(unit1SpikeTimes(j))-150:double(unit1SpikeTimes(j))+1500,c);
    end
end

meanSpMUA = squeeze(nanmean(sptrigMUA,2));
xabs = [-50:500];
yvec = [16-3:16+4];
fint_data = filterCSD(meanSpMUA(:,16-3:16+3)');
h = figure();
imagesc(xabs,yvec,fint_data)
xticks(linspace(-50,500,6))
%imagesc(meanSpMUA')
colorbar;
title(sprintf('Spike triggered MUA %s unit 1 on contact %s',nevfilename(1:8),el), 'interpreter', 'none')
plotdir= 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\analysis\plots\';
saveas(gcf,strcat(plotdir, sprintf('spiketrig_MUA_unit1_chan1319_%s_%s.png',nevfilename(1:8),el)));
saveas(gcf,strcat(plotdir, sprintf('spiketrig_MUA_unit1_chan1319_%s_%s.svg',nevfilename(1:8),el)));

%2) spike trigger to unit 2 (looks more like a relay cell spike)
muaChunk = fieldnames(mua.stackedMUA);
relMua = mua.stackedMUA.(muaChunk{1}).MUA(:,:); %take surrounding channels (-3 +3)
sptrigMUA = nan(1651, length(unit2SpikeTimes(2:end-2)),size(relMua,2));
for c = 1:size(relMua,2)
    for j =2:length(unit2SpikeTimes(2:end-2))
        sptrigMUA(1:1651,j,c) = relMua(double(unit2SpikeTimes(j))-150:double(unit2SpikeTimes(j))+1500,c);
    end
end

meanSpMUA = squeeze(nanmean(sptrigMUA,2));
xabs = [-50:500];
yvec = [1:24];
fint_data = filterCSD(meanSpMUA');
h = figure();
imagesc(xabs,yvec,fint_data)
xticks(linspace(-50,500,6))
%imagesc(meanSpMUA')
colorbar;
title(sprintf('Spike triggered MUA %s unit 2 on contact %s',nevfilename(1:8),el), 'interpreter', 'none')
plotdir= 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\analysis\plots\';
saveas(gcf,strcat(plotdir, sprintf('spiketrig_MUA_unit2_%s_%s.png',nevfilename(1:8),el)));
saveas(gcf,strcat(plotdir, sprintf('spiketrig_MUA_unit2_%s_%s.svg',nevfilename(1:8),el)));


