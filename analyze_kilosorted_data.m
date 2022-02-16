%%This script intends to load the spikes data sorted by kilosort
%%to :
%- get spike features
%- spike trigger to MUA
%Loic Daumail 12/06/2021

%%1. Look at spike features


%modulxFilenames =  {'x160616_I','x160623_I','x160628_I', ...
%    'x180827_I','x181212_B','x190119_B' ,'x190124_B' , 'x190205_B'}; %
allfilename = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\all_cont_lgn_data_pooled_09292021';
selected_data = load(strcat(allfilename, '.mat'));
selected_data =selected_data.selected_data_pooled;
xFilenames = fieldnames(selected_data);

sessions = [];
for j =1:length(xFilenames)
    session= xFilenames{j};
    sessions{j} = session(2:9);
end
[selectSess, SelectIdx] = unique(sessions);
modulxFilenames =selectSess;
 waveforms = struct();
for i = 1:length(modulxFilenames)
    selectDate = modulxFilenames{i};
    session = selectDate(1:end);
    xselectDate = strcat('x',selectDate);
    dataDir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\kilosort_binary\all_sessions', '\', sprintf('%s',session), '\');
    rezNm = dir(strcat(dataDir,'\','rez.mat'));
    if sum(size(rezNm.name,1)) == 1
    %Load binary data
      fs          = dir(fullfile(dataDir, '*.bin')) ;
      fbinary = fullfile(dataDir, fs(1).name);
      fid = fopen(fbinary);
      binDat = fread(fid, 'int16');
      fclose(fid);
      %%should load channel map .npy to figure out number of channels
      fkiloInfo = fullfile(dataDir, 'rez.mat');%
      kiloInfo = load(fkiloInfo);%
      nchans = length(kiloInfo.rez.ops.chanMap);

      %fkiloInfo = fullfile(dataDir, 'channel_map.npy');
      %kiloInfo = readNPY(fkiloInfo );
     % nchans = length( kiloInfo );
      %reshBinDat = reshape(binDat, [nchans,length(binDat)/nchans]);
    %Load Spike times
      fspikeTimes = fullfile(dataDir, 'spike_times.npy');
      spikeTimes = readNPY(fspikeTimes);
      
    %spike templates
     fspikeTemplates = fullfile(dataDir, 'spike_templates.npy');
     spikeTemplates = readNPY(fspikeTemplates);
     
     %get spike waveforms
     TemplatesID = unique(spikeTemplates);
     
     %get templates best channel
     %fbestChan = fullfile(dataDir, 'channel_positions.npy'); not good

     %bestChan = readNPY(fbestChan);
     
      gwfparams.dataDir = dataDir;    % KiloSort/Phy output folder
      gwfparams.fileName = fs(1).name;         % .dat file containing the raw 
      gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
      gwfparams.nCh = nchans;                      % Number of channels that were streamed to disk in .dat file
      gwfparams.wfWin = [-40 41];              % Number of samples before and after spiketime to include in waveform
      gwfparams.nWf = 1000;                    % Number of waveforms per unit to pull out
     %
     %get spike waveforms
    
     gwfparams.spikeTimes =[];
     gwfparams.spikeClusters = [];
     for t =1:length(TemplatesID)
         if (length(spikeTimes(spikeTemplates == TemplatesID(t))) >= 1000) && (kiloInfo.rez.good(t) == 1) %if we have at least 1000 spikes per template and the quality of isolation is good
            
            templateSpikeTimes = spikeTimes(spikeTemplates == TemplatesID(t));
            gwfparams.spikeTimes =  [ gwfparams.spikeTimes; templateSpikeTimes] ; % Vector of cluster spike times (in samples) same length as .spikeClusters
            gwfparams.spikeClusters = [ gwfparams.spikeClusters ;TemplatesID(t)]; % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes
            
         end
     end
     waveforms.(xselectDate).('templatesWF') =  getWaveForms(gwfparams); %get waveform
     peakval = max(abs(waveforms.(xselectDate).('templatesWF').waveFormsMean),[],3);
     [~,idx] =max(peakval,[],2);
     waveforms.(xselectDate).('templatesWF').bestChannel = idx;
     for c =1:length(idx)
     waveforms.(xselectDate).('spikeWidth')(c) = computeSpikeWidths( squeeze(waveforms.(xselectDate).('templatesWF').waveFormsMean(:,idx(c),:)), 1); %[1:length(waveforms.(selectDate).('templatesWF').unitIDs)]
     end
     %figure();
    % plot(squeeze(abs(waveforms.(selectDate).('templatesWF').waveFormsMean(1,:,:))))
    end
end  


%%plot spike widths histogram:
widths =[];
wfs =[];
cnt = 0;
for i =1:length(modulxFilenames)
    selectDate = modulxFilenames{i};
    xselectDate = strcat('x',selectDate);
    if isfield(waveforms,(xselectDate))
        if isfield(waveforms.(xselectDate),('spikeWidth'))
            cnt = cnt+1;
            widths = [widths,  waveforms.(xselectDate).('spikeWidth').*1000./(30000)];
            chan = waveforms.(xselectDate).('templatesWF').bestChannel;
            for n = 1:length(chan)
                %if size(squeeze(mean(waveforms.(xselectDate).('templatesWF').waveFormsMean,2))',1) > size(squeeze(mean(waveforms.(xselectDate).('templatesWF').waveFormsMean,2))',2)
                    wf = squeeze(waveforms.(xselectDate).('templatesWF').waveFormsMean(n,chan(n),:));
                    normWf = wf./max(abs(wf));
                    wfs =  [wfs , wf];
               % elseif size(squeeze(mean(waveforms.(xselectDate).('templatesWF').waveFormsMean,2))',1) < size(squeeze(mean(waveforms.(xselectDate).('templatesWF').waveFormsMean,2))',2)
                 %   wf = squeeze(waveforms.(xselectDate).('templatesWF').waveFormsMean(n,chan(n),:));
                  %  normWf = wf./max(abs(wf));
                  %  wfs =  [wfs , wf];
                %end
                figure()
                subplot(1,2,1)
                plot(wf)
                subplot(1,2,2)
                histogram(abs(waveforms.(xselectDate).('spikeWidth').*1000./(30000)),100)
            end
        end
    end
end

figure();
subplot(1,2,1)
x = linspace(-82*1000/(30000*2),82*1000/(30000*2),82);
plot(x,wfs)
xlabel('Time(ms)')
ylabel('Voltage (microVolt)')

set(gca,'box','off')
set(gca,'linewidth',2)
subplot(1,2,2)
histogram(widths,100)
%xlim([0 0.5])
xlabel('Spike width(ms)')
ylabel('N units')
set(gca,'box','off')
set(gca,'linewidth',2)
sgtitle('Kilosorted spike waveforms from all sessions')
plotdir= 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\analysis\plots\';
saveas(gcf,strcat(plotdir, 'all_units_wf_rawspikewidths','.png'));
saveas(gcf,strcat(plotdir, 'all_units_wf_rawspikewidths','.svg'));

%xlim([1 15])
