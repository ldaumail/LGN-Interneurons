function [rectMUA] = rectifiedMUA(modulxFilenames, chan)
%Loic Daumail 
%this function aims to extrac the rectified multi-unit activity from our
%data. data is first high pass filtered at 750 Hz, then low pass filtered
%at 5000 Hz, rectified,and finally low pass filtered at 200 Hz.


rectMUA = struct();
for i = 1:length(modulxFilenames)
    selectxFile = modulxFilenames{i}; %try with one file
    selectDate = selectxFile(2:9);
    dataDir = strcat('C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\ns6_bino_modul\',selectDate,'\');
    ns6Filenames = dir(strcat(dataDir, '*cinterocdrft*.ns6'));
    
    curChan = char(chan(i));
    stackedMUA = struct();
    for fn = 1:length({ns6Filenames.name})
        ns6_filename = ns6Filenames(fn).name;
        clear ext NS_header banks neural
        % Read in NS Header
        ext          = 'ns6';
        NS_Header    = openNSx(strcat(dataDir,ns6_filename),'noread');
        
        el  = curChan(1:2);
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
                NS        = openNSx(strcat(dataDir,ns6_filename),electrode,'read','uV');
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
            Key   = curChan(1:2);
            Str(strfind(Str, '%02d')) = [];
            
            Index = strfind(Str, Key);
            idx(1, i) = sscanf(Str(Index(1) + length(Key):end), '%g', 1);
            
        end
        Srt_MUA = nan(length(MUA(:,1)), length(NeuralLabels));
        %apply sorting
        Srt_MUA(:,idx) = MUA(:,:);
        sortedLabels = NeuralLabels(idx);
        
        mua_filename = strcat('x',ns6_filename(1:end-4));
        stackedMUA.(mua_filename).MUA = Srt_MUA;
        stackedMUA.(mua_filename).sortedLabels = sortedLabels;
        
    end
muaDir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\mua_bino_modul\';
save(strcat(muaDir,sprintf('%s_cinterocdrft_750hzhighpass_1khz_MUA.mat',selectxFile(2:13))), 'stackedMUA','-v7.3');
 rectMUA.(selectxFile) = stackedMUA; 
end
end