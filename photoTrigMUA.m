function [all_photoMUA] = photoTrigMUA(penetrations, dataDir)


ns6Dir = strcat(dataDir, 'ns6_bino_modul\');
muaDir = strcat(dataDir, 'mua_bino_modul\'); 
photoTrigMuaDir = strcat(dataDir, 'photoReTrig_mua\');

for i = 1: length(penetrations) %number of sessions
    selectFile = penetrations{i}; %try with one file, but use 'i' index otherwise !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    selectDate = selectFile(1:8);
    peneDir = strcat(ns6Dir,selectDate,'\');
    
    textFilenames = dir(strcat(peneDir, '*.gCINTEROCDRFTGrating_di'));
    for n = 1:length({textFilenames.name})
        BRdatafile = textFilenames(n).name(1:24); %try with first file use 'n' index !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        filename   = [peneDir BRdatafile];
        
        
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
        nevFilenames = dir(strcat(peneDir, '*cinterocdrft*.nev'));
        nevfilename = strcat(peneDir, nevFilenames(n).name); %use 'n' index !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
       
        stacked_MUA = load(strcat(muaDir,sprintf('%s_cinterocdrft_750hzhighpass_1khz_MUA.mat',selectFile)));
        muaFiles = fieldnames(stacked_MUA.stackedMUA);
        pre   = -600;
        post  = 1300;
        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        STIM.aMUA = trigData(stacked_MUA.stackedMUA.(muaFiles{n}).MUA,floor(STIM.onsets./30),-pre,post); % use 'n' index !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        %% VIII. re-trigger MUA to photodiode = use photoRetriggerSTIM..?
        TP = [STIM.trstart STIM.trend];
        [newTP, trigger, STIM]= photoReTrigger(TP, filename,STIM.ypos,{}, STIM);
        
        photoDiodeMUA.(muaFiles{n}) = STIM;
        
        save(strcat(photoTrigMuaDir,sprintf('%s_cinterocdrft_750hzhighpass_1khz_photoReTrigMUA.mat',selectFile)), 'photoDiodeMUA','-v7.3');
    end
    all_photoMUA.(strcat('x',selectFile)) =  photoDiodeMUA;
end
end