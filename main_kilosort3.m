%%this script runs Kilosort 3
%Loic Daumail 12/05/2021
%% you need to change most of the paths in this block
addpath(genpath('C:\Users\daumail\OneDrive - Vanderbilt\Documents\MATLAB\Kilosort-main')) % path to kilosort folder
addpath('C:\Users\daumail\OneDrive - Vanderbilt\Documents\MATLAB\npy-matlab-master') % for converting to Phy

%modulxFilenames =  {'x160616_I','x160623_I','x160628_I', ...
%    'x180827_I','x181212_B', 'x190119_B' ,'x190124_B' ,'x190205_B'};
%nchan = [24    24    24    32    24    24    24    24];
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
nchan =  [24    24    24    24    24    24    24    24    24    24    24    24    24    24 ...
    24    24    24    24    24    24    24    24    24    24    24    24     3     3 ...
     3     3     3     3    32    32    32    24    24    24    24    24    24    24 ...
    24    24    24    24];
selectNchan = nchan(SelectIdx);
for i = 25:length(selectSess)
    selectDate = selectSess{i};
    session = selectDate(1:end);%(2:end);
    
    dataDir = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\kilosort_binary\all_sessions';
    rootZ = strcat(dataDir, '\',sprintf('%s',session)); % the raw data binary file is in this folder
    rootH = strcat(dataDir,'\', sprintf('%s',session),'\','temp'); % path to temporary binary file (same size as data, should be on fast SSD)
    mkdir(rootH)
    pathToYourConfigFile = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\loic_code_042021\LGN-Interneurons\kilosort_code\configFiles'; % take from Github folder and put it somewhere else (together with the master_file)
    
    %  create a channel map file  
    Nchannels = selectNchan(i);
    connected = true(Nchannels, 1);
    chanMap   = 1:Nchannels;
    chanMap0ind = chanMap - 1;
    xcoords   = ones(Nchannels,1);
    ycoords   = [1:100:Nchannels*100]'; %version 3
    kcoords   = ones(Nchannels,1); % grouping of channels (i.e. tetrode groups)
    fs = 30000; % sampling frequency
    pathToChanMap = fullfile(strcat(rootZ,'\','chanMapFiles'));
    mkdir(strcat(rootZ,'\','chanMapFiles'));
    save(strcat(pathToChanMap,'\','chanMap.mat'), ...
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')
    
    chanMapFile = 'chanMap.mat';
    
    ops.trange    = [0 Inf]; % time range to sort
    ops.NchanTOT  = Nchannels; % total number of channels in your recording
    
    run(fullfile(pathToYourConfigFile, 'StandardConfig_MOVEME.m'))
    %create temporary data file 
    fid = fopen(strcat(rootH,'\','temp_wh.dat'), 'w');
    fwrite(fid, 0, 'int16');
    fclose(fid);
    ops.fproc   = fullfile(rootH, 'temp_wh.dat'); % proc file on a fast SSD
    ops.chanMap = fullfile(pathToChanMap, chanMapFile);
    %% this block runs all the steps of the algorithm
    fprintf('Looking for data inside %s \n', rootZ)
    
    % main parameter changes from Kilosort2 to v2.5
    ops.sig        = 20;  % spatial smoothness constant for registration
    ops.fshigh     = 300; % high-pass more aggresively
    ops.nblocks    = 5; % blocks for registration. 0 turns it off, 1 does rigid registration. Replaces "datashift" option.
    
    % main parameter changes from Kilosort2.5 to v3.0
    ops.Th       = [9 9];
    
    % is there a channel map file in this folder?
    fs = dir(fullfile(rootZ, 'chan*.mat'));
    if ~isempty(fs)
        ops.chanMap = fullfile(rootZ, fs(1).name);
    end
    
    % find the binary file
    fs          = [dir(fullfile(rootZ, '*.bin')) dir(fullfile(rootZ, '*.dat'))];
    ops.fbinary = fullfile(rootZ, fs(1).name);
    
    rez                = preprocessDataSub(ops);
    rez                = datashift2(rez, 1);
    
    [rez, st3, tF]     = extract_spikes(rez);
    
    rez                = template_learning(rez, tF, st3);
    
    [rez, st3, tF]     = trackAndSort(rez);
    
    rez                = final_clustering(rez, tF, st3);
    
    rez                = find_merges(rez, 1);
    
    %rootZ = fullfile(rootZ, 'kilosort3');
    %mkdir(rootZ)
    rezToPhy2(rez, rootZ);
    save(strcat(rootZ,'\','rez.mat'), 'rez','-v7.3');
end
%% 
