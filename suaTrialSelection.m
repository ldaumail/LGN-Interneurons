function [selectData] = suaTrialSelection(unitsData, filenames)
%This function allows to select:
%    -the peak locations of the smoothed (low-pass filtered) data
%    - as well as the trials data they correspond to.
%    - the binary data of those trials
%Last edited by Loic Daumail on 09/27/2021
%Location of multiunits data isolated in the previous section

%select trials and isolate peak locations of smoothed data ==> use code from %%"get_clean_peaks_and_data.m",
% use multiple contrast levels
% select trials the same way as in the monocular adaptation
% analysis
 %% find peak locations of smoothed data, to further allow us to isolate peak values of unfiltered data in order to analyze them on R and fit a LMER
 %find the multiple contrast levels present in the data
 
 %lets create contrast limits (bins to pool different contrast levels)
 
 % option 1: more stringent = more bins
 contLims = [0,0.1,0.3,0.5,0.7,1]; 
 contPairs = combvec(contLims(1:end-1), contLims(1:end-1)); %all possible combinations of lower contrast boundary levels (exclude highest contrast, as it is the highest boundary, it won't work in the for loop)
 
 channum = 1: length(unitsData.new_data);
 xabs = -199:1300;
 nyq = 500;
 selectData = struct(); %store: 1) Convolved trial data, 2)Binary spike data 3) filtered data peak locations used to isolate peak values of unfiltered data 4) channel number
 
 for n = 1:size(contPairs,2)+5
     clear i
     for i = channum
         if ~isempty(filenames{i})
             filename = filenames{i};
             
             NDETrCt = unitsData.new_data(i).channel_data.fixedc; %trial by trial contrasts in NDE
             DETrCt = unitsData.new_data(i).channel_data.contrast; %trial by trial contrasts in DE
             
             blankcontrast = DETrCt ==  0 & NDETrCt ==  0; %get logicacal indices of trials with 0 contrast in both eyes
             
             if n < size(contPairs,2)+1
             %if (find(contLims == contPairs(1,n))+1 <= length(contLims)) & (find(contLims == contPairs(2,n))+1 <= length(contLims))
             lowNDECont = contPairs(1,n); %lower NDE contrast boundary
             highNDECont = contLims(find(contLims == contPairs(1,n))+1); %higher NDE contrast boundary
             
             lowDECont = contPairs(2,n); %lower DE contrast boundary
             highDECont = contLims(find(contLims == contPairs(2,n))+1); %higher DE contrast boundary
             
             contrastBin = (NDETrCt >  lowNDECont &  NDETrCt <= highNDECont)& (DETrCt > lowDECont & DETrCt <= highDECont);
              else
                  if n >= size(contPairs,2)+1 %for the monocular condition
                      lowDECont = contLims(n-size(contPairs,2)); %lower DE contrast boundary
                      highDECont = contLims(n - size(contPairs,2)+1); %higher DE contrast boundary
                      contrastBin =  (NDETrCt ==  0) & (DETrCt > lowDECont & DETrCt <= highDECont);
                  end
             end
             trialidx = 1:length(unitsData.new_data(i).channel_data.sdftr_chan(1,:)); %trial number of each trial for a given unit
             noFiltBs = nan(length(xabs), length(trialidx)); %to store the baseline corrected unfiltered data
             filtBs = nan(length(xabs), length(trialidx));
             origin_data = nan(length(xabs)+401, length(trialidx)); %this is the convolved SUA we will eventually keep
             
             powerstim = nan(length(trialidx),1025);
             freqstim = nan(length(trialidx),1025);
             fourhzpowerstim =nan(length(trialidx),1);
             mean_wnd1 = nan(1,length(trialidx));
             
             all_pks = nan(4,length(unitsData.new_data(i).channel_data.sdftr_chan(1,contrastBin)));
             
             for tridx = trialidx
                 
                 origin_data(:,tridx) = unitsData.new_data(i).channel_data.sdftr_chan(:,tridx);
                 all_data = unitsData.new_data(i).channel_data.sdftr_chan(401:1900,tridx);
                 noFiltBs(:,tridx) = all_data(1:end)- mean(all_data(1:200));
                 
                 lpc       = 4.5; %low pass cutoff
                 lWn       = lpc/nyq;
                 [bwb,bwa] = butter(4,lWn,'low');
                 lpdSUA      = filtfilt(bwb,bwa, noFiltBs(:,tridx));
                 
                 filtBs(:,tridx) = lpdSUA;
                 mean_wnd1(tridx) = mean(lpdSUA(201:480)); %compute mean spiking response over 280ms following stimulus onset. 250+30
                 
                 %%% power
                 [powerstim(tridx,:), freqstim(tridx,:)] = calcFFT(all_data(200:1350)); %fourrier transform
                 %find the index of the frequency vector closest to 4hz and point to the
                 %power value of this index for every trial, and store the value in
                 %fourhzpower
                 [val,index] = min(abs(4-freqstim(tridx,:))); %find index closest to 4Hz
                 fourhzpowerstim(tridx,1) = powerstim(tridx,index); %get power at that index, assumed to be 4Hz
                 
             end
             
             %%%%%%%%%%% %keep trials with response within 280ms after stim onset > baseline && reject trials below Mean + 1.96*STD in the blank condition %%%%%%
             %power related variables
             power0 = fourhzpowerstim(blankcontrast); %power of responses in blank condition
             powerResp = fourhzpowerstim(contrastBin); %power of responses with contrast stimulus >0 in DE and 0 contrast in NDE
             
             %spiking activity related variables
             mean_wnd1_Resp =mean_wnd1(contrastBin);
             origin_data_Resp = origin_data(:,contrastBin);
             bsl_origin_data_Resp = noFiltBs(:,contrastBin);
             filtered_dSUA_Resp = filtBs(:,contrastBin);
             sua_bsl =  mean(filtered_dSUA_Resp(1:200,:),1);
             trialnb = find(contrastBin);
             %trials = unitsData.new_data(i).channel_data.trial(contrastBin);
             
             clear tr
             for tr = 1:length(powerResp)
                 if mean_wnd1_Resp(tr) > mean(sua_bsl)+1.96*std(sua_bsl)  && powerResp(tr) > mean(power0)+1.96*std(power0) %/sqrt(length(sua_bsl)) /sqrt(length(power0))
                     
                     filtered_dSUA_Resp(:,tr) = filtered_dSUA_Resp(:,tr);
                     origin_data_Resp(:,tr) = origin_data_Resp(:,tr);
                     bsl_origin_data_Resp(:,tr) = bsl_origin_data_Resp(:,tr);
                 else
                     
                     filtered_dSUA_Resp(:,tr) = nan(length(filtered_dSUA_Resp(:,tr)),1);
                     origin_data_Resp(:,tr) =  nan(length(origin_data_Resp(:,tr)),1);
                     bsl_origin_data_Resp(:,tr) = nan(length(bsl_origin_data_Resp(:,tr)),1);
                 end
             end
             
             %%%%%%%%%%%determine the first peak location for each trial of a given single unit %%%%%%%%
             all_locsdSUA_trials = nan(6,length(filtered_dSUA_Resp(1,:)));
             clear trial
             for trial = 1:length(filtered_dSUA_Resp(1,:))  
                 for ln = 1:550
                     if filtered_dSUA_Resp(200+ln,trial) < filtered_dSUA_Resp(200+ln+1,trial) && ~all(isnan(filtered_dSUA_Resp(:,trial)))
                         locsdSUA_trial = findpeaks(filtered_dSUA_Resp(200+ln:1499,trial));
                         
                         %if peak1 is too small, peak2 becomes peak1
                         if filtered_dSUA_Resp(locsdSUA_trial.loc(1)+200+ln,trial) >= 0.4*filtered_dSUA_Resp(locsdSUA_trial.loc(2)+200+ln)
                             %store first peak location
                             all_locsdSUA_trials(1:length(locsdSUA_trial.loc),trial) = locsdSUA_trial.loc(1:end)+200+ln;
                         else
                             all_locsdSUA_trials(1:length(locsdSUA_trial.loc(2:end)),trial) = locsdSUA_trial.loc(2:end)+200+ln;
                             
                         end
                         
                         break
                     end
                 end
                 
                 if nnz(~isnan(all_locsdSUA_trials(:,trial))) >= 4 && ~all(isnan(all_locsdSUA_trials(:,trial)))
                     %adjust location to the first data point of lpsu (+ln),
                     
                     all_pks(:,trial) = filtered_dSUA_Resp(all_locsdSUA_trials(1:4,trial), trial);
                     filtered_dSUA_Resp(:,trial) = filtered_dSUA_Resp(:,trial);
                     all_locsdSUA_trials(:,trial) = all_locsdSUA_trials(:,trial);
                     origin_data_Resp(:,trial) = origin_data_Resp(:,trial);
                     bsl_origin_data_Resp(:,trial) = bsl_origin_data_Resp(:,trial);
                 else
                     filtered_dSUA_Resp(:,trial) = nan(length(filtered_dSUA_Resp(:,trial)),1);
                     all_locsdSUA_trials(:,trial) = nan(size(all_locsdSUA_trials(:,trial)));
                     origin_data_Resp(:,trial) =  nan(length(origin_data_Resp(:,trial)),1);
                     bsl_origin_data_Resp(:,trial) =  nan(length(bsl_origin_data_Resp(:,trial)),1);
                     
                 end
                 
                 if ~all(isnan(all_locsdSUA_trials(:,trial))) && (all_locsdSUA_trials(4,trial) ~= 1500) %remove trials for which the pk4 is the last data point (= not a peak if this happens)
                     %adjust location to the first data point of lpsu (+ln),
                     
                     all_pks(:,trial) = filtered_dSUA_Resp(all_locsdSUA_trials(1:4,trial), trial);
                     filtered_dSUA_Resp(:,trial) = filtered_dSUA_Resp(:,trial);
                     all_locsdSUA_trials(:,trial) = all_locsdSUA_trials(:,trial);
                     origin_data_Resp(:,trial) = origin_data_Resp(:,trial);
                     bsl_origin_data_Resp(:,trial) = bsl_origin_data_Resp(:,trial);
                 else
                     all_pks(:,trial) = nan(length(all_pks(:,trial)),1);
                     filtered_dSUA_Resp(:,trial) = nan(length(filtered_dSUA_Resp(:,trial)),1);
                     all_locsdSUA_trials(:,trial) = nan(size(all_locsdSUA_trials(:,trial)));
                     origin_data_Resp(:,trial) =  nan(length(origin_data_Resp(:,trial)),1);
                     bsl_origin_data_Resp(:,trial) =  nan(length(bsl_origin_data_Resp(:,trial)),1);
                     
                     
                 end
             end
             %{
        figure(); plot(-199:1300, filtered_dSUA_Resp(1:1500,:))
        hold on
        plot(all_locsdSUA_trials(1:4,1)-200, all_pks(:,1))
        set(gca,'box','off')
             %}
             %%% reject outlier peaks and the corresponding trials in
             %%% filtered_dSUA_high
             
             
             %%%%%%%%%%reject if there is a peak 1 outlier, if the max peak value in the baseline is an outlier %%%%%%%%%
             
             % First find peaks before stimulus onset
             
             bsl_peaks = nan(1, length(filtered_dSUA_Resp(1,:)));
             clear tr
             for tr = 1:length(filtered_dSUA_Resp(1,:))
                 
                 for loc = 1:200
                     if filtered_dSUA_Resp(loc,tr) < filtered_dSUA_Resp(loc+1,tr) && ~all(isnan(filtered_dSUA_Resp(:,tr)))
                         if length(filtered_dSUA_Resp(loc:200,tr)) >= 3
                             if ~isempty(findpeaks(filtered_dSUA_Resp(loc:200,tr)))
                                 bsl_peak_locs = findpeaks(filtered_dSUA_Resp(loc:200,tr));
                                 bsl_peaks(1,tr) = max(filtered_dSUA_Resp(bsl_peak_locs.loc+loc,tr));
                             else
                                 bsl_peaks(1,tr) = NaN;
                             end
                         end
                         break
                     end
                 end
             end
             
             % bsl peaks above outlier threshold
            % out_bsl_peaks = isoutlier(bsl_peaks);
             out_bsl_peaks = bsl_peaks > nanmean(bsl_peaks)+ 1.96*std(bsl_peaks,[],'omitnan');
             
             %p1outliers = isoutlier(all_pks(1,:));
             p1outliers = all_pks(1,:) > nanmean(all_pks)+ 1.96*std(all_pks,[],'omitnan');
             clear tr
             for tr = 1:length(filtered_dSUA_Resp(1,:))
                 %exclude trials
                 if p1outliers(tr) == 0 && ~all(isnan(all_pks(:,tr))) && out_bsl_peaks(tr) ==0
                     
                     filtered_dSUA_Resp(:,tr) = filtered_dSUA_Resp(:, tr);
                     all_pks(:, tr) = all_pks(:,tr);
                     all_locsdSUA_trials(:,tr) = all_locsdSUA_trials(:,tr);
                     origin_data_Resp(:,tr) = origin_data_Resp(:, tr);
                     bsl_origin_data_Resp(:,tr) =  bsl_origin_data_Resp(:,tr);
                     
                     
                 else
                     filtered_dSUA_Resp(:,tr) = nan(length(filtered_dSUA_Resp(:,tr)),1);
                     all_pks(:,tr) = nan(length(all_pks(:,tr)),1);
                     all_locsdSUA_trials(:,tr) = nan(size(all_locsdSUA_trials(:,tr)));
                     origin_data_Resp(:,tr) = nan(length(origin_data_Resp(:,tr)),1);
                     bsl_origin_data_Resp(:,tr) =  nan(length(bsl_origin_data_Resp(:,tr)),1);
                     
                 end
             end
             trialnb = trialnb(~all(isnan(origin_data_Resp)));
             binary_data_Resp = unitsData.new_data(i).channel_data.spk_bin_chan(:,trialnb); %get the binary spikes data
             origin_data_Resp = origin_data_Resp(:,~all(isnan(origin_data_Resp)));
             all_locsdSUA_trials =  all_locsdSUA_trials(:,~all(isnan(all_locsdSUA_trials)));
             
             %filtered_dMUA_high = filtered_dMUA_high(:,~all(isnan(filtered_dMUA_high))); % for nan - cols
             %all_pks = all_pks(:, ~all(isnan(all_pks)));
             %bsl_origin_data_high = bsl_origin_data_high(:,~all(isnan(bsl_origin_data_high)));
             
             
             xfilename = strcat('x',filename);
             if n < size(contPairs,2)+1
                 binNb = sprintf('lowcontDE%dNDE%d', lowDECont*100, lowNDECont*100);
             else
                 if n >= size(contPairs,2)+1
                     binNb = sprintf('lowcontDE%dNDE0Mono', lowDECont*100);
                 end
             end
             if  length(origin_data_Resp(1,:)) >=10 %first bin == high contrast monocular condition will serve as an indicator of the minimum number of trials required for the analysis
                 
                 selectData.(xfilename).NoFiltMultiContSUA.(binNb) = origin_data_Resp;
                 selectData.(xfilename).peakLocs.(binNb) = all_locsdSUA_trials; %create dynamical peak locations structures
                 selectData.(xfilename).binSpkTrials.(binNb) = binary_data_Resp;
                 selectData.(xfilename).chan = unitsData.new_data(i).channel_data.chan;
                                 
             elseif  length(origin_data_Resp(1,:)) <10
                 selectData.(xfilename).NoFiltMultiContSUA.(binNb) = [];
                 selectData.(xfilename).peakLocs.(binNb) = [];
                 selectData.(xfilename).binSpkTrials.(binNb) = [];
             end
             
             
             %data_peaks(i).namelist = all_pks(:,~all(isnan(all_pks)));
             %all_pks = all_pks(:,~all(isnan(all_pks)));
             %channelfilename = [unitsDir 'su_peaks_03032020_corrected\individual_units\' filename 'multiContrast'];
             %save(strcat(channelfilename, '.mat'), 'peakLocs');
             %end
         end
     end
 end

