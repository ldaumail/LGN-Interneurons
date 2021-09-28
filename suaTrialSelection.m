function [peakLocs, NoFiltMultiContSUA] = suaTrialSelection(channelsdir, filenames)
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
contLims = [0,0.1,0.3,0.5,0.7,1];  
contPairs = combvec(contLims, contLims); %all possible combinations of contrast limit levels
channum = 1: length(unitsData.new_data);
xabs = -199:1300;
nyq = 500;

%mean_filtered_dSUA = struct();


 FiltMultiContSUA =  struct();
 NoFiltMultiContSUA = struct();
 BsNoFiltMultiContSUA = struct();
%data_peaks = struct();
peakLocs = struct(); %store filtered data peak locations used to isolate peak values of unfiltered data

for n = 1:size(contPairs,2)
    clear i
     for i = channum
        if ~isempty(filenames{i})
            filename = filenames(i);
            
            NDETrCt = unitsData.new_data(i).channel_data.fixedc; %trial by trial contrasts in NDE
            DETrCt = unitsData.new_data(i).channel_data.contrast; %trial by trial contrasts in DE

            blankcontrast = DETrCt ==  0 & NDETrCt ==  0; %get logicacal indices of trials with 0 contrast in both eyes
                        
            lowNDECont = contPairs(1,n); %lower NDE contrast boundary
            highNDECont = contLims(find(contLims == contPairs(1,n))+1); %higher NDE contrast boundary
            
            lowDECont = contPairs(2,n); %lower DE contrast boundary
            highDECont = contLims(find(contLims == contPairs(2,n))+1); %higher DE contrast boundary
  
            contrastBin = (NDETrCt >  lowNDECont &  NDETrCt <= highNDECont)& (DETrCt > lowDECont & DETrCt <= highDECont);
            
            trialidx = 1:length(unitsData.new_data(i).channel_data.sdftr_chan(1,:)); %trial number of each trial for a given unit
            noFiltBs = nan(length(xabs), length(trialidx)); %to store the baseline corrected unfiltered data
            filtBs = nan(length(xabs), length(trialidx));
            origin_data = nan(length(xabs)+401, length(trialidx));
            %all_norm_lpdSUA= nan(length(xabs),length(trialidx));
            
            powerstim = nan(length(trialidx),1025);
            freqstim = nan(length(trialidx),1025);
            fourhzpowerstim =nan(length(trialidx),1);
            % bsl = nan(1, length(trialidx));
            mean_wnd1 = nan(1,length(trialidx));
            
            all_pks = nan(4,length(unitsData.new_data(i).channel_data.sdftr_chan(1,contrastBin)));

         for tridx = trialidx

                all_data = unitsData.new_data(i).channel_data.sdftr_chan(401:1900,tridx);
                origin_data(:,tridx) = unitsData.new_data(i).channel_data.sdftr_chan(:,tridx);
                noFiltBs(:,tridx) = all_data(1:end)- mean(all_data(1:200));
                
               
                lpc       = 4.5; %low pass cutoff
                lWn       = lpc/nyq;
                [bwb,bwa] = butter(4,lWn,'low');
                lpdSUA      = filtfilt(bwb,bwa, noFiltBs(:,tridx));


                filtBs(:,tridx) = lpdSUA;
                %all_norm_lpdSUA(:,tridx) = (lpdSUA - min(lpdSUA))/(max(lpdSUA)- min(lpdSUA));
                mean_wnd1(tridx) = mean(lpdSUA(201:480)); %compute mean spiking response over 280ms following stimulus onset. 250+30

                %%% power


                [powerstim(tridx,:), freqstim(tridx,:)] = calcFFT(all_data(200:1350)); %fourrier transform

                %find the index of the frequency vector closest to 4hz and point to the
                %power value of this index for every trial, and store the value in
                %fourhzpower
                [val,index] = min(abs(4-freqstim(tridx,:))); %find index closest to 4Hz
                fourhzpowerstim(tridx,1) = powerstim(tridx,index); %get power at that index, assumed to be 4Hz

         end

      %%%%%%%%%%% %reject trials below Mean + 1.96*STD in the blank condition %%%%%%
       %power related variables
       power0 = fourhzpowerstim(blankcontrast); %power of responses in blank condition
       powerDE = fourhzpowerstim(contrastBin); %power of responses with contrast stimulus >0 in DE and 0 contrast in NDE

       %spiking activity related variables
       mean_wnd1_DE =mean_wnd1(contrastBin);
       filtered_dSUA_high = filtBs(:,contrastBin);
       %filtered_dSUA_blank = filtBs(:,blankcontrast);
       origin_data_high = origin_data(:,contrastBin);
       %origin_data_blank = origin_data(:,blankcontrast);
       bsl_origin_data_high = noFiltBs(:,contrastBin);
       trials = unitsData.new_data(i).channel_data.trial(contrastBin);
  
       sua_bsl =  mean(filtered_dSUA_high(1:200,:),1);
       
 

        %%%%%%%%%%%determine the first peak location for each trial of a given single unit %%%%%%%%
        all_locsdSUA_trials = nan(6,length(filtered_dSUA_high(1,:)));
        clear trial
        for trial = 1:length(filtered_dSUA_high(1,:))
            
            for ln = 1:550
                if filtered_dSUA_high(200+ln,trial) < filtered_dSUA_high(200+ln+1,trial) && ~all(isnan(filtered_dSUA_high(:,trial)))
                    [~,locsdSUA_trial] = findpeaks(filtered_dSUA_high(200+ln:1499,trial));
                     
                    %if peak1 is too small, peak2 becomes peak1
                    if filtered_dSUA_high(locsdSUA_trial(1)+200+ln,trial) >= 0.4*filtered_dSUA_high(locsdSUA_trial(2)+200+ln)
                        %store first peak location
                        all_locsdSUA_trials(1:length(locsdSUA_trial),trial) = locsdSUA_trial(1:end)+200+ln;
                    else
                        all_locsdSUA_trials(1:length(locsdSUA_trial(2:end)),trial) = locsdSUA_trial(2:end)+200+ln;
                        
                    end
                    
                    break
                end
            end
            
            if nnz(~isnan(all_locsdSUA_trials(:,trial))) >= 4 && ~all(isnan(all_locsdSUA_trials(:,trial)))
                %adjust location to the first data point of lpsu (+ln),
                
                all_pks(:,trial) = filtered_dSUA_high(all_locsdSUA_trials(1:4,trial), trial);
                filtered_dSUA_high(:,trial) = filtered_dSUA_high(:,trial);
                all_locsdSUA_trials(:,trial) = all_locsdSUA_trials(:,trial);
                origin_data_high(:,trial) = origin_data_high(:,trial);
                bsl_origin_data_high(:,trial) = bsl_origin_data_high(:,trial);
            else
                filtered_dSUA_high(:,trial) = nan(length(filtered_dSUA_high(:,trial)),1);
                all_locsdSUA_trials(:,trial) = nan(size(all_locsdSUA_trials(:,trial)));
                origin_data_high(:,trial) =  nan(length(origin_data_high(:,trial)),1);
                bsl_origin_data_high(:,trial) =  nan(length(bsl_origin_data_high(:,trial)),1);
                
            end
            
            if ~all(isnan(all_locsdSUA_trials(:,trial))) && (all_locsdSUA_trials(4,trial) ~= 1500) %remove trials for which the pk4 is the last data point (= not a peak if this happens)
                %adjust location to the first data point of lpsu (+ln),
                
                all_pks(:,trial) = filtered_dSUA_high(all_locsdSUA_trials(1:4,trial), trial);
                filtered_dSUA_high(:,trial) = filtered_dSUA_high(:,trial);
                all_locsdSUA_trials(:,trial) = all_locsdSUA_trials(:,trial);
                origin_data_high(:,trial) = origin_data_high(:,trial);
                bsl_origin_data_high(:,trial) = bsl_origin_data_high(:,trial);
            else
                all_pks(:,trial) = nan(length(all_pks(:,trial)),1);
                filtered_dSUA_high(:,trial) = nan(length(filtered_dSUA_high(:,trial)),1);
                all_locsdSUA_trials(:,trial) = nan(size(all_locsdSUA_trials(:,trial)));
                origin_data_high(:,trial) =  nan(length(origin_data_high(:,trial)),1);
                bsl_origin_data_high(:,trial) =  nan(length(bsl_origin_data_high(:,trial)),1);
                
                
            end
        end
        %{
        figure(); plot(-199:1300, filtered_dSUA_high(1:1500,:))
        hold on
        plot(all_locsdSUA_trials(1:4,1)-200, all_pks(:,1))
        set(gca,'box','off')
        %}
        %%% reject outlier peaks and the corresponding trials in
        %%% filtered_dSUA_high
        
        
        %%%%%%%%%%reject if there is a peak 1 outlier, if the max peak value in the baseline is an outlier %%%%%%%%%
        
        % First find peaks before stimulus onset

        bsl_peaks = nan(1, length(filtered_dSUA_high(1,:)));
        clear tr
        for tr = 1:length(filtered_dSUA_high(1,:))
            
            for loc = 1:200
                if filtered_dSUA_high(loc,tr) < filtered_dSUA_high(loc+1,tr) && ~all(isnan(filtered_dSUA_high(:,tr)))
                    if length(filtered_dSUA_high(loc:200,tr)) >= 3
                        if ~isempty(findpeaks(filtered_dSUA_high(loc:200,tr)))
                            [~, bsl_peak_locs] = findpeaks(filtered_dSUA_high(loc:200,tr));
                            bsl_peaks(1,tr) = max(filtered_dSUA_high(bsl_peak_locs+loc,tr));
                        else
                            bsl_peaks(1,tr) = NaN;
                        end
                    end
                    break
                end
            end
        end

        out_bsl_peaks = isoutlier(bsl_peaks);

        p1outliers = isoutlier(all_pks(1,:));
        clear tr
        for tr = 1:length(filtered_dSUA_high(1,:))
            %exclude trials
            if p1outliers(tr) == 0 && ~all(isnan(all_pks(:,tr))) && out_bsl_peaks(tr) ==0 

                filtered_dSUA_high(:,tr) = filtered_dSUA_high(:, tr);
                all_pks(:, tr) = all_pks(:,tr);
                all_locsdSUA_trials(:,tr) = all_locsdSUA_trials(:,tr);
                origin_data_high(:,tr) = origin_data_high(:, tr);
                bsl_origin_data_high(:,tr) =  bsl_origin_data_high(:,tr);


            else 
                filtered_dSUA_high(:,tr) = nan(length(filtered_dSUA_high(:,tr)),1);
                all_pks(:,tr) = nan(length(all_pks(:,tr)),1);
                all_locsdSUA_trials(:,tr) = nan(size(all_locsdSUA_trials(:,tr)));
                origin_data_high(:,tr) = nan(length(origin_data_high(:,tr)),1);
                bsl_origin_data_high(:,tr) =  nan(length(bsl_origin_data_high(:,tr)),1);

            end
        end
       all_locsdSUA_trials =  all_locsdSUA_trials(:,~all(isnan(all_locsdSUA_trials))); % remove nan - cols
       all_pks = all_pks(:, ~all(isnan(all_pks)));
       selTrials = trials(~all(isnan(origin_data_high)));
      
       filtered_dSUA_high = filtered_dSUA_high(:,~all(isnan(filtered_dSUA_high))); 
       origin_data_high = origin_data_high(:,~all(isnan(origin_data_high)));
       bsl_origin_data_high = bsl_origin_data_high(:,~all(isnan(bsl_origin_data_high)));
        
       filename = sprintf('x%s',char(filename));
       binNb = sprintf('bin%d', n);
       if n ==1 && length(filtered_dSUA_high(1,:)) >=10 %first bin == high contrast monocular condition will serve as an indicator of the minimum number of trials required for the analysis
           
           %eval(['peakLocs.' num2str(i) '.bin' num2str(n) ' = all_locsdSUA_trials;'])
          
           peakLocs.(filename).(binNb) = all_locsdSUA_trials; %create dynamical peak locations structures
           FiltMultiContSUA.(filename).(binNb) =  filtered_dSUA_high;
           NoFiltMultiContSUA.(filename).(binNb).neuralDat = origin_data_high;
           BsNoFiltMultiContSUA.(filename).(binNb) = bsl_origin_data_high;
           
           NoFiltMultiContSUA.(filename).(binNb). trials = selTrials; %save trial indices
           NoFiltMultiContSUA.(filename).(binNb). peaklocs =all_locsdSUA_trials;
           
           FiltMultiContSUA.(filename).cellclass = cellClass{i}; %get the cell class of selected units
           NoFiltMultiContSUA.(filename).cellclass = cellClass{i};
           BsNoFiltMultiContSUA.(filename).cellclass = cellClass{i};
           
       elseif n == 1 && length(filtered_dSUA_high(1,:)) <10
           
           all_pks(:,:) = [];
           FiltMultiContSUA.(filename).(binNb) =  [];
           NoFiltMultiContSUA.(filename).(binNb) = [];
           peakLocs.(filename).(binNb) = [];
        
       elseif n > 1 && length(filtered_dSUA_high(1,:)) >=10 
           peakLocs.(filename).(binNb) = all_locsdSUA_trials; %create dynamical peak locations structures
           FiltMultiContSUA.(filename).(binNb) =  filtered_dSUA_high;
           % FiltMultiContSUA.(filename).bin0 =  filtered_dSUA_blank;
           NoFiltMultiContSUA.(filename).(binNb).neuralDat = origin_data_high;
           BsNoFiltMultiContSUA.(filename).(binNb) = bsl_origin_data_high;
           
           NoFiltMultiContSUA.(filename).(binNb). trials = selTrials;
           NoFiltMultiContSUA.(filename).(binNb). peaklocs =all_locsdSUA_trials;
           %NoFiltMultiContSUA.(filename).bin0 = origin_data_blank;
           FiltMultiContSUA.(filename).cellclass = cellClass{i}; %get the cell class of selected units
           NoFiltMultiContSUA.(filename).cellclass = cellClass{i};
           BsNoFiltMultiContSUA.(filename).cellclass = cellClass{i};
       end
       

     %data_peaks(i).namelist = all_pks(:,~all(isnan(all_pks)));
     %all_pks = all_pks(:,~all(isnan(all_pks)));
    %channelfilename = [unitsDir 'su_peaks_03032020_corrected\individual_units\' filename 'multiContrast'];
    %save(strcat(channelfilename, '.mat'), 'peakLocs');
        end
     end
end

