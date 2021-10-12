function [all_sigs95] = binocularModulationPowerTest(suaData, binoConds, monoConds)
%function that performs tests of comparisons between monocular and multiple
%written on 09-28-2021

%binocular conditions
%suaData = selected_data; 
xfilenames = fieldnames(suaData);

Ses = struct();
%bs_data = struct();
channum = 1: length(fieldnames(suaData));
%mean_S = nan(1174,38, length(channum));

%xabs = -100:1301;

%filtered_dMUA = nan(length(xabs), length(channum));
%dim 2 = channel, dim3 = trials
Fs = 1000;
movingwin       = [.256 .001]; % length of moving window in seconds (should be to the power of 2) + length of sliding window
params.tapers   = [2 3];
params.Fs       = Fs;
params.fpass    = [1 150];

all_sigs95 = nan(length(fieldnames(suaData)),length(monoConds));
all_sigs90 = nan(length(fieldnames(suaData)),length(monoConds));

for i = 1:length(channum)
    xfilename = xfilenames{i};
    for n =1:length(monoConds) %could be 1:size(binoConds,2) as well
        monoBinNb = sprintf('lowcontDE%dNDE0Mono', monoConds(n));
        binoBinNb = sprintf('lowcontDE%dNDE%d',binoConds(1,n), binoConds(2,n));

        if ~isempty(suaData.(xfilename).NoFiltMultiContSUA.(binoBinNb)) & ~isempty(suaData.(xfilename).NoFiltMultiContSUA.(monoBinNb))
            
            binoData = suaData.(xfilename).NoFiltMultiContSUA.(binoBinNb)(401:1900,:);
            binoBsl = mean(binoData(1:200,:));
            bino_mean_bs = binoData(72:end, :) - binoBsl;
            
            monoData = suaData.(xfilename).NoFiltMultiContSUA.(monoBinNb)(401:1900,:);
            monoBsl = mean(monoData(1:200,:));
            mono_mean_bs = monoData(72:end, :) - monoBsl;
            
            %{
            figure();
            plot(mean(bino_mean_bs,2))
            hold on
            plot(mean(mono_mean_bs,2))
            legend('Bino','Mono')
            %}
            
            %bs_data(i).(xfilename).(binoBinNb) = bino_mean_bs;
            clear binoS monoS;
            [binoS,t,f]  = mtspecgramc(bino_mean_bs(:,:) ,movingwin, params);
            Ses(i).(xfilename).(binoBinNb) = binoS;
            %mean_S(:,:,i) = nanmean(S,3);
            [monoS,t,f]  = mtspecgramc(mono_mean_bs(:,:) ,movingwin, params);
            Ses(i).(xfilename).(monoBinNb) = monoS;
            tvec = t*1000 -129;
    
    %time_adj = 1:128;
    %tvec = cat(2, time_adj , t*1000) ;
    %we can also store tvec and f in a struct, but they are all identical
        
 
        part1 = nanmean(squeeze(Ses(i).(xfilename).(binoBinNb)(:,1,:)), 1);
        part2 = nanmean(squeeze(Ses(i).(xfilename).(monoBinNb)(:,1,:)), 1);
        
        if nanmean(part1) > nanmean(part2)
            cond1               = part1;
            cond2               = part2;
        else
            cond1               = part2;
            cond2               = part1;
        end
        
        [X,Y,T,AUC]      = perfcurve([ones(length(cond1),1); repmat(2,length(cond2),1)],[cond1 cond2],1); %returns X and Y coordinates of ROC curve on our mono vs bino samples + AUC, given the labels (cond1 cond2) and values associated (or scores) 
        NP                    = length(cond2);
        PR                    = length(cond1);
        catdat               = [cond1 cond2];
        reps   = 10000;
        for r = 1:reps
            clear shufNP shufPR
            shufPR         = catdat(randperm(length(catdat),PR));
            shufNP         = catdat(randperm(length(catdat),NP));
            [~,~,~,...
                shufAUC(r)]    = perfcurve([ones(PR,1); repmat(2,NP,1)],[shufPR shufNP],1);
        end
        critT95         = quantile(shufAUC,.95);
        %critT90         = quantile(shufAUC,.90);
        
        if (AUC > critT95) && (nanmean(part1) > nanmean(part2))%if AUC superior to 95% quantile, likely not by chance, reject null hypothesis
            sig95          = 1;
        else
            if (AUC > critT95) && (nanmean(part1) < nanmean(part2))
                sig95          = -1;
            else
                sig95          = 0;
            end
        end
        
        all_sigs95(i,n) = sig95;
        %{
        if AUC > critT90
            sig90        = 1;
        else
            sig90         = 0;
        end
        all_sigs90(i,n) = sig90;
        %}
        end
    end
end