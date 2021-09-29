%% First: use the list of single units file names that were selected in the adaptation analysis with
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
 
 
%% Second: Get selected trials data

selected_data = suaTrialSelection(unitsData, filenames);%get 
allfilename = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\lgn_interneuron_suppression\all_cont_lgn_data_09282021';
save(strcat(allfilename, '.mat'), 'selected_data');

%% Test for binocular modulation

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

%% Try with more pooling of contrast levels to save trials
selected_data = suaTrialSelectionContPool(unitsData, filenames);%get 

lowContNDE = [50];
lowContDE = [0];

%binoconds:
binoConds = combvec(lowContDE,lowContNDE);
monoConds = binoConds(1,:);

sigResults = binocularModulationPowerTest(selected_data, binoConds, monoConds);






 
