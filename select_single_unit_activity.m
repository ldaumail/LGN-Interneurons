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

%% contrast levels encountered in the dominant eye

 allContLevels =[];
 for i =1:71
     if ~isempty(filenames{i})
         contLevels = unique(unitsData.new_data(i).channel_data.fixedc);
         allContLevels = unique([allContLevels; contLevels]);
     end
 end
 
 
%% SECOND: 

allfilename = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\all_units\all_locs_data_95CI_05022021';
save(strcat(allfilename, '.mat'), 'peakLocs');


%allfilename = [newdatadir 'su_peaks_03032020_corrected\all_units\all_data_peaks'];
 %save(strcat(allfilename, '.mat'), 'data_peaks');
 allfilename = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\all_units\FiltMultiContSUA_05022021';
 save(strcat(allfilename, '.mat'), 'FiltMultiContSUA');
% allfilename = [newdatadir 'su_peaks_03032020_corrected\all_units\clean_SUA_locs'];
% save(strcat(allfilename, '.mat'), 'peaks_locs');
 allfilename = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\all_units\NoFiltMultiContSUA_06212021';
 save(strcat(allfilename, '.mat'), 'NoFiltMultiContSUA');
 
 allfilename = 'C:\Users\daumail\OneDrive - Vanderbilt\Documents\LGN_data_042021\single_units\binocular_adaptation\all_units\BsNoFiltMultiContSUA_05022021';
 save(strcat(allfilename, '.mat'), 'BsNoFiltMultiContSUA');
 
