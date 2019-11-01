%% combine array data's that exist within a folder structure

    single_pulse = 0;
    if(single_pulse)
        highest_folderpath{1} = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\StimRecData\Duncan\SinglePulse\';
        highest_folderpath{2} = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\StimRecData\Han\SinglePulse\';
        search_word = 'array_data';
    else 
        highest_folderpath{1} = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\StimRecData\Duncan\DblPulse_trains\';
        highest_folderpath{2} = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\StimRecData\Han\DblPulse_trains\all_processed\';
        search_word = 'arrayData';
    end
    
    monkey_name{1} = 'Duncan';
    monkey_name{2} = 'Han';

    array_data_all = {};
    input_data_all = {};
    
    pwd = cd;
    for h = 1:numel(highest_folderpath)
        cd(highest_folderpath{h})
        files = dir([highest_folderpath{h} '**\*',search_word,'*']);

        for f = 1:numel(files)
            % load files(f)
            load([files(f).folder,'\' files(f).name]);
            
            
            % get stim chan from filename
            chan_idx = strfind(files(f).name,'chan');
            underscore_idx = strfind(files(f).name,'_');
            underscore_idx = underscore_idx(find(underscore_idx > chan_idx,1,'first'));
            
            stim_chan = num2str(files(f).name(chan_idx+4:underscore_idx-1));
            % put arrayData into array_data_all;
            if(single_pulse)
                num_units = numel(arrayData);
            else
                num_units = numel(array_data);
            end
            
            for arr_idx = 1:num_units
                if(single_pulse)
                    array_data_all{end+1} = arrayData{arr_idx};
                else
                    array_data_all{end+1} = array_data{arr_idx};
                    array_data_all{end}.stimData = array_data_all{end}.trial_num;
                    array_data_all{end}.numStims = array_data_all{end}.num_stims;
                    input_data_all{end+1} = input_data;
                end                
                
                array_data_all{end}.CHAN_LIST = stim_chan;
                
                array_data_all{end}.monkey = monkey_name{h};
            end
        end
    end

    
    if(single_pulse)
        arrayData = array_data_all; % rename to match other scripts
    else
        array_data = array_data_all; % match with train script 
    end
   