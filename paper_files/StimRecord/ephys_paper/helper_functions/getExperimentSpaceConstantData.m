function [ array_data_all ] = getExperimentSpaceConstantData( input_data )


% combine array data's that exist within a folder structure

    if(input_data.home_computer == 1)
        highest_folderpath{1} = 'D:\Lab\Data\KarthikModelData\SinglePulseNonStimChan\';
    else
        highest_folderpath{1} = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\StimRecData\SinglePulseNonStimChan\';
    end
    search_word = 'arrayData';
    
    array_data_all = {};
    input_data_all = {};
    chan_rec_list = [];
    
    pwd = cd;
    for folder_num = 1:numel(highest_folderpath)
        cd(highest_folderpath{folder_num})
        files = dir([highest_folderpath{folder_num} '**\*',search_word,'*']);

        for file_num = 1:numel(files)
            % load files(f)
            load([files(file_num).folder,'\' files(file_num).name]);
            
            underscore_idx = strfind(files(file_num).name,'_');

            % put arrayData into array_data_all;
            num_units = numel(arrayData);

            for arrayDataIdx = 1:num_units
                
                monkey_underscore_idx = underscore_idx(1);
                array_data_all{end+1} = arrayData{arrayDataIdx};
                chan_rec_list(end+1) = arrayData{arrayDataIdx}.CHAN_REC;
                array_data_all{end}.monkey = files(file_num).name(1:monkey_underscore_idx-1);
                
            end
        end
    end

   


end

