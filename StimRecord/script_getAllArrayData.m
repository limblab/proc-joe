%% combine array data's that exist within a folder structure

    highest_folderpath = 'C:\Users\jts3256\Desktop\Han_stim_data\';
    search_word = 'array_data';


    array_data_all = {};

    pwd = cd;
    cd(highest_folderpath)
    files = dir([highest_folderpath '**\*',search_word,'*']);

    for f = 1:numel(files)
        % load files(f)
        load([files(f).folder,'\' files(f).name]);

        % put arrayData into array_data_all;
        array_data_all{end+1} = arrayData{1};


    end


