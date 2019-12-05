%% load in a ns5


    folderpath = 'E:\Data\Joseph\Han_stim_data\Han_20191114_longTrains_dukeGen2\chan4\';
        
    cd(folderpath);
    file_list = dir('*.ns5');
    
    analog_pin_idx = 97;
    sync_idx = 98;
    artifact_data = {};
    pwd = cd;
    window = [-2,8]; % ms
    pulse_width_1 = zeros(numel(file_list),1); % us
    pulse_width_2 = zeros(size(pulse_width_1)); % us
    interphase = 53; % us
    
    window_idx = window*30; % convert to data points

    
    for file_num = 1%:numel(file_list)
        disp(file_list(file_num).name);
        NS5 = openNSx([folderpath,file_list(file_num).name],'uV');
    
        artifact_data{file_num} = NS5.Data(analog_pin_idx,:);
        sync_line_data{file_num} = NS5.Data(sync_idx,:);
%         
%         % get pulse widths
%         pw1_idx = strfind(file_list(file_num).name,'PW1');
%         pw2_idx = strfind(file_list(file_num).name,'PW2');
%         amp1_idx = strfind(file_list(file_num).name,'A1');
%         amp2_idx = strfind(file_list(file_num).name,'A2');
%         underscore_idx = strfind(file_list(file_num).name,'_');
%         
%         pulse_width_1(file_num) = str2num(file_list(file_num).name(pw1_idx+4:underscore_idx(find(underscore_idx > pw1_idx,1,'first'))-1));
%         pulse_width_2(file_num) = str2num(file_list(file_num).name(pw2_idx+4:underscore_idx(find(underscore_idx > pw2_idx,1,'first'))-1));
%         amp_1(file_num) = str2num(file_list(file_num).name(amp1_idx+3:underscore_idx(find(underscore_idx > amp1_idx,1,'first'))-1));
%         amp_2(file_num) = str2num(file_list(file_num).name(amp2_idx+3:underscore_idx(find(underscore_idx > amp2_idx,1,'first'))-1));
    end
    cd(pwd);
% get stim on for each file
    freq_deliver = {};
    for file_num = 1%:numel(file_list)
        stim_on=find(diff(sync_line_data{file_num}-mean(sync_line_data{file_num})>3)>.5);
        unique_IPI = unique(diff(stim_on))/30;
        disp(unique_IPI(unique_IPI < 210));
    end