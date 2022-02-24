%% summarize recovery data from multiple sessions and monkeys in one place


    folderpath = 'D:\Lab\Data\StimPDs\RecoveryDataAll\';

    monkey_name_list = {'Han','Pop'};
    array_list = [0,1]; % 0 = S1, 1 = M1;

    
    cd(folderpath);
    fnames = dir('*recovData*');
    
    % load in each file, get amp and array data from arrayData file name,
    % also get which elecs were stimmed
    
    array = [];
    monkey_list = {};
    chan_stim= {};
    num_chan_stim = [];
    amp = [];
    
    tau_list = [];
    ts_list = [];
    dir_list = [];
    mag_list = [];
    
    for i_file = 1:numel(fnames)
        load([folderpath,fnames(i_file).name]);
        
        % get meta data for each file
        num_entries = numel(recov_tau);
        
        for i_arr = 1:numel(array_data_files)
            underscore_idx = strfind(array_data_files(i_arr).name,'_');
            monk_name = array_data_files(i_arr).name(1:underscore_idx(1)-1);
            monk_idx = find(strcmpi(monkey_name_list,monk_name));
            
            monkey_list{end+1} = monk_name;
            array(end+1) = array_list(monk_idx);
            
            and_idx = strfind(array_data_files(i_arr).name,'and');
            chan_idx = strfind(array_data_files(i_arr).name,'chan');
            stim_idx = strfind(array_data_files(i_arr).name,'stim');
            
            if(~isempty(and_idx)) % two channels
                chan_stim{end+1} = [str2num(array_data_files(i_arr).name(chan_idx(1)+4:and_idx-1)),str2num(array_data_files(i_arr).name(and_idx+3:stim_idx(1)-1))];
                num_chan_stim(end+1) = 2;
            else % one channel
                chan_stim{end+1,1} = str2num(array_data_files(i_arr).name(chan_idx(1)+4:stim_idx(1)-1));
                num_chan_stim(end+1) = 1;
            end
            
            
            uA_idx = strfind(array_data_files(i_arr).name,'uA');
            uscore_idx = find(underscore_idx < uA_idx,1,'last');
            amp(end+1) = str2num(array_data_files(i_arr).name(underscore_idx(uscore_idx)+1:uA_idx-1));
        end
        % store data
        tau_list = [tau_list,recov_tau];
        ts_list = [ts_list,recov_ts];
        
        dir_list = [dir_list;dir_proj];
        mag_list = [mag_list;mag_proj];
    end


%% plot recovery time for single channel stim for each amplitude and array loc
    
    f_recov = figure('Position',[2217 422 409 420]); hold on;
    offset = 1;
    amp_list = unique(amp);

    boxplot_params = [];
    boxplot_params.use_same_color_for_all = 1;
    boxplot_params.box_width = 1.5;
    for i_amp = 1:numel(amp_list)
        x_loc = amp_list(i_amp);
        data_mask = amp == amp_list(i_amp) & num_chan_stim == 1;
        % m1 first (1 == m1)
        array_mask = array == 1;
        boxplot_params.master_color = getColorFromList(1,2);
        boxplot_wrapper(x_loc-offset,ts_list(data_mask & array_mask),boxplot_params);
       
        % S1 second (0 == s1)
        array_mask = array == 0;
        boxplot_params.master_color = getColorFromList(1,3);
        boxplot_wrapper(x_loc+offset,ts_list(data_mask & array_mask),boxplot_params);
    end
    formatForLee(gcf);
    set(gca,'fontsize',14)
    xlabel('Amplitude (\muA)');
    ylabel('Recovery time constant (ms)');
    ylim([-100,600])

    ax=gca;
    ax.XTick = amp_list;


%% plot projection dir and magnitude for both cortices and amplitudes
    amp_list = unique(amp);
    color_list = [getColorFromList(1,4);getColorFromList(1,5)];
    
    idx_plot = 2; % 1,2 is stim, 3 is right after stim
    % projection direction for each amplitude
    for i_arr = 0:1 % 0 == S1, 1 == M1
        f=figure();
        keep_mask = array == i_arr & num_chan_stim == 1;
        thetas = dir_list(keep_mask==1,idx_plot);
        rhos = mag_list(keep_mask ==1,idx_plot);
        
        amps = amp(keep_mask==1);
        thetas = [zeros(numel(thetas),1),thetas];
        rhos = [zeros(numel(rhos),1),ones(numel(rhos),1)];
        
        color_idx = mod(find(amps == amp_list'),2);
        color_idx(color_idx == 0) = 2;
        
        x=polarplot(thetas',rhos');
        for i_line = 1:numel(x)
            x(i_line).Color = color_list(color_idx(i_line),:);
        end
    end

%% get paired conditions (same array, same chan stim, different amp)
    pair_idx = []; pair_arr = [];
    chan_list = unique([chan_stim{num_chan_stim==1}]);
    
    for i_arr = unique(array)
        for i_chan = 1:numel(chan_list)
            % check for 25uA and 50uA idx
            chan_mask = checkChanListEquality(chan_stim,chan_list(i_chan));
            amp1_mask = amp == amp_list(1) & chan_mask' & array == i_arr;
            amp2_mask = amp == amp_list(2) & chan_mask' & array == i_arr;

            % if both, store in pair_idx
            if(sum(amp1_mask) >0 && sum(amp2_mask) > 0)
                pair_idx(end+1,:) = [find(amp1_mask),find(amp2_mask)];
                pair_arr(end+1,1) = i_arr;
            end
        end
    end
    
% projection magnitude for each amp -- compare within a channel pair
    f=figure(); hold on;
    % 1 = m1 = getColorFromList(1,2);
    % 0 = s1 = getColorFromList(1,3);
    idx_plot = [1,2,3]
    for i_pair = 1:size(pair_idx,1)
        plot([25,50],[mean(mag_list(pair_idx(i_pair,1),idx_plot)),mean(mag_list(pair_idx(i_pair,2),idx_plot))])
    end






