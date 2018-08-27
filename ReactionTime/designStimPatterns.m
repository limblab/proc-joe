%% use this script to design biomimetic stimulation patterns


%% determine filename and input data
    input_data.folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\Psychophysics\Han_20180801_COBump\';
%     inputData.folderpath = 'D:\Lab\Data\ReactionTime\Han_20180427_training\';
    input_data.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';

    input_data.task='taskCObump';
    input_data.ranBy='ranByJoseph'; 
    input_data.array1='arrayLeftS1'; 
    input_data.monkey='monkeyHan';
    input_data.labnum = 6;
    
    pwd=cd;
    cd(input_data.folderpath)
    fileList = dir('*s.nev*');
    cd(pwd)
    
%% load cds, convert to td, compute PDs for all units, determine if units are well tuned
    cds = commonDataStructure();
    cds.file2cds([input_data.folderpath fileList(1).name],input_data.task,input_data.ranBy,...
        input_data.monkey,input_data.labnum,input_data.array1,input_data.mapFileName,'recoverPreSync');
    cd(pwd);
    
% convert into td
    params.event_list = {'goCueTime';'tgtDir'};
    params.trial_results = {'R'};
    params.extra_time = [1,2];
    params.include_ts = 0;
    params.exclude_units = [1,255];
    td_all = parseFileByTrial(cds,params);
    td_all = removeBadTrials(td_all);
    td_all = getMoveOnsetAndPeak(td_all);
    td_all = removeBadTrials(td_all);
    [td_all] = removeBadNeurons(td_all,params);
    
% compute PD for all units
    params = [];
    params.array = 'LeftS1';
    params.window = {'idx_movement_on',0;'idx_movement_on',12};
    clear tcs
    [tcs,confBounds,fr,covar] = getNeuronTuning(td_all,'regress',params);
    % tcs = mean, modulation depth, pd (b0 + b1*cos(x-b2))
    tcs = tcs(:,1:3);
    
% determine if units are well tuned (conf bounds less than +- 45 deg)
    pd_idx = 3;
    conf_diff = angleDiff(confBounds{pd_idx}(:,1),tcs(:,pd_idx),1,0);
    well_tuned = abs(conf_diff) < pi/2;
    
%% randomly choose a direction
    dir_model = rand()*2*pi-pi;
    pd_elec = tcs(:,3);
    
%% determine pattern based on that direction and a number of electrodes
    input_data.num_elec = 20;
    input_data.well_tuned = well_tuned;
    input_data.pd_elec = pd_elec;
    input_data.dir_model = dir_model;
    input_data.dir_model_all = input_data.dir_model;
    input_data.chan_list = td_all(1).LeftS1_unit_guide(:,1);
%     input_data = rmfield(input_data,'freq_all_norm');
    
    pattern = {};
    pattern{1} = getBiomimeticPattern(td_all,input_data);
%     input_data.dir_model = input_data.dir_model - pi;
%     input_data.dir_model_all(2) = input_data.dir_model;
%     pattern{2} = getBiomimeticPattern(td_all,input_data);
    
%     input_data.freq_all_norm = pattern{1}.biomimetic_freq_norm(pattern{1}.biomimetic_freq_norm ~= 0);
%% build N patterns with the same number of electrodes  
    input_data.min_freq_norm = 16/330; % 16 is min cerestim, 330 is current max freq
%     input_data.freq_all_norm = linspace(input_data.min_freq_norm,1,15); % make empty if no freq binning needed
    input_data.chan_list = td_all(1).LeftS1_unit_guide(:,1);

    input_data.num_elec = 64;
    input_data.well_tuned = well_tuned;
    input_data.pd_elec = tcs(:,3);
    
    pattern = {};
    
    for i = 1:3
        input_data.dir_model = rand()*2*pi-pi;
        pattern{i} = getBiomimeticPattern(td_all,input_data);
        input_data.dir_model_all(i) = input_data.dir_model;
    end
    
    input_data = rmfield(input_data,'dir_model');
    
%% convert patterns into wave_mappings and freq_all for code

    wave_mappings = {}; % chan_num, wave_freq_norm, wave_num
    for pattern_idx = 1:numel(pattern)
        for bio = 1:2
            wave_mappings{end+1} = [];
            if(bio == 1) % biomimetic
                field_name = 'biomimetic_freq_norm';
            else % nonbiomimetic
                field_name = 'nonbiomimetic_freq_norm';
            end
            
            chan_ = pattern{pattern_idx}.chan_num(pattern{pattern_idx}.(field_name)>0);
            freq_ = pattern{pattern_idx}.(field_name)(pattern{pattern_idx}.(field_name)>0);
            if(isfield(input_data,'freq_all_norm'))
                for freq_idx = 1:numel(freq_)
                    freq_all_idx = find(freq_(freq_idx) == input_data.freq_all_norm);
                    wave_mappings{end}(freq_idx,1) = chan_(freq_idx);
                    wave_mappings{end}(freq_idx,2) = input_data.freq_all_norm(freq_all_idx);
                    wave_mappings{end}(freq_idx,3) = freq_all_idx;
                end
            else
                for freq_idx = 1:numel(freq_)
                    if(freq_(freq_idx) > 0)
                        wave_mappings{end}(freq_idx,2) = freq_(freq_idx);
                        wave_mappings{end}(freq_idx,3) = -1;
                        wave_mappings{end}(freq_idx,1) = chan_(freq_idx);
                    end
                end
            end
            
        end
    end


%% save wave_mapping and etc pattern
    save([input_data.folderpath,'Han_20180726_patterns_32elecAll'],...
        'wave_mappings','pattern','input_data');
    
      
    
    
%% convert wave_mappings into an array for > 16 electrode stimulation
    input_data.max_freq = 250;
%     input_data.max_freq_all = [150,250,330,500];
    input_data.train_length = 225; % in ms
    input_data.dt = 2.5;
    
    stim_array = {};
%     for t = 1:4
%         input_data.max_freq = input_data.max_freq_all(t);

        for i = 1:numel(wave_mappings)
            stim_array{i} = convertWaveMappings(wave_mappings{i}, input_data);
            stim_array{i}.stim_pattern = squeeze(stim_array{i}.stim_pattern);
             % plot ISI distribution for each channel
            pattern = squeeze(stim_array{i}.stim_pattern);
            ISI = {};
            [row,col] = find(pattern);
            figure();
            for chan = 1:size(pattern,1)
                mask = row == chan;
                ISI{chan} = diff(col(mask))*input_data.dt;
                ISI_desired = 1000/(wave_mappings{i}(chan,2)*input_data.max_freq);
                if(~isempty(ISI{chan}))
                    plot(chan,ISI{chan},'.','markersize',12,'color','b')
                    hold on
                    plot(chan,ISI_desired,'o','markersize',12,'color','r')
                    plot(chan,mean(ISI{chan}),'.','markersize',12,'color','k')
                end
            end
            xlabel('Chan index')
            ylabel('ISI (ms)')
            formatForLee(gcf)
            set(gca,'fontsize',14);
        end
%     end
    
    
   
    
    
    
    
    
    
    
%% plot heatmap of each pattern
    map_data = loadMapFile(input_data.mapFileName(8:end));
    color_list = [(1:1:64)'/64 zeros(64,1) zeros(64,1)];
    figure();
    plot([1,11,11,1,1],[1,1,11,11,1],'k','linewidth',1.5);
    hold on
    for c = 1:numel(chan_num)
        map_idx = find(map_data.chan == chan_num(c));
        pos = [map_data.row(map_idx),map_data.col(map_idx)];
        color_idx = max(1,floor(biomimetic_freq_norm(c)*size(color_list,1)));
        rectangle('Position',[pos(1),pos(2),1,1],'EdgeColor',color_list(color_idx,:),'FaceColor',color_list(color_idx,:));
    end
    
%     figure();
%     plot([1,11,11,1,1],[1,1,11,11,1],'k','linewidth',1.5);
%     hold on
%     for c = 1:numel(chan_num)/2
%         map_idx = find(map_data.chan == chan_num(c));
%         pos = [map_data.row(map_idx),map_data.col(map_idx)];
%         color_idx = max(1,floor(nonbiomimetic_freq_norm(c)*size(color_list,1)));
%         rectangle('Position',[pos(1),pos(2),1,1],'EdgeColor',color_list(color_idx,:),'FaceColor',color_list(color_idx,:));
%     end


    
    
   
    
    
%% plot heatmap of PDs
    pd_plot = ((pd_elec(elec_idx)))*180/pi; % 0 to 360 deg
    
    map_data = loadMapFile(input_data.mapFileName(8:end));
    c_data = zeros(10,10);
    for c = 1:numel(chan_num)
        map_idx = find(map_data.chan == chan_num(c));
        pos = [map_data.row(map_idx),11-map_data.col(map_idx)];
        c_data(pos(2),pos(1)) = pd_plot(c);
    end

    figure();
    imagesc(c_data);
    b = colorbar;
    b.Limits = [-180,180];
    colormap(jet)
    
    axis tight
    set(gca,'visible','off')
    


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    