%% set initial parameters

    input_data.folderpath = 'D:\Lab\Data\CObumpmove\Han_20200608\';
%     mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp';
    mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
%     mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Pop_18E3\Array Map Files\6250-002085\SN 6250-002085.cmp';
    
    input_data.date = '20200608';
    input_data.array = 'arrayLeftS1';
    input_data.monkey = 'monkeyHan';
    input_data.ranBy = 'ranByJoe';
    input_data.lab = 6;
    input_data.mapFileName = strcat('mapFile',mapFileName);
    input_data.task = 'taskCObump';

    pwd=cd;
    input_data.fileNum = 1;
  
    input_data.center_y = -33;
    input_data.center_x = 3;
    input_data.num_bins = 8;
    
%% make cds
    cd(input_data.folderpath)

    file_name = dir('*nev*');
    
    params.event_list = {'goCueTime';'bumpTime';'bumpDir'};
    params.trial_results = {'R'};
    params.extra_time = [1,2];
    params.include_ts = 0;
    params.exclude_units = [255];

    move_onset_params.pre_move_thresh = 1000;
    move_onset_params.min_s = 3;
    move_onset_params.start_idx_offset = 0;
    move_onset_params.max_rt_offset = 400;
    

    cds = commonDataStructure();
    cds.file2cds(strcat(input_data.folderpath,file_name(1).name),input_data.array,input_data.monkey,input_data.ranBy,...
        input_data.lab,input_data.mapFileName,input_data.task,'recoverPreSync','ignoreJumps','ignoreFilecat');
    
    
% REMOVE ID FROM UNITS and make trial data
   
    td_all = parseFileByTrial(cds,params);
    td_all = stripSpikeSorting(td_all);
    td_all = getNorm(td_all,struct('signals',{'vel'}));
    td_all = removeBadTrials(td_all);
    
    if(td_all(1).bin_size < 0.01)
        % set it to 50ms
        td_all = binTD(td_all,ceil(0.05/td_all(1).bin_size));
    end
  
    td_all = getSpeed(td_all);
    td_all = getMoveOnsetAndPeak(td_all,move_onset_params);
    td_all = removeBadTrials(td_all);
    td_all(isnan([td_all.idx_movement_on])) = [];
    
%% for each neuron, plot the mean firing rate for each unique BD (like a PD plot, but no cosine tuning)
    bump_dir = [90,270]; % relative to target direction
    tgt_dir = [0,90,180,270]; 
    mean_fr_bump_move = [];
    mean_fr_move = [];
    std_fr_move = [];
    z_score = [];
    chan_num = [];
    array_name = 'LeftS1';
    
    move_to_bump = [td_all.idx_bumpTime] - [td_all.idx_movement_on];
    move_to_bump(isnan(move_to_bump) | move_to_bump < 0 ) = [];
    move_to_bump_avg = floor(mean(move_to_bump));
    bin_size = td_all(1).bin_size;

    for i_unit = 1:size(td_all(1).([array_name,'_unit_guide']),1)
        chan_num(i_unit,1) = td_all(1).([array_name,'_unit_guide'])(i_unit,1);
%         i_unit = find(td_all(1).([array_name,'_unit_guide'])(:,1) == temp);
%         
        for i_tgt_dir = 1:numel(tgt_dir)
            for i_bump_dir = 1:numel(bump_dir)
                trial_mask_tgt = isEqual(round([td_all.target_direction]*180/pi),tgt_dir(i_tgt_dir));
                trial_mask_bump = isEqual(round([td_all.bumpDir]), bump_dir(i_bump_dir)+round([td_all.target_direction]*180/pi)) & trial_mask_tgt;
                trial_mask_bump_move = trial_mask_bump & [td_all.idx_bumpTime] > [td_all.idx_goCueTime];
                trial_mask_move = isnan([td_all.bumpDir]) & trial_mask_tgt;
                
                
                fr_data_bump_move = [];
                fr_data_move = [];
                fr_data_pre_bump = [];
                for i_trial = 1:numel(td_all)
                    if(trial_mask_bump_move(i_trial) == 1)
                        window_bump = td_all(i_trial).idx_bumpTime + floor([0,0.150]/bin_size); % [s,s]
                        window_pre = td_all(i_trial).idx_bumpTime + floor([-0.15,-0.1]/bin_size);
                        fr_data_bump_move(end+1,1) = sum(td_all(i_trial).([array_name,'_spikes'])(window_bump(1):window_bump(2),i_unit))/(diff(window_bump)*td_all(i_trial).bin_size);
                        fr_data_pre_bump(end+1,1) = sum(td_all(i_trial).([array_name,'_spikes'])(window_pre(1):window_pre(2),i_unit))/(diff(window_pre)*td_all(i_trial).bin_size);
                    elseif(trial_mask_move(i_trial) == 1)
                        window_move = td_all(i_trial).idx_movement_on + move_to_bump_avg + floor([0.,0.150]/bin_size);
                        fr_data_move(end+1,1) = sum(td_all(i_trial).([array_name,'_spikes'])(window_move(1):window_move(2),i_unit))/(diff(window_move)*td_all(i_trial).bin_size);
                    end
                end
                mean_fr_bump_move(i_unit,i_tgt_dir,i_bump_dir) = mean(fr_data_bump_move);
                mean_fr_pre(i_unit,i_tgt_dir,i_bump_dir) = mean(fr_data_pre_bump);
                mean_fr_move(i_unit,i_tgt_dir) = mean(fr_data_move);
                std_fr_move(i_unit,i_tgt_dir) = std(fr_data_move);
                z_score(i_unit,i_tgt_dir,i_bump_dir) = (mean_fr_bump_move(i_unit,i_tgt_dir,i_bump_dir) - mean_fr_move(i_unit,i_tgt_dir))./(std_fr_move(i_unit,i_tgt_dir)+eps);
            end
           
        end
%         figure()
%         plot(tgt_dir*pi/180,z_score(i_unit,:,1),'b.','markersize',14)
%         hold on
%         plot(tgt_dir*pi/180,z_score(i_unit,:,2),'r.','markersize',14)
        
    end
 


%% plot z-score of responses during a tgt and bump dir
    map_filename = input_data.mapFile(8:end);
    map_data = loadMapFile(map_filename);
    min_max_z_score = 1.5;


    numColors = 1000;
    color_1 = [153,51,255]/256;
    color_2 = [255,153,51]/256;
    colors_1 = [(0:numColors-1)'/numColors,(0:numColors-1)'/numColors,(0:numColors-1)'/numColors].*color_1;
    colors_2 = [1-(0:numColors-1)'/numColors,1-(0:numColors-1)'/numColors,1-(0:numColors-1)'/numColors].*color_2;
    colors = [colors_2;colors_1];

    for i_tgt_dir = 1:size(mean_fr_bump_move,2)
        for i_bump_dir = 1:size(mean_fr_bump_move,3)
            fr_z_score = z_score(:,i_tgt_dir,i_bump_dir);
            chan_num_z_score = chan_num;

            % remove inf entries
            chan_num_z_score(abs(fr_z_score) > 10000 | isnan(fr_z_score)) = [];
            fr_z_score(abs(fr_z_score) > 10000 | isnan(fr_z_score)) = [];
            disp(max(abs(fr_z_score)))
            figure();
            hold on
            for elec = 1:size(chan_num_z_score,1)
                map_idx = find(chan_num_z_score(elec) == map_data.chan);
                color_idx = max(1,min(size(colors,1),ceil((fr_z_score(elec)+min_max_z_score)/(2*min_max_z_score)*size(colors,1))));

                rectangle('Position',[map_data.col(map_idx),map_data.row(map_idx),1,1],'FaceColor',colors(color_idx,:), 'EdgeColor','none');
                hold on
            end
            plot([1,11,11,1,1],[1,1,11,11,1],'k','linewidth',1.5)

            f=gcf;
            set(gca,'visible','off')
            xlim([1,11])
            ylim([1,11])
            axis square

            b = colorbar;
            colormap(colors);
            set(gca,'Visible','off');
            b.FontSize = 14;
            b.Ticks = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];
            b.TickDirection = 'out';
            b.Ticks = min_max_z_score*[-1,-0.5,0,0.5,1];
            b.Ticks = [0,0.25,0.5,0.75,1];
            b.TickLabels = {['<-',num2str(min_max_z_score)],'','0','',['>',num2str(min_max_z_score)]};
            b.Label.String = 'Normalized change in firing rate';
            b.Label.FontSize = 16;

            f.Name = ['Han_20190215_centerHoldBumpHeatmap_tgt',num2str(tgt_dir(i_tgt_dir)),'deg_bump',num2str(bump_dir(i_bump_dir)),'deg'];
        %     saveFiguresLIB(f,input_data.folderpath,f.Name);

        end
    end

%% design stim pattern 
    min_MSD_percentile = 75; % percentage
    num_scrambles = 10000;
    
    pattern_data = [];
    pattern_data.num_patterns = 8; % must be even (pairs of bio and nonbio)
    pattern_data.num_chans = 40;
    
    pattern_MSD_list = [];
    
    pattern_data.pred_tgt_dir = [];
    pattern_data.pred_bump_dir = [];
    pattern_data.is_biomimetic = [];
    
    pattern_data.max_z_score = 2;
    pattern_data.z_score_cutoff = 0.2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pattern_data.bump_z_score_min_diff = 0.2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pattern_data.pattern = [];

    % remove channels from chan_num which don't respond in any direction

    chan_mask = any(any(z_score > pattern_data.z_score_cutoff,2),3); 
    diff_is_large = abs(z_score(:,:,2)-z_score(:,:,1)) > pattern_data.bump_z_score_min_diff;
    chans_use_all = chan_num(chan_mask);

    % pick num_patterns/2 directions from the unique list of tgt_dirs used
    tgt_dirs = unique([td_all.target_direction]); % in deg
    bump_dirs = [90,270];
    t_dirs = [0,0,180,180]*pi/180; %datasample(tgt_dirs,pattern_data.num_patterns/2,'Replace',true);
    b_dirs = [90,270,90,270]; %datasample(bump_dirs,pattern_data.num_patterns/2);
    
    for d = 1:numel(t_dirs) % make two patterns per direction (1 non bio, 1 bio)
        td_idx = find(t_dirs(d) == tgt_dirs);
        bd_idx = find(b_dirs(d) == bump_dirs);
        
        chans_use_dir = chans_use_all(abs(z_score(chan_mask,td_idx,bd_idx)) < 1000 & ~isnan(z_score(chan_mask,td_idx,bd_idx)));
        chans = datasample(chans_use_dir,pattern_data.num_chans,'Replace',false); % actual channel numbers
        
        % d is direction, chans represents chan numbers. Make bio pattern
        bio_pattern = []; bio_pattern.stim_norm = zeros(numel(chans),1); bio_pattern.chans = chans;
        bio_pattern.stim_norm_all_dirs = zeros(numel(chans),size(z_score,2)*size(z_score,3));
        
        for c = 1:numel(chans) % for each channel, find z_score in direction and map to amplitude
            z_score_idx = find(chan_num == chans(c));
            z = z_score(z_score_idx,td_idx,bd_idx);
            z_all = reshape(z_score(z_score_idx,:,:),size(z_score,2)*size(z_score,3),1);
            
            if(~diff_is_large(z_score_idx,td_idx))
                z = -3;
            end
            
            stim_norm = min(1,max(0,z/pattern_data.max_z_score));
            stim_norm_all_dirs = min(1,max(0,z_all/pattern_data.max_z_score));
            
            bio_pattern.stim_norm(c) = stim_norm;
            bio_pattern.stim_norm_all_dirs(c,:) = stim_norm_all_dirs;
            
            if(bio_pattern.stim_norm(c) <= 0)
                bio_pattern.stim_norm(c) = 0;
            end
        end
        
        pattern_data.is_biomimetic(end+1) = 1;
        pattern_data.pred_tgt_dir(end+1) = t_dirs(d);
        pattern_data.pred_bump_dir(end+1) = b_dirs(d);
        pattern_data.pattern{end+1} = bio_pattern;
    end
    
    
   % make nonbio_pattern by scrambling amps on same channels
    for i_pattern = 1:numel(pattern_data.pattern)
        % first scramble nonbio pattern a ton to generate distribution of
        % MSD
        min_MSD_list = zeros(num_scrambles,1);
        for i_scramble = 1:num_scrambles
            scramble_idx = randperm(numel(pattern_data.pattern{i_pattern}.chans));
            scramble_stim_norm = pattern_data.pattern{i_pattern}.stim_norm(scramble_idx);

            min_MSD = 1000000;
            for i_zscore = 1:size(z_all,2)
                MSD = mean((scramble_stim_norm-pattern_data.pattern{i_pattern}.stim_norm_all_dirs(:,i_zscore)).^2);
                if(MSD < min_MSD)
                    min_MSD = MSD;
                end
            end
            min_MSD_list(i_scramble) = min_MSD;
        end
        threshold_min_MSD = prctile(min_MSD_list,min_MSD_percentile);
        
        % scramble pattern until MSD > threshold
        curr_min_MSD = threshold_min_MSD-1;
        pattern_counter = 0;
        best_min_MSD = 0;
        best_scramble_idx = 0;
        while(curr_min_MSD < threshold_min_MSD && pattern_counter < 100)
            scramble_idx = randperm(numel(pattern_data.pattern{i_pattern}.chans));
            scramble_stim_norm = pattern_data.pattern{i_pattern}.stim_norm(scramble_idx);
            curr_min_MSD = 1000000;
            for i_zscore = 1:numel(z_all)
                MSD = mean((scramble_stim_norm-pattern_data.pattern{i_pattern}.stim_norm_all_dirs(:,i_zscore)).^2);
                if(MSD < curr_min_MSD)
                    curr_min_MSD = MSD;
                end
            end
            
            if(curr_min_MSD > best_min_MSD)
                best_scramble_idx = scramble_idx;
                best_min_MSD = curr_min_MSD;
            end
            pattern_counter = pattern_counter + 1;
        end
        
        nonbio_pattern = pattern_data.pattern{i_pattern};
        nonbio_pattern.stim_norm = pattern_data.pattern{i_pattern}.stim_norm(best_scramble_idx);
        
        pattern_data.is_biomimetic(end+1) = 0;
        pattern_data.pred_tgt_dir(end+1) = pattern_data.pred_tgt_dir(i_pattern);
        pattern_data.pred_bump_dir(end+1) = pattern_data.pred_bump_dir(i_pattern);
        pattern_data.pattern{end+1} = nonbio_pattern;
    end

%% plot heatmaps
    heatmap_input_data = [];
    heatmap_input_data.map_filename = input_data.mapFileName;
    heatmap_input_data.num_colors = 100;
    heatmap_input_data.base_color = [153,51,255]/256; % max stim value is this color
    
    plotPatternDataHeatmap(pattern_data, heatmap_input_data);

    total_charge = zeros(numel(pattern_data.pattern),1);
    num_chans = zeros(numel(pattern_data.pattern),1);
    for i_pattern = 1:numel(pattern_data.pattern)
        total_charge(i_pattern) = sum(pattern_data.pattern{i_pattern}.stim_norm)*330*30; % *freq*stim_amp 
        num_chans(i_pattern) = sum(pattern_data.pattern{i_pattern}.stim_norm > 0);
    end
    disp(total_charge)
    disp(num_chans)
    
  
%% convert pattern data to stim_array for the cerestim
    stim_params = [];
    stim_params.IPI = 2.5/1000; % s
    stim_params.max_freq = 330; % Hz
    stim_params.min_freq = 0; % Hz
    stim_params.train_length = 0.25; % s
    stim_array_data = makeStimArrayWrapper(pattern_data,stim_params);
    % make stim_array which is formatted properly
    stim_array = cell(numel(stim_array_data.stim_array),1);
    for i_pattern = 1:numel(stim_array)
        stim_array{i_pattern}.stim_pattern = stim_array_data.stim_array{i_pattern};
        stim_array{i_pattern}.chans = stim_array_data.chans{i_pattern};
    end
    
    
%% plot mean IPI and list of IPIs with desired IPI for each electrode
    figure();
    for i_pattern = 1:numel(stim_array_data.stim_array)
        subplot(4,2,i_pattern); hold on
        for i_elec = 1:size(stim_array_data.stim_array{i_pattern},1)
            IPI_list = diff(find(stim_array_data.stim_array{i_pattern}(i_elec,:)))*stim_params.IPI;
            if(~isempty(IPI_list))
                mean_IPI = mean(IPI_list);
                IPI_list = unique(IPI_list);

                plot(i_elec,mean_IPI,'r.','markersize',12);
                plot(i_elec,IPI_list,'bo','markersize',8);
            end
            
            % get desired IPI
            chan_idx = find(pattern_data.pattern{i_pattern}.chans == stim_array_data.chans{i_pattern}(i_elec));
            desired_IPI = 1/(pattern_data.pattern{i_pattern}.stim_norm(chan_idx)*stim_params.max_freq);
            plot(i_elec,desired_IPI,'k*','markersize',8);
        end
    end
    


    
%% convert patterns into wave_mappings and freq_all for code

    wave_mappings = {}; % chan_num, wave_freq_norm, wave_num
    freq_all_norm = pattern_data.freqs/pattern_data.max_freq;
    
    for i_patt = 1:numel(pattern_data.pattern)
        wave_mappings{end+1} = [];
        field_name = 'wave_num';
        chan_ = pattern_data.pattern{i_patt}.chans(pattern_data.pattern{i_patt}.(field_name)>0);
        wave_ = pattern_data.pattern{i_patt}.(field_name)(pattern_data.pattern{i_patt}.(field_name)>0);
            
        for i_wave = 1:numel(wave_)
            wave_mappings{end}(i_wave,1) = chan_(i_wave);
            wave_mappings{end}(i_wave,2) = freq_all_norm(wave_(i_wave));
            wave_mappings{end}(i_wave,3) = wave_(i_wave);
        end
        
    end 