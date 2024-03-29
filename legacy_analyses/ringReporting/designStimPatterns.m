%% use this script to design biomimetic stimulation patterns


%% determine filename and input data
    input_data.folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\CObump\Han_20200127_CObumpmove\';
%     input_data.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp';
    input_data.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';

    input_data.task='taskCObump';
    input_data.ranBy='ranByJoseph'; 
    input_data.array1='arrayLeftS1'; 
    input_data.monkey='monkeyHan';
    input_data.labnum = 6;
    
    pwd=cd;
    cd(input_data.folderpath)
    fileList = dir('*.nev*');
    cd(pwd)
    
%% load cds, convert to td, compute PDs for all units, determine if units are well tuned
    cds = commonDataStructure();
    cds.file2cds([input_data.folderpath fileList(1).name],input_data.task,input_data.ranBy,...
        input_data.monkey,input_data.labnum,input_data.array1,input_data.mapFileName,'recoverPreSync');
    cd(pwd);
    
  %% covert to td
    params.event_list = {'goCueTime';'tgtDir';'bumpDir';'bumpTime'};
    params.trial_results = {'R'};
    params.extra_time = [1,2];
    params.include_ts = 0;
    params.exclude_units = [];
    params.remove_nan_idx = true;
    params.nan_idx_names = {'idx_goCueTime','idx_startTime','idx_endTime'};
    td = parseFileByTrial(cds,params);
    td = stripSpikeSorting(td);
    td = removeBadTrials(td,params);

    params.start_idx = 'idx_goCueTime';
    params.end_idx = 'idx_endTime';
    td = getMoveOnsetAndPeak(td,params);
    [td] = removeBadNeurons(td,params);
   
    td_all = td([td.idx_goCueTime] > [td.idx_bumpTime]);
    
    for tr = 1:numel(td_all)
        td_all(tr).bumpDir = mod(td_all(tr).bumpDir,360);
    end
%% for each neuron, get the mean firing rate for each unique BD and during baseline
    bump_dir = unique([td_all.bumpDir]);
    bump_dir(isnan(bump_dir)) = [];
    
    mean_fr_bump = [];
    mean_fr_baseline = [];
    mean_fr_move = [];
    std_fr_baseline_all_dir = [];
    mean_fr_baseline_all_dir = [];
    z_score = [];
    chan_num = [];
    array_name = 'LeftS1';
    
    
    for unit = 1:size(td_all(1).([array_name,'_unit_guide']),1)
        for bd = 1:numel(bump_dir)
            trial_mask = isEqual(round([td_all.bumpDir]), bump_dir(bd));
            fr_data_bump = [];
            fr_data_baseline = [];
            fr_data_move = [];
            td_temp = td_all(trial_mask);
            for t = 1:numel(td_temp)
                window_bump = td_temp(t).idx_bumpTime + [0,5];
                window_baseline = td_temp(t).idx_bumpTime - [40,20];
%                 window_move = td_temp(t).idx_movement_on + [5,25];
                
%                 window_baseline = [td_all(t).idx_startTime, td_all(t).idx_bumpTime-2];
                fr_data_bump(t) = sum(td_temp(t).([array_name,'_spikes'])(window_bump(1):window_bump(2),unit))/(diff(window_bump)*td_temp(t).bin_size);
                fr_data_baseline(t) = sum(td_temp(t).([array_name,'_spikes'])(window_baseline(1):window_baseline(2),unit))/(diff(window_baseline)*td_temp(t).bin_size);
%                 fr_data_move(t) = sum(td_temp(t).([array_name,'_spikes'])(window_move(1):window_move(2),unit))/(diff(window_move)*td_temp(t).bin_size);
            end
            mean_fr_bump(unit,bd) = mean(fr_data_bump);
            mean_fr_baseline(unit,bd) = mean(fr_data_baseline);
            
%             mean_fr_move(unit,bd) = mean(fr_data_move);
            std_fr_baseline(unit,bd) = std(fr_data_baseline);
%             z_score(unit,bd) = (mean_fr_bump(unit,bd) - mean_fr_baseline(unit,bd))./(std_fr_baseline(unit,bd)+eps);
        end
        chan_num(unit,1) = td_all(1).([array_name,'_unit_guide'])(unit,1);
        std_fr_baseline_all_dir(unit,1) = std(mean_fr_baseline(unit,:));

        mean_fr_baseline_all_dir(unit,1) = mean(mean_fr_baseline(unit,:));
        z_score(unit,:) = (mean_fr_bump(unit,:) - mean_fr_baseline_all_dir(unit))./(std_fr_baseline_all_dir(unit)+eps);
        
%         figure()
%         plot(bump_dir,z_score(unit,:),'.','markersize',14)
    end

% mean_fr = mean_fr_bump - mean(mean_fr_bump,2);
% mean_fr_baseline = mean_fr_baseline - mean(mean_fr_baseline,2);

%% plot z-score of responses during a bump dir
    map_filename = input_data.mapFileName(8:end);
    map_data = loadMapFile(map_filename);
    min_max_z_score = 4;

    numColors = 1000;
    color_1 = [153,51,255]/256;
    color_2 = [255,153,51]/256;
    colors_1 = [(0:numColors-1)'/numColors,(0:numColors-1)'/numColors,(0:numColors-1)'/numColors].*color_1;
    colors_2 = [1-(0:numColors-1)'/numColors,1-(0:numColors-1)'/numColors,1-(0:numColors-1)'/numColors].*color_2;
    colors = [colors_2;colors_1];

    for bd = 1:size(mean_fr_bump,2)
        fr_z_score = z_score(:,bd);
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

            rectangle('Position',[map_data.row(map_idx),map_data.col(map_idx),1,1],'FaceColor',colors(color_idx,:), 'EdgeColor','none');
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

        f.Name = ['Duncan_20190410_RRBumpHeatmap_',num2str(bump_dir(bd)),'deg'];
    %     saveFiguresLIB(f,input_data.folderpath,f.Name);

    end


%% design stim pattern
        
    pattern_data = [];
    pattern_data.num_patterns = 10000; % must be even (pairs of bio and nonbio)
    pattern_data.num_chans = 32;
    pattern_data.frequency = 200;
    pattern_data.max_amp = 60;
    
    num_amps = 15;
    pattern_data.amps = floor(linspace(pattern_data.max_amp/num_amps,pattern_data.max_amp,num_amps));

    pattern_MSD_list = [];
    
    pattern_data.pred_dirs = [];
    pattern_data.is_biomimetic = [];
    pattern_data.max_z_score = 4;
    pattern_data.z_score_cutoff = 1;
    pattern_data.pattern = [];

    % remove channels from chan_num which don't respond in any direction

    chan_mask = any(z_score > pattern_data.z_score_cutoff,2); 
    chans_use_all = chan_num(chan_mask);

    % pick num_patterns/2 directions from the unique list of tgt_dirs used
%     bump_dirs = unique([td_all.bumpDir]); % in deg
    bump_dirs = [0,90,180,270]; % only do these....
    dirs = datasample(bump_dirs,pattern_data.num_patterns/2,'Replace',true);
    for d = dirs % make two patterns per direction (1 non bio, 1 bio)
        bd_idx = find(bump_dirs == d);

        pattern_data.pred_dirs(end+1:end+2) = d;
        chans_use_dir = chans_use_all(abs(z_score(chan_mask,bd_idx)) < 1000 & ~isnan(z_score(chan_mask,bd_idx)));
        chans = datasample(chans_use_dir,pattern_data.num_chans,'Replace',false); % actual channel numbers
        
        % d is direction, chans represents chan numbers. Make bio pattern
        bio_pattern = []; bio_pattern.amp_ideal = zeros(numel(chans),1); bio_pattern.chans = chans;
        
        for c = 1:numel(chans) % for each channel, find z_score in direction and map to amplitude
            z_score_idx = find(chan_num == chans(c));
            z = z_score(z_score_idx,bd_idx);
            
            norm_amp = min(1,max(0,z/pattern_data.max_z_score));
            bio_pattern.amp_ideal(c) = norm_amp*pattern_data.max_amp;
            bio_pattern.wave_num(c) = round(norm_amp*num_amps);
            if(bio_pattern.wave_num(c) == 0)
                bio_pattern.amp(c) = 0;
            else
                bio_pattern.amp(c) = pattern_data.amps(bio_pattern.wave_num(c));
            end
        end
        
        % make nonbio_pattern by scrambling amps on same channels
        nonbio_pattern = bio_pattern;
        scramble_idx = randperm(numel(bio_pattern.amp_ideal));
        nonbio_pattern.amp_ideal = nonbio_pattern.amp_ideal(scramble_idx);
        nonbio_pattern.wave_num = nonbio_pattern.wave_num(scramble_idx);
        nonbio_pattern.amp = nonbio_pattern.amp(scramble_idx);
        
        pattern_data.pattern{end+1} = bio_pattern;
        pattern_data.is_biomimetic(end+1) = 1;
        pattern_data.pattern{end+1} = nonbio_pattern;
        pattern_data.is_biomimetic(end+1) = 0;
        
        pattern_MSD_list(end+1,1) = mean((nonbio_pattern.amp-bio_pattern.amp).^2)/(pattern_data.max_amp.^2);
    end

    for p = 1:pattern_data.num_patterns
        total_current(p) = sum(pattern_data.pattern{p}.amp);
    end
%     total_current
    
 %% heatmap of patterns using pattern_data
    
    map_filename = input_data.mapFileName(8:end); map_data = loadMapFile(map_filename);
    numColors = 100;
    base_color = [153,51,255]/256;
    colors = [(0:numColors-1)'/numColors,(0:numColors-1)'/numColors,(0:numColors-1)'/numColors].*base_color;
    % chan_num_use = chan_num(randperm(96,64));

    for pattern_idx = 1:pattern_data.num_patterns
        figure();
        hold on
        
        for chan_idx = 1:pattern_data.num_chans
            map_idx = find(pattern_data.pattern{pattern_idx}.chans(chan_idx) == map_data.chan);
            
            color_idx = max(1,min(size(colors,1),ceil(pattern_data.pattern{pattern_idx}.amp(chan_idx)*size(colors,1)/pattern_data.max_amp)));

            rectangle('Position',[map_data.row(map_idx),map_data.col(map_idx),1,1],'FaceColor',colors(color_idx,:), 'EdgeColor','none');
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
        b.Ticks = [-1,-0.5,0,0.5,1];
        b.Ticks = [0,0.25,0.5,0.75,1];
        b.TickLabels = {'0','','0.5','','1'};
        b.Label.String = 'Normalized stimulation amplitude';
        b.Label.FontSize = 16;

    %     f.Name = ['Han_20171030_COactpas_bumpHeatmap_',num2str(bump_dir(bd)),'deg'];
    %     f.Name = ['Han_20171030_COactpas_centerHoldHeatmap'];
        f.Name = ['Duncan_20190410_RRheatmap_pattern',num2str(pattern_idx)];
    %     saveFiguresLIB(f,input_data.folderpath,f.Name);
    end
    
    
    
    
    
    
    
%% similarity between patterns
    similarity_metric = zeros(numel(pattern_data.pattern),numel(tgt_dir));
    for patt = 1:numel(pattern_data.pattern)
        for bd = 1:numel(tgt_dir)
            % build pattern in tgt_dir(bd) for chans in
            % pattern_data.pattern{patt}.chans
            amp = [];
            amp_ideal = [];
            wave_num = [];
            for c = 1:numel(pattern_data.pattern{patt}.chans)
                z_score_idx = find(chan_num == pattern_data.pattern{patt}.chans(c));
                z = z_score(z_score_idx,bd);
                
                norm_amp = min(1,max(0,z/pattern_data.max_z_score));
                amp_ideal(c) = norm_amp*pattern_data.max_amp;
                wave_num = round(norm_amp*numel(pattern_data.amps));
                if(wave_num == 0)
                    amp(c) = 0;
                else
                    amp(c) = pattern_data.amps(wave_num);
                end 
            end
            
            % get similarity metric (mean difference in amps on those chans
            similarity_metric(patt,bd) = mean((amp-pattern_data.pattern{patt}.amp).^2);
        end
        
    end
    
    
    
    
    
    
    
    

%% plot heatmap of stim pattern
    make_nonbiomimetic = 1;

    map_filename = input_data.mapFileName(8:end);
    map_data = loadMapFile(map_filename);
    max_freq = 4;
    numColors = 1000;
    color_1 = [153,51,255]/256;
    color_2 = [255,153,51]/256;
    colors_1 = [(0:numColors-1)'/numColors,(0:numColors-1)'/numColors,(0:numColors-1)'/numColors].*color_1;
    colors_2 = [1-(0:numColors-1)'/numColors,1-(0:numColors-1)'/numColors,1-(0:numColors-1)'/numColors].*color_2;
    colors = colors_1;
    rng('shuffle','twister')
    % chan_num_use = chan_num(randperm(96,64));

    for bd = 1%:size(mean_fr_bump,2)
        figure();
        hold on
        chan_mapping = [];
        for elec = 1:size(chan_num_use,1)
            map_idx = find(chan_num_use(elec) == map_data.chan);
    %         freq_norm = (mean_fr_bump(chan_num_use(elec),bd) - mean_fr_baseline(chan_num_use(elec),bd))./std_fr_baseline(chan_num_use(elec),bd);
            freq_norm = (mean_fr_bump(chan_num_use(elec),bd) - mean_fr_baseline(chan_num_use(elec),bd))./max(mean_fr_bump(chan_num_use(elec),:) - mean_fr_baseline(chan_num_use(elec),:));      
    %         freq_norm = freq_norm/max_freq;
            freq_norm = max(freq_norm,0);
            if(freq_norm >= 0 && freq_norm < 100 && ~isnan(freq_norm))
                chan_mapping(end+1,:) = [chan_num_use(elec),freq_norm];
            end
        end

        if(make_nonbiomimetic)
            chan_mapping(:,1) = chan_mapping(randperm(size(chan_mapping,1)),1);
        end

        for elec = 1:size(chan_mapping,1)
            map_idx = find(chan_mapping(elec,1) == map_data.chan);
            if(chan_mapping(elec,2) >= 0 && chan_mapping(elec,2) < 100 && ~isnan(chan_mapping(elec,2)))
                color_idx = max(1,min(size(colors,1),ceil(chan_mapping(elec,2)*size(colors,1))));

                rectangle('Position',[map_data.row(map_idx),map_data.col(map_idx),1,1],'FaceColor',colors(color_idx,:), 'EdgeColor','none');
                hold on
            end
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
        b.Ticks = [-1,-0.5,0,0.5,1];
        b.Ticks = [0,0.25,0.5,0.75,1];
        b.TickLabels = {'0','','0.5','','1'};
        b.Label.String = 'Normalized stimulation amplitude';
        b.Label.FontSize = 16;

    %     f.Name = ['Han_20171030_COactpas_bumpHeatmap_',num2str(bump_dir(bd)),'deg'];
    %     f.Name = ['Han_20171030_COactpas_centerHoldHeatmap'];
        f.Name = ['Duncan_20190215_centerHoldStimHeatmapBIO_',num2str(tgt_dir(bd)),'deg'];
    %     saveFiguresLIB(f,input_data.folderpath,f.Name);
    end

%% put desired axis into wave_mappings
    axis_idx = 3;
    wave_mappings = {};
    wave_mappings{1}(:,:) = [chan_num,mean_fr_norm(:,axis_idx)];
    wave_mappings{2}(:,:) = [chan_num,mean_fr_norm(:,axis_idx+4)];

    wave_mappings{1}(wave_mappings{1}(:,2)<0.15,:) = [];
    wave_mappings{2}(wave_mappings{2}(:,2)<0.15,:) = [];

%% dissimilarity metric
    % euclidean distance between electrode and nearest electrode multiplied
    % by normalized freq. 


