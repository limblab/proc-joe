%% set initial parameters

    input_data.folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\CObump\grill\';
%     mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp';
    mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
%     mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Pop_18E3\Array Map Files\6250-002085\SN 6250-002085.cmp';
    
    input_data.date = '20191015';
    input_data.array = 'arrayLeftS1';
    input_data.monkey = 'monkeyHan';
    input_data.ranBy = 'ranByJoe';
    input_data.lab = 6;
    input_data.mapFile = strcat('mapFile',mapFileName);
    input_data.task = 'taskCObump';

    pwd=cd;
    input_data.fileNum = 1;
  
    input_data.center_y = -33;
    input_data.center_x = 3;
    input_data.num_bins = 8;
    
%% make cds
    cd(input_data.folderpath)
    td_all = []
    file_name = dir('*nev*');
    for f = 1:numel(file_name)
        params.event_list = {'goCueTime';'bumpTime';'bumpDir'};
        params.trial_results = {'R'};
        params.extra_time = [1,2];
        params.include_ts = 0;
        params.exclude_units = [255];

        move_onset_params.pre_move_thresh = 1000;
        move_onset_params.min_s = 3;
        move_onset_params.start_idx_offset = -10;
        move_onset_params.max_rt_offset = 400;


        cds = commonDataStructure();
        cds.file2cds(strcat(input_data.folderpath,file_name(f).name),input_data.array,input_data.monkey,input_data.ranBy,...
            input_data.lab,input_data.mapFile,input_data.task,'recoverPreSync','ignoreJumps','ignoreFilecat');

        td_temp = parseFileByTrial(cds,params);
        td_all = [td_all,td_temp];
    end
    %% REMOVE ID FROM UNITS and make trial data
   
    td_all = parseFileByTrial(cds,params);
    td_all = stripSpikeSorting(td_all);
    td_all = getSpeed(td_all);
    td_all = removeBadTrials(td_all);
    
    if(td_all(1).bin_size < 0.01)
        % set it to 50ms
        td_all = binTD(td_all,ceil(0.05/td_all(1).bin_size));
    end
  
    
    td_all = getMoveOnset(td_all,move_onset_params);
%     td_all = removeBadTrials(td_all);


%% for each neuron, plot the mean firing rate for each unique BD (like a PD plot, but no cosine tuning)
    bump_dir = [90,270]; % relative to target direction
    tgt_dir = [0,45,90,135,180,225,270,315]; 
    mean_fr_bump = [];
    mean_fr_move = [];
    std_fr_move = [];
    z_score = [];
    chan_num = [];
    array_name = 'LeftS1';
    for unit = 1:size(td_all(1).([array_name,'_unit_guide']),1)
        chan_num(unit,1) = td_all(1).([array_name,'_unit_guide'])(unit,1);
        
        for td = 1:numel(tgt_dir)
            for bd = 1:numel(bump_dir)
                trial_mask_tgt = isEqual(round([td_all.target_direction]*180/pi),tgt_dir(td));
                trial_mask_bump = isEqual(round([td_all.bumpDir]), bump_dir(bd)+round([td_all.target_direction]*180/pi)) & trial_mask_tgt;
                trial_mask_move = isnan([td_all.bumpDir]) & trial_mask_tgt;
                
                
                fr_data_bump = [];
                fr_data_move = [];
                
                for t = 1:numel(td_all)
                    if(trial_mask_bump(t) == 1)
                        window_bump = td_all(t).idx_bumpTime + [0,4];
                        fr_data_bump(end+1,1) = sum(td_all(t).([array_name,'_spikes'])(window_bump(1):window_bump(2),unit))/(diff(window_bump)*td_all(t).bin_size);
                    elseif(trial_mask_move(t) == 1)
                        window_move = td_all(t).idx_movement_on + [1,5];
                        fr_data_move(end+1,1) = sum(td_all(t).([array_name,'_spikes'])(window_move(1):window_move(2),unit))/(diff(window_move)*td_all(t).bin_size);
                    end
                end
                mean_fr_bump(unit,td,bd) = mean(fr_data_bump);
                mean_fr_move(unit,td,bd) = mean(fr_data_move);
                std_fr_move(unit,td,bd) = std(fr_data_move)/sqrt(numel(fr_data_move));
                z_score(unit,td,bd) = (mean_fr_bump(unit,td,bd) - mean_fr_move(unit,td,bd))./(std_fr_move(unit,td,bd)+eps);
            end
           
        end
%         figure()
%         plot(tgt_dir,mean_fr_bump(unit,:,1) - mean_fr_move(unit,:,1),'k.','markersize',14)
%         hold on
%         plot(tgt_dir,mean_fr_bump(unit,:,2) - mean_fr_move(unit,:,2),'r.','markersize',14)
    end
 


%% plot z-score of responses during a tgt and bump dir
    map_filename = input_data.mapFile(8:end);
    map_data = loadMapFile(map_filename);
    min_max_z_score = 8;


    numColors = 1000;
    color_1 = [153,51,255]/256;
    color_2 = [255,153,51]/256;
    colors_1 = [(0:numColors-1)'/numColors,(0:numColors-1)'/numColors,(0:numColors-1)'/numColors].*color_1;
    colors_2 = [1-(0:numColors-1)'/numColors,1-(0:numColors-1)'/numColors,1-(0:numColors-1)'/numColors].*color_2;
    colors = [colors_2;colors_1];

    for td = 1:size(mean_fr_bump,2)
        for bd = 1:size(mean_fr_bump,3)
            fr_z_score = z_score(:,td,bd);
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

            f.Name = ['Han_20190215_centerHoldBumpHeatmap_tgt',num2str(tgt_dir(td)),'deg_bump',num2str(bump_dir(bd)),'deg'];
        %     saveFiguresLIB(f,input_data.folderpath,f.Name);

        end
    end


    
%% design stim pattern
        
    pattern_data = [];
    pattern_data.num_patterns = 10; % must be even (pairs of bio and nonbio)
    pattern_data.num_chans = 64;
    pattern_data.frequency = 200;
    pattern_data.max_amp = 60;
    
    num_amps = 15;
    pattern_MSD_list = [];
    
    pattern_data.amps = floor(linspace(pattern_data.max_amp/num_amps,pattern_data.max_amp,num_amps));

    pattern_data.pred_tgt_dir = [];
    pattern_data.pred_bump_dir = [];
    pattern_data.is_biomimetic = [];
    pattern_data.max_z_score = 10;
    pattern_data.z_score_cutoff = 0.5;
    pattern_data.pattern = [];

    % remove channels from chan_num which don't respond in any direction

    chan_mask = any(any(z_score > pattern_data.z_score_cutoff,2),3); 
    chans_use_all = chan_num(chan_mask);

    % pick num_patterns/2 directions from the unique list of tgt_dirs used
    tgt_dirs = unique([td_all.target_direction]); % in deg
    bump_dirs = [90,270];
    t_dirs = datasample(tgt_dirs,pattern_data.num_patterns/2,'Replace',true);
    b_dirs = datasample(bump_dirs,pattern_data.num_patterns/2);
    
    for d = 1:numel(t_dirs) % make two patterns per direction (1 non bio, 1 bio)
        td_idx = find(t_dirs(d) == tgt_dirs);
        bd_idx = find(b_dirs(d) == bump_dirs);

        pattern_data.pred_tgt_dir(end+1:end+2) = t_dirs(d);
        pattern_data.pred_bump_dir(end+1:end+2) = b_dirs(d);
        
        chans_use_dir = chans_use_all(abs(z_score(chan_mask,td_idx,bd_idx)) < 1000 & ~isnan(z_score(chan_mask,td_idx,bd_idx)));
        chans = datasample(chans_use_dir,pattern_data.num_chans,'Replace',false); % actual channel numbers
        
        % d is direction, chans represents chan numbers. Make bio pattern
        bio_pattern = []; bio_pattern.amp_ideal = zeros(numel(chans),1); bio_pattern.chans = chans;
        
        for c = 1:numel(chans) % for each channel, find z_score in direction and map to amplitude
            z_score_idx = find(chan_num == chans(c));
            z = z_score(z_score_idx,td_idx,bd_idx);
            z_all = reshape(z_score(z_score_idx,:,:),size(z_score,2)*size(z_score,3),1);
            
            norm_amp = min(1,max(0,z/pattern_data.max_z_score));
            norm_amp_all = min(1,max(0,z_all/pattern_data.max_z_score));
            
            bio_pattern.amp_ideal(c) = norm_amp*pattern_data.max_amp;
            bio_pattern.wave_num(c) = round(norm_amp*num_amps);
            bio_pattern.amp_all_dirs(:,c) = norm_amp_all*pattern_data.max_amp;
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
        
        pattern_MSD_list(end+1,1) = min(mean((nonbio_pattern.amp-bio_pattern.amp_all_dirs).^2,2)/(pattern_data.max_amp.^2));
    end

    for p = 1:pattern_data.num_patterns
        total_current(p) = sum(pattern_data.pattern{p}.amp);
    end
%     total_current
    
 %% heatmap of patterns using pattern_data
    
    map_filename = input_data.mapFile(8:end); map_data = loadMapFile(map_filename);
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
    
    