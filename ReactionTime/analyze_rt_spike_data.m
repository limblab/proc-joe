%% assumes we have a td with idx_goCueTime and idx_movement_on (td_reward from analyze_reactionTime)

    map_data = loadMapFile(inputData.mapFileName(8:end));
%% look at bump data (does # of spikes correlate with RT during bumps)
    
    td_bump = td(~isnan([td.bumpDir])); 
    rt = []; num_spikes = [];
    
    
    for tr = 1:numel(td_bump)
        window = td_bump(tr).idx_goCueTime + [0,15];
        rt(tr) = td_bump(tr).idx_movement_on - td_bump(tr).idx_goCueTime;
        num_spikes(tr,:) = sum(sum(td_bump(tr).LeftS1_spikes(window(1):window(2),:)));

    end

    remove_idx = isnan(rt);
    rt(remove_idx) = []; 
    num_spikes(remove_idx,:) = [];

    plot(rt,num_spikes,'.','markersize',20)
    
%% look at stim data
    neighbor_size = 5;
    stim_chan = 44;
    td_stim = td(~isnan([td.stimCode]) & [td.stimCode] == 4); 
    
    map_idx_stim_chan = find(map_data.chan == stim_chan);
    stim_chan_pos = [map_data.row(map_idx_stim_chan),map_data.col(map_idx_stim_chan)];
    
    dists = [];
    for chan = 1:size(td_stim(1).LeftS1_unit_guide,1)
        map_idx_rec_chan = find(map_data.chan == td_stim(1).LeftS1_unit_guide(chan,1));
        rec_chan_pos = [map_data.row(map_idx_rec_chan),map_data.col(map_idx_rec_chan)];
        dists(chan) = sqrt(sum((stim_chan_pos-rec_chan_pos).^2));
    end
    
    dist_mask = dists <= neighbor_size;
    
    rt = []; num_spikes = [];
    for tr = 1:numel(td_stim)
        window = td_stim(tr).idx_goCueTime + [0,20];
        rt(tr) = td_stim(tr).idx_movement_on - td_stim(tr).idx_goCueTime;
        
        
        num_spikes(tr) = sum(sum(td_stim(tr).LeftS1_spikes(window(1):window(2),dist_mask)));

    end

    remove_idx = isnan(rt);
    rt(remove_idx) = []; 
    num_spikes(remove_idx) = [];

    plot(rt,num_spikes,'.','markersize',20)



%% fails and rewards for the same stim code

    td_use = td_all([td_all.result] == 'R' | [td_all.result] == 'F');
    td_use = td_use(~isnan([td_use.stimCode]));
    
    % get movement onset
    params.field_idx = 1;
    params.start_idx_offset = 10*10;
    params.be_aggressive = 1;
    params.which_field = 'acc';
    params.max_rt_offset = 400;
    % Han's parameters
    params.threshold_acc = 35; % absolute threshold on acceleration, using this instead of threshold_mult
    params.min_s = 1;
    params.pre_move_thresh = 50;

    params.use_emg = 0;
    
    td_use = getMoveOnset(td_use,params);
  
%%

    neighbor_size = 20;
    stim_chan = 23;
    stim_code_use = 1;

    map_idx_stim_chan = find(map_data.chan == stim_chan);
    stim_chan_pos = [map_data.row(map_idx_stim_chan),map_data.col(map_idx_stim_chan)];
    
    dists = [];
    td_stim_code = td_use([td_use.stimCode] == stim_code_use);
    is_rewarded = [];
    num_spikes = [];
    
    dists = [];
    for chan = 1:size(td_stim_code(1).LeftS1_unit_guide,1)
        map_idx_rec_chan = find(map_data.chan == td_stim_code(1).LeftS1_unit_guide(chan,1));
        rec_chan_pos = [map_data.row(map_idx_rec_chan),map_data.col(map_idx_rec_chan)];
        dists(chan) = sqrt(sum((stim_chan_pos-rec_chan_pos).^2));
    end
    
    dist_mask = dists <= neighbor_size;
    
    
    for tr = 1:numel(td_stim_code)
        
        if(td_stim_code(tr).result == 'R')
            marker_use = '.';
            is_rewarded(tr) = 1;
        else
            marker_use = 'x';
            is_rewarded(tr) = 0;
        end
        window = td_stim_code(tr).idx_goCueTime + [0,200];
        num_spikes(tr,:) = (sum(td_stim_code(tr).LeftS1_spikes(window(1):window(2),dist_mask)));
        
    end

    %% build classifier and do leave one out cross validation
    % make data_table
    data_table = table(is_rewarded','VariableNames',{'is_rewarded'});
    
    for i = 1:size(num_spikes,2)
        data_table.(['var',num2str(i)]) = num_spikes(:,i); 
    end
    
    actual = is_rewarded; pred = [];
    for i = 1:size(data_table,1)
        % leave ith entry out
        data_use = data_table; data_use(i,:) = [];
        mdl = fitcsvm(data_use,'is_rewarded');
%         mdl = fitcdiscr(data_use,'is_rewarded');
        pred(i) = predict(mdl,data_table(i,2:end));
    end
    
    [sum(actual == pred)/numel(pred), sum(actual==1)/numel(actual)]
