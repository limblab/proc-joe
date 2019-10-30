%% choice direction stuff
%% start with simple neurometric curves -- plot firing rate as function of bump direction for each axis
    input_data.make_plot = 0;
    input_data.window = [0,0.301]; % s, will round down even if whole number is correct.
    neurometric_data = makeNeurometricCurves(td_all_use,psych_data,input_data);
    
    
%% see if distribution of FR for single neurons, or the population, is different for each choice
    % uses data in neurometric_data
    make_neuron_plots = 0; % plots FR distribution for each neuron in first window
    input_data.make_plot = 0; 
    bump_dirs_use = [85];
    
    windows = [0,0.03; 0,0.05; 0,0.08; 0,0.1; 0,0.2; 0,0.3];
    bin_edges = [0:5:200];
    axis_idx = 1;
    AUC_data = zeros(size(neurometric_data{axis_idx}.spike_data,2),1);
    
    for window_idx = 1:size(windows,1)
        input_data.window = windows(window_idx,:); % s
        neurometric_data = makeNeurometricCurves(td_all_use,psych_data,input_data);
        
        for nn = 1:size(neurometric_data{axis_idx}.spike_data,2)
            if(window_idx == 1 && make_neuron_plots == 1)
                figure();
                plot(bin_edges(1:end-1)+mode(diff(bin_edges))/2,bin_counts_0_deg,'k','linewidth',1.5);
                hold on
                plot(bin_edges(1:end-1)+mode(diff(bin_edges))/2,bin_counts_180_deg,'r','linewidth',1.5);
            end
            
            % only use bump dirs in bump_dirs_use
            bump_dir_mask = sum(neurometric_data{axis_idx}.bump_dirs == bump_dirs_use,2) > 0;
            % bin data based on choice direction
            FR_data = sum(neurometric_data{axis_idx}.spike_data(bump_dir_mask,nn,:),3)/diff(input_data.window);
            bin_counts_0_deg = histcounts(FR_data(neurometric_data{axis_idx}.reached_0_deg(bump_dir_mask) == 1),bin_edges,'normalization','probability');
            bin_counts_180_deg = histcounts(FR_data(neurometric_data{axis_idx}.reached_0_deg(bump_dir_mask) == 0),bin_edges,'normalization','probability');

            [~,~,~,AUC_data(nn,window_idx)] = perfcurve(neurometric_data{axis_idx}.reached_0_deg(bump_dir_mask),FR_data,1);
        end
    end
    
    figure();
    plot(windows(:,2),AUC_data')
    
%% get PD around bumps
    time_post_bump = 0.4; % s
    td_bump = trimTD(td_all_use,{'idx_bumpTime',0},{'idx_bumpTime',floor(time_post_bump/td_all_use(1).bin_size)});
    
    pd_params = [];
    pd_params.out_signals = 'LeftS1_spikes';
    pd_params.bootForTuning = 0;
    pd_params.num_boots = 1;
    pd_params.move_corr = 'vel';
    
    pd_bump = getTDPDs(td_bump,pd_params);
    
%% compare AUC result with bump PD
    plot(pd_bump.velPD,AUC_data(:,end),'.')

%% predict correct target based on neural data and compare to monkey's behavior
    num_pca_dims = 6;
    axis_idx = 1;
    time_post_bump = 0.4; % s

    unique_target_dirs = unique([td_all_use.target_direction]);
    params_pca = [];
    params_pca.algorithm = 'pca';
    params_pca.signals = [input_data.array1(6:end),'_spikes'];
    params_pca.num_dims = num_pca_dims;
    
    
    td_bump = trimTD(td_all_use,{'idx_bumpTime',0},{'idx_bumpTime',floor(time_post_bump/td_all_use(1).bin_size)});

    % get only target direction in the specified direction
    td_bump = td_bump([td_bump.target_direction] == unique_target_dirs(axis_idx));
        
    [td_bump,info_out] = dimReduce(td_bump,params_pca);
    
    % put input and output signals together
    neural_pca_data = zeros(numel(td_bump),numel(td_bump(1).LeftS1_pca));
    kin_data = zeros(numel(td_bump),size(td_bump(1).vel,1)*3*2); % 3 kin variables, 2 directions
    
    for tr = 1:numel(td_bump)
        neural_pca_data(tr,:) = reshape(td_bump(tr).LeftS1_pca,1,numel(td_bump(tr).LeftS1_pca));
        kin_data(tr,:) = [reshape(td_bump(tr).pos,1,numel(td_bump(tr).pos)), ...
            reshape(td_bump(tr).vel,1,numel(td_bump(tr).vel)),...
            reshape(td_bump(tr).acc,1,numel(td_bump(tr).acc))];
    end
    target_dir = ([td_bump.bumpDir] < 90 | [td_bump.bumpDir] > 270)';
    choice_dir = neurometric_data{axis_idx}.reached_0_deg;
    
%% get LDA and do cross validation
    Y = target_dir;
    X = neural_pca_data;
    
    gamma = 0.75; % regularization term
    num_cross_val = 10;

    indices = crossvalind('Kfold',Y,num_cross_val);
    lda_predictions = zeros(size(Y));
    lda_training_acc = zeros(num_cross_val,1);
    
    for cv_idx = 1:num_cross_val
        test_mask = indices == cv_idx;
        
        x_fit = X(test_mask == 0,:); x_test = X(test_mask == 1,:);
        y_fit = Y(test_mask == 0); y_test = Y(test_mask == 1);
        
        lda = fitcdiscr(x_fit,y_fit,'Gamma',gamma);
        lda_cv_pred = lda.predict(x_test);
        
        lda_predictions(test_mask) = lda_cv_pred;
        lda_training_acc(cv_idx) = sum(lda.predict(x_fit)==y_fit)/numel(y_fit);
        
    end
    
    sum(lda_predictions == Y)/numel(Y)
    
%% plot prediction accuracy for each bump direction
    
    bump_dirs = unique(neurometric_data{axis_idx}.bump_dirs);
    monkey_acc = zeros(size(bump_dirs));
    lda_acc = zeros(size(bump_dirs));
    
    for bump_idx = 1:numel(bump_dirs)
        bump_mask = [td_bump.bumpDir] == bump_dirs(bump_idx);
        
        % get monkey_acc
        monkey_acc(bump_idx) = sum(neurometric_data{axis_idx}.reached_0_deg(bump_mask) == 1)/sum(bump_mask);
        if(bump_dirs(bump_idx) > 90)
            monkey_acc(bump_idx) = 1 - monkey_acc(bump_idx);
        end
        
        % get lda_acc
        lda_acc(bump_idx) = sum(lda_predictions(bump_mask) == Y(bump_mask))/sum(bump_mask);
        
    end

    figure();
    plot(bump_dirs,monkey_acc,'k');
    hold on
    plot(bump_dirs,lda_acc,'r');
    







%% perform analysis done by (Ingaki, 2019) -- coding direction
% find n-dimensional vector (n = num neurons) where each entry is the
% average difference between spikes during correct left and right trials

    input_data.sample_rate = 0.5;
    input_data.bump_range = 90 + 50*[-1,1];
    for axis = 1:numel(neurometric_data)
        cd_data{axis} = getCodingDirection(neurometric_data{axis},input_data);
    end
    
%% predict decision and visualize projection onto cd over time
    input_data.bump_range = 90 + 90*[-1,1];

    for axis = 1:numel(neurometric_data)
        cd_pred_data{axis} = getCodingDirectionPredictions(neurometric_data{axis},cd_data{axis},input_data);
    end
    


