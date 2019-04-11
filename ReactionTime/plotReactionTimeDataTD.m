function [outputData,plots] = plotReactionTimeDataTD(td_reward,td_all,opts)

    %% configure opts
    opts = configureOpts(opts);
    plots = [];
    f_bump = [];
    f_stim = [];
    g_bump = [];
    g_stim = [];
    learn_mdl_bump = [];
    learn_mdl_stim = [];
    %% get indexes for each cue type
    if(isempty(opts.BUMP_MAGS))
        bumpList = unique([td_reward.bumpMagnitude]);
        bumpList = bumpList(~isnan(bumpList)); % bump mag of 0 means no bump
    else
        bumpList = opts.BUMP_MAGS;
    end
    
    if(isempty(opts.STIM_CODES))
        stimCodeList = unique([td_reward.stimCode]);
        stimCodeList = stimCodeList(~isnan(stimCodeList)); 
        stimCodeList = 0:max(stimCodeList);
        stimCodeList(end+1) = -1; % representing no stim
    else
        stimCodeList = opts.STIM_CODES;
    end
    
    
    cueIdx = 1;
    cueInfo = [];
    
%     num_remove = floor(size(td_reward,2)*0.05);
%     rt_all = [td_reward.idx_movement_on]-[td_reward.idx_goCueTime];
%     [~,rt_sort_idx] = sort(rt_all);
%     rt_remove = [rt_sort_idx(1:num_remove),rt_sort_idx(end-num_remove+1:end)];
%     trial_id_remove = [td_reward(rt_remove).trial_id];
%     td_reward(rt_remove) = [];
%     for i = 1:numel(trial_id_remove)
%         td_all_idx = find([td_all.trial_id] == trial_id_remove(i));
%         td_all(td_all_idx) = [];
%     end
    %% remove chains of fail from td_all -- this is when the monkey is not paying attention
    fail_idx = find([td_all.result] == 'F');
    td_all_remove = fail_idx(find(fail_idx(opts.NUM_SEQUENCE_FAILS:end) - fail_idx(1:end-opts.NUM_SEQUENCE_FAILS+1) == opts.NUM_SEQUENCE_FAILS-1));
    if(~isempty(td_all_remove))
        td_all_remove = unique([td_all_remove;repmat(td_all_remove,opts.NUM_SEQUENCE_FAILS-1,1) + repmat([1:opts.NUM_SEQUENCE_FAILS-1]',1,size(td_all_remove,2))]);
        td_all(td_all_remove) = [];
    end
    
    %% for each cue, store indexes in td_reward, get reaction time
    for b = 1:numel(bumpList)
        for s = 1:numel(stimCodeList)
            % get idx in td_reward
            if(stimCodeList(s) == -1) % handle no stim, just a bump case
                cueInfo(cueIdx).td_idx_reward = find(isEqual(bumpList(b),[td_reward.bumpMagnitude]) & ~[td_reward.isStimTrial] & [td_reward.result]=='R');
                cueInfo(cueIdx).td_idx_all = find(isEqual(bumpList(b),[td_all.bumpMagnitude]) & ~[td_all.isStimTrial]);
            else % handle bump and stim (note, bump == 0 means no bump), so the stim only cases are covered here as well
                cueInfo(cueIdx).td_idx_reward = find(isEqual(bumpList(b),[td_reward.bumpMagnitude]) & isEqual(stimCodeList(s),[td_reward.stimCode]) & [td_reward.isStimTrial] & [td_reward.result]=='R');
                cueInfo(cueIdx).td_idx_all = find(isEqual(bumpList(b),[td_all.bumpMagnitude]) & isEqual(stimCodeList(s),[td_all.stimCode]) & [td_all.isStimTrial]);
            end
            
            % find reaction time
            cueInfo(cueIdx).rt = mode([td_reward.bin_size])*[(td_reward(cueInfo(cueIdx).td_idx_reward).idx_movement_on)] - [td_reward(cueInfo(cueIdx).td_idx_reward).goCueTime];
            
            % find % responsive
            cueInfo(cueIdx).percent_respond = numel(cueInfo(cueIdx).td_idx_reward)/(numel(cueInfo(cueIdx).td_idx_all));
            if(numel(cueInfo(cueIdx).td_idx_all) > 0 && isnan(cueInfo(cueIdx).percent_respond))
                cueInfo(cueIdx).percent_respond = 0;
            end
            % store cue info
            cueInfo(cueIdx).bumpMag = bumpList(b);
            cueInfo(cueIdx).stimCode = stimCodeList(s);
            
            cueIdx = cueIdx + 1;

        end
    end
        
    %% remove outliers
        % anything above 1.5 Inter quartile range goes bye bye
    for c = 1:numel(cueInfo)
        Q1 = quantile(cueInfo(c).rt,0.25);
        Q3 = quantile(cueInfo(c).rt,0.75);
        outlier_mask = cueInfo(c).rt > Q3+1.5*(Q3-Q1);

%         cueInfo(c).td_idx_reward(outlier_mask) = [];
%         cueInfo(c).td_idx_all(outlier_mask) = [];
        cueInfo(c).rt(outlier_mask) = [];
    end
        
    %% find fastest bump data (across bump magnitudes)
    fastest_bump_data.mean = 10000;
    fastest_bump_data.std_err = [];
    
    
    for cueIdx = 1:numel(cueInfo)
        if(cueInfo(cueIdx).bumpMag~=0 && cueInfo(cueIdx).stimCode == -1)
            if(mean(cueInfo(cueIdx).rt) < fastest_bump_data.mean && numel(cueInfo(cueIdx).rt)>3)
                fastest_bump_data.mean = mean(cueInfo(cueIdx).rt);
                fastest_bump_data.std_err = std(cueInfo(cueIdx).rt)/sqrt(numel(cueInfo(cueIdx).rt));
            end
        end
    end
        
    %% get visual data (if it exists)
    visual_data = [];
    cueInfo_idx = find([cueInfo.bumpMag] == 0 & [cueInfo.stimCode] == -1);
    if(~isempty(cueInfo_idx) && ~isempty(cueInfo(cueInfo_idx).rt))
        visual_data.mean = mean(cueInfo(cueInfo_idx).rt);
        visual_data.std_err = std(cueInfo(cueInfo_idx).rt)/sqrt(numel(cueInfo(cueInfo_idx).rt));
    end
    
    %% make a histogram
%     % get bin counts
%     bE = opts.MIN_BIN:opts.BIN_SIZE:opts.MAX_BIN;
%     %  each cue
%     for cueIdx = 1:numel(cueInfo)
%         bC(cueIdx,:) = histcounts(cueInfo(cueIdx).rt,bE)/numel(cueInfo(cueIdx).rt);
%         cueIdx = cueIdx + 1;
%     end
%     % make plot
%     plots{end+1} = figure();
%     plots{end}.Name = strcat(opts.FIGURE_PREFIX,'_rtHistogram');
%     hold on
%     for cueIdx = 1:size(bC,1)
%         plot(bE(1:end-1)+mode(diff(bE))/2,bC(cueIdx,:),'-','color',opts.COLORS{cueIdx},'linewidth',opts.LINE_WIDTH)
%     end
%     l=legend(num2str([cueInfo.bumpMag]'));
%     set(l,'box','off','location','best');
%     xlabel('RT (s)');
%     ylabel('Proportion of trials');
%     formatForLee(gcf);
    
    %% plot mean rt as a function of bump magnitude
    
    if(opts.PLOT_BUMP)
        
        % extract data
        plot_data.x = [];
        plot_data.y = [];
        
        if(~isempty(visual_data))
            plot_data.visual_data = visual_data;
        end
        
        for cueIdx = 1:numel(cueInfo)
            if(cueInfo(cueIdx).bumpMag ~= 0 && ~isempty(cueInfo(cueIdx).rt))
                plot_data.x = [plot_data.x;ones(numel(cueInfo(cueIdx).rt),1)*cueInfo(cueIdx).bumpMag];
                plot_data.y = [plot_data.y;cueInfo(cueIdx).rt'];
            end
        end
        % call plotting function
        figure_name = strcat(opts.FIGURE_PREFIX,'_bump_rtMeans');
        bump_flag = 1;
        [plots{end+1}, f_bump] = plotRTvsCue(plot_data,figure_name,bump_flag,opts);

    end
%% plot mean rt as a function of stim code
    if(opts.PLOT_STIM)
        
        % extract data
        plot_data.x = [];
        plot_data.y = [];
        
        plot_data.bump_data = fastest_bump_data;
        if(~isempty(visual_data))
            plot_data.visual_data = visual_data;
        end
        
        stim_param_idx = 1;
        for cueIdx = 1:numel(cueInfo)
            if(cueInfo(cueIdx).stimCode ~= -1 && (isempty(opts.STIM_CODES) || ~isempty(find(opts.STIM_CODES == cueInfo(cueIdx).stimCode)))...
                    && cueInfo(cueIdx).bumpMag == 0)
                plot_data.x = [plot_data.x;ones(numel(cueInfo(cueIdx).rt),1)*opts.STIM_PARAMS(stim_param_idx)];
                plot_data.y = [plot_data.y;cueInfo(cueIdx).rt'];
                stim_param_idx = stim_param_idx + 1;                
            end
        end
        
        % call plotting function
        figure_name = strcat(opts.FIGURE_PREFIX,'_stim_rtMeans');
        bump_flag = 0;
        [plots{end+1}, f_stim] = plotRTvsCue(plot_data,figure_name,bump_flag,opts);
    end

    %% plot psychometric curve for bump data
    
    if(opts.PLOT_BUMP)
        % extract data
        plot_data.x = [];
        plot_data.y = [];
        
        for cueIdx = 1:numel(cueInfo)
            if(cueInfo(cueIdx).bumpMag ~= 0 && ~isempty(cueInfo(cueIdx).rt))
                plot_data.x = [plot_data.x;cueInfo(cueIdx).bumpMag];
                plot_data.y = [plot_data.y;cueInfo(cueIdx).percent_respond];
            end
        end
        % call plotting function
        figure_name = strcat(opts.FIGURE_PREFIX,'_bump_detectionCurve');
        bump_flag = 1;
        [plots{end+1}, f_bump] = plotDetectionCurve(plot_data,figure_name,bump_flag,opts);
        
    end
    %% plot psychometric curve for stim data
    
    if(opts.PLOT_STIM)
        % extract data
        plot_data.x = [];
        plot_data.y = [];
        
        stim_param_idx = 1;
        for cueIdx = 1:numel(cueInfo)
            if(cueInfo(cueIdx).stimCode ~= -1 && cueInfo(cueIdx).bumpMag == 0 && ~isnan(cueInfo(cueIdx).percent_respond))
                plot_data.x = [plot_data.x;opts.STIM_PARAMS(stim_param_idx)];
                plot_data.y = [plot_data.y;cueInfo(cueIdx).percent_respond];
                stim_param_idx = stim_param_idx + 1;
            end
        end
        % call plotting function
        figure_name = strcat(opts.FIGURE_PREFIX,'_stim_detectionCurve');
        bump_flag = 0;
        [plots{end+1}, f_bump] = plotDetectionCurve(plot_data,figure_name,bump_flag,opts);
       
        
    end
    %% plot BUMP reaction time vs trial idx
    
    if(opts.PLOT_BUMP)
        plots{end+1} = figure();
        plots{end}.Name = strcat(opts.FIGURE_PREFIX,'_bump_learning');
        if(isfield(td_reward,'isBumpTrial'))
            mask = [td_reward.isBumpTrial] == 1;
        else
            mask = [td_reward.bumpMagnitude]~=0;
        end
        td_rt = td_reward(mask);
        rt = [td_rt.idx_movement_on] - [td_rt.idx_goCueTime];
        rt = rt*mode([td_rt.bin_size]);
        plot(rt,'k','linewidth',1.5)
        xlabel('Trial')
        ylabel('Reaction time (s)');
        ax = gca;
        ax.YLim(1) = 0;
        if(opts.FIT)
            rt_learn = [];
            xFit = 1:numel(rt);
            yFit = rt;
            tbl = table(xFit',yFit','VariableNames',{'trial','rt'});
            learn_mdl_bump = fitlm(tbl,'rt ~ trial');
%             [rt_learn.fitObj, rt_learn.gof] = fit(xFit',yFit','a*x+b');
            hold on
            plot(xFit, learn_mdl_bump.Coefficients.Estimate(2)*xFit+learn_mdl_bump.Coefficients.Estimate(1),'r--','linewidth',1.5);
        end
        formatForLee(gcf);
        set(gca,'fontsize',opts.FONT_SIZE);

    end
    
    %% plot STIM reaction time vs trial idx
    if(opts.PLOT_STIM)
        plots{end+1} = figure();
        plots{end}.Name = strcat(opts.FIGURE_PREFIX,'_stim_learning');
        td_rt = td_reward([td_reward.isStimTrial]==1);
        rt = [td_rt.idx_movement_on] - [td_rt.idx_goCueTime];
        rt = rt*mode([td_rt.bin_size]);
        plot(rt,'k','linewidth',1.5)
        xlabel('Trial')
        ylabel('Reaction time (s)');
        ax = gca;
        ax.YLim(1) = 0;
        if(opts.FIT)
            rt_learn = [];
            xFit = 1:numel(rt);
            yFit = rt;
            tbl = table(xFit',yFit','VariableNames',{'trial','rt'});
            learn_mdl_stim = fitlm(tbl,'rt ~ trial');
%             [rt_learn.fitObj, rt_learn.gof] = fit(xFit',yFit','a*x+b');
            hold on
            plot(xFit, learn_mdl_stim.Coefficients.Estimate(2)*xFit+learn_mdl_stim.Coefficients.Estimate(1),'r--','linewidth',1.5);
        end
        formatForLee(gcf);
        set(gca,'fontsize',opts.FONT_SIZE);

    end
    
    %% plot reaction time as a function of hold period
    plots{end+1} = figure();
    plots{end}.Name = strcat(opts.FIGURE_PREFIX,'_holdPeriod_rt');
    hold_period = ([td_rt.idx_goCueTime] - [td_rt.idx_tgtOnTime])*mode([td_rt.bin_size]);
    rt = ([td_rt.idx_movement_on] - [td_rt.idx_goCueTime])*mode([td_rt.bin_size]);
    plot(hold_period,rt,'k.','markersize',opts.MARKER_SIZE)
    xlabel('Hold period (s)')
    ylabel('RT (s)')
    formatForLee(gcf)
    set(gca,'fontsize',opts.FONT_SIZE);
    
    %% deal with saving figures
    if(opts.SAVE_FIGURES && strcmp(opts.FOLDER_PATH,'')==0)
        for p = 1:numel(plots)
            saveFiguresLIB(plots{p},opts.FOLDER_PATH,plots{p}.Name);
        end
    end
    %% setup output data
    outputData.rt_fit_stim = f_stim;
    outputData.rt_fit_bump = f_bump;
    outputData.psychometric_fit_stim = g_stim;
    outputData.psychometric_fit_bump = g_bump;
    outputData.cueInfo = cueInfo;
    outputData.stim_learn_fit = learn_mdl_stim;
    outputData.bump_learn_fit = learn_mdl_bump;
end



function [fig,fit_out] = plotRTvsCue(plot_data,figure_name,bump_flag,opts)

    fig = figure();
    fig.Name = figure_name;
    hold on;
    % plot mean and std error
    unique_x_vals = unique(plot_data.x);
    mean_y = zeros(size(unique_x_vals));
    std_err_y = zeros(size(unique_x_vals));
    
    % if fit, fit with a decaying exponential
    fit_out = [];
    if(opts.FIT)
        [fit_out.fitObj,fit_out.gof] = fit(plot_data.x,plot_data.y,'a*exp(b*x)+c','startPoint',[0,0,0]);

        x_data_fit = linspace(min(plot_data.x*0.9),max(plot_data.x)*1.1,100);
        y_data_fit = fit_out.fitObj.a*exp(fit_out.fitObj.b*x_data_fit)+fit_out.fitObj.c;
        plot(x_data_fit,y_data_fit,'k--','linewidth',opts.LINE_WIDTH,'color',opts.COLOR);
    end

    ax1 = gca;

    if(isfield(plot_data,'bump_data') && plot_data.bump_data.mean < 100 && ~bump_flag)
        bump_bar = fill([min(plot_data.x)*0.5,max(plot_data.x)*1.1,max(plot_data.x)*1.1,min(plot_data.x)*0.5],...
            plot_data.bump_data.mean(1)+plot_data.bump_data.std_err*[-1,-1,1,1],opts.COLOR,'linestyle','none');
        alpha(opts.ALPHA);
        plot([min(plot_data.x)*0.5,max(plot_data.x)*1.1],plot_data.bump_data.mean+[0,0],'-','linewidth',1,'color','k');
        uistack(bump_bar,'bottom');
    end
    
    if(isfield(plot_data,'visual_data') && plot_data.visual_data.mean < 100 && ~bump_flag)
        visual_bar = fill([min(plot_data.x)*0.5,max(plot_data.x)*1.1,max(plot_data.x)*1.1,min(plot_data.x)*0.5],...
            plot_data.visual_data.mean(1)+plot_data.visual_data.std_err*[-1,-1,1,1],opts.COLOR,'linestyle','none');
        alpha(opts.ALPHA);
        plot([min(plot_data.x)*0.5,max(plot_data.x)*1.1],plot_data.visual_data.mean+[0,0],'--','linewidth',1,'color','k');
        uistack(visual_bar,'bottom');
    end
    
    x_diff = max(unique_x_vals)-min(unique_x_vals);
    for i = 1:numel(unique_x_vals)
        if(~isempty(opts.STIM_COLOR_IDX))
            color = getColorFromList(opts.COLOR_LIST,opts.STIM_COLOR_IDX(i));
            if(~isempty(opts.STIM_COLOR_ALPHA))
                color = ((1-opts.STIM_COLOR_ALPHA(i))*1+(opts.STIM_COLOR_ALPHA(i)*color));
            end
        else
            color = opts.COLOR;
        end
        plot_data_mask = plot_data.x == unique_x_vals(i);
        mean_y(i) = mean(plot_data.y(plot_data_mask));
        std_err_y(i) = std(plot_data.y(plot_data_mask))/sqrt(sum(plot_data_mask));

        % plot std error
        errorbar(unique_x_vals(i),mean_y(i),std_err_y(i),...
            'marker','.','markersize',opts.MARKER_SIZE,...
            'linewidth',opts.LINE_WIDTH,'color',color);
%         errorbar_LIB(unique_x_vals(i),mean_y(i),std_err_y(i),...
%             'marker','.','markersize',opts.MARKER_SIZE,...
%             'linewidth',opts.LINE_WIDTH,'color',color,'caplength',opts.CAPLENGTH_MULTIPLIER*x_diff);
            
        % plot means
%         plot(unique_x_vals(i),mean_y(i),'.','color',color,'markersize',opts.MARKER_SIZE)
%         plot(unique_x_vals(i)+[0,0],mean_y(i)+[-1,1]*std_err_y(i),'-','color',opts.COLOR,'linewidth',opts.LINE_WIDTH);
        % plot horizontal bar on std error line
%         plot(unique_x_vals(i)+opts.HORZ_BAR_LENGTH*[-1,1],mean_y(i)+std_err_y(i)+[0,0],'-','color',opts.COLOR,'linewidth',opts.LINE_WIDTH);
%         plot(unique_x_vals(i)+opts.HORZ_BAR_LENGTH*[-1,1],mean_y(i)-std_err_y(i)+[0,0],'-','color',opts.COLOR,'linewidth',opts.LINE_WIDTH);

        % plot all dots
        plot_data_idx = find(plot_data.x == unique_x_vals(i));
        scatter(plot_data.x(plot_data_idx),plot_data.y(plot_data_idx),opts.MARKER_AREA,'markeredgecolor','none','markerfacecolor',color);
        alpha(opts.ALPHA);
    end
    
%     % plot all dots
%     scatter(plot_data.x,plot_data.y,opts.MARKER_AREA,'markeredgecolor','none','markerfacecolor',opts.COLOR);
%     alpha(opts.ALPHA);
    
    
    if(bump_flag)
        xlabel('Bump Magnitude (N)');
    else
        xlabel(opts.STIM_LABEL,'Interpreter','tex');
        if(~isempty(opts.STIM_X_LABEL))
            [~,sort_idx] = sort(opts.STIM_PARAMS);
            set(gca,'XTick',opts.STIM_PARAMS(sort_idx),'XTickLabel',opts.STIM_X_LABEL(sort_idx));
            set(gca,'XMinorTick','off');
        end
    end
    ylabel('RT (s)');
    formatForLee(gcf);
    set(gca,'fontsize',opts.FONT_SIZE);

end

function [fig, psych_fit] = plotDetectionCurve(plot_data,figure_name,bump_flag,opts)
    psych_fit = [];
    fig = figure();
    fig.Name = figure_name;
    hold on

    plot(plot_data.x,plot_data.y,'.','markersize',opts.MARKER_SIZE,'color',opts.COLOR);
   
    % fit psychometric curve
    if(opts.FIT)
        x_data = linspace(0,max(plot_data.x*1.1),100);

        if(opts.USE_ML_FIT)
            mlFit_fun = @(params)(sum((plot_data.y-(params(1)+params(2)*erf(params(3)*(plot_data.x-params(4))))).^2));
            min_options = optimset('MaxIter',10000,'MaxFunEvals',10000);
            psych_fit = fminsearch(mlFit_fun, [.45 .4 .05 1],min_options);

            y_data = psych_fit(1) + psych_fit(2)*erf(psych_fit(3)*(x_data-psych_fit(4))); 
        else
            s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', [.5 .5 .1 90], 'Lower', [0 0 0 0], 'Upper', [1 1 10 1.5],...
                'MaxFunEvals',10000,'MaxIter',1000,'TolFun',10E-8,'TolX',10E-8);
            ft = fittype('a+b*(erf(c*(x-d)))','options',s);
            [psych_fit.fitObj,psych_fit.gof] = fit(plot_data.x', plot_data.y', ft);

            y_data = feval(psych_fit.fitObj,x_data);
        end

        plot(x_data,y_data,'k--','linewidth',opts.LINE_WIDTH);
    end

    ax = gca;

    if(bump_flag)
        xlabel('Bump Magnitude (N)');
    else
        xlabel(opts.STIM_LABEL);
        if(~isempty(opts.STIM_X_LABEL))
            [~,sort_idx] = sort(opts.STIM_PARAMS);
            set(gca,'XTick',opts.STIM_PARAMS(sort_idx),'XTickLabel',opts.STIM_X_LABEL(sort_idx));
            set(gca,'XMinorTick','off');
        end
    end
    ylabel('Proportion detected');
    formatForLee(gcf);
    set(gca,'fontsize',opts.FONT_SIZE);

    
end


function [opts] = configureOpts(optsInput)

    opts = [];
    
    opts.MIN_BIN = 0.1;
    opts.MAX_BIN = 0.6;
    opts.BIN_SIZE = 0.02;
    
    opts.LINE_WIDTH = 1.5;
    opts.MARKER_SIZE = 24; % mean data points ...
    opts.MARKER_AREA = 20; % opaque raw data...
    opts.FIT = 1;
    opts.USE_ML_FIT = 1;
    opts.FONT_SIZE = 16;
    opts.ALPHA = 0.5;
    opts.HORZ_BAR_LENGTH = 0.; % horz bar on std error bars
    
    opts.COLORS = [228,26,28;55,126,184;77,175,74;152,78,163;255,127,0;255,255,51;166,86,40;247,129,191;153,153,153]/255;
    opts.COLOR = 'k';
    opts.COLOR_LIST = 1;
    
    opts.SAVE_FIGURES = 0;
    opts.FIGURE_PREFIX = '';
    opts.FOLDER_PATH = '';
    
    opts.BUMP_MAGS = [];
    opts.STIM_CODES = [];
    opts.STIM_PARAMS = [];
    opts.STIM_X_LABEL = {};
    
    opts.PLOT_BUMP = 1;
    opts.PLOT_STIM = 1;
    
    opts.NUM_SEQUENCE_FAILS = 3;
    opts.STIM_LABEL = 'Stim amp (\muA)';
    opts.STIM_COLOR_IDX = [];
    opts.STIM_COLOR_ALPHA = [];
    
    opts.CAPLENGTH_MULTIPLIER = 0.01;
    %% check if in optsSave and optsSaveInput, overwrite if so
    try
        inputFieldNames = fieldnames(optsInput);
        for fn = 1:numel(inputFieldNames)
           if(isfield(opts,inputFieldNames{fn}))
               opts.(inputFieldNames{fn}) = optsInput.(inputFieldNames{fn});
           end
        end
    catch
        % do nothing, [] was inputted which means use default setting
        error('could not parse opts');
    end
    

end