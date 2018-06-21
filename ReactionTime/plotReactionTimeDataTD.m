function [outputData,plots] = plotReactionTimeDataTD(td_reward,td_all,opts)

    %% configure opts
    opts = configureOpts(opts);
    plots = [];
    f_means_bump = [];
    f_all_bump = [];
    f_all_stim = [];
    f_means_stim = [];
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
    
    %% remove fastest 5% and slowest 5%
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
        td_all_remove = unique([td_all_remove,reshape(td_all_remove + [1:opts.NUM_SEQUENCE_FAILS-1]',1,numel(td_all_remove)*(opts.NUM_SEQUENCE_FAILS-1))]);
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
        plots{end+1} = figure();
        plots{end}.Name = strcat(opts.FIGURE_PREFIX,'_bump_rtMeans');
        fitData.x = [];
        fitData.y = [];
        for cueIdx = 1:numel(cueInfo)
            if(cueInfo(cueIdx).bumpMag ~= 0 && ~isempty(cueInfo(cueIdx).rt))
                plot(cueInfo(cueIdx).bumpMag,mean(cueInfo(cueIdx).rt),'k.','markersize',opts.MARKER_SIZE)
                hold on
                plot(cueInfo(cueIdx).bumpMag + [0,0],mean(cueInfo(cueIdx).rt) + [-std(cueInfo(cueIdx).rt), std(cueInfo(cueIdx).rt)],'k')
                fitData.x(end+1,1) = cueInfo(cueIdx).bumpMag;
                fitData.y(end+1,1) = mean(cueInfo(cueIdx).rt);
            end
        end

        % if fit, fit with a decaying exponential
        f_means_bump = [];
        if(opts.FIT)
            [f_means_bump.fitObj,f_means_bump.gof] = fit(fitData.x,fitData.y,'a*exp(b*x)+c','startPoint',[0,0,0.15]);

            xData = linspace(min(fitData.x*0.9),max(fitData.x)*1.1,100);
            yData = f_means_bump.fitObj.a*exp(f_means_bump.fitObj.b*xData)+f_means_bump.fitObj.c;
            plot(xData,yData,'k--','linewidth',opts.LINE_WIDTH);
        end

        ax1 = gca;
        ax1.XLim(1) = 0;
        ax1.YLim(1) = 0;
        xlabel('Bump Magnitude (N)');
        ylabel('RT (s)');
        formatForLee(gcf);
        set(gca,'fontsize',opts.FONT_SIZE);
        % plot all rt's for bumps
        plots{end+1} = figure();
        plots{end}.Name = strcat(opts.FIGURE_PREFIX,'_bump_rtAllDots');
        hold on
        fitData.x = [];
        fitData.y = [];
        for cueIdx = 1:numel(cueInfo)
            if(cueInfo(cueIdx).bumpMag ~= 0 && ~isempty(cueInfo(cueIdx).rt))
                plot(cueInfo(cueIdx).bumpMag,cueInfo(cueIdx).rt,'k.','markersize',opts.MARKER_SIZE)
                fitData.x = [fitData.x,cueInfo(cueIdx).bumpMag+zeros(size(cueInfo(cueIdx).rt))];
                fitData.y = [fitData.y,cueInfo(cueIdx).rt];
            end
        end

        
        if(opts.FIT)
            [f_all_bump.fitObj,f_all_bump.gof] = fit(fitData.x',fitData.y','a*exp(b*x)+c','startPoint',[0,0,0.15]);

            xData = linspace(min(fitData.x*0.9),max(fitData.x)*1.1,100);
            yData = f_all_bump.fitObj.a*exp(f_all_bump.fitObj.b*xData)+f_all_bump.fitObj.c;
            plot(xData,yData,'k--','linewidth',opts.LINE_WIDTH);
        end

        ax2 = gca;
        ax2.XLim(1) = 0;
        ylim([0,max(ax1.YLim(2),ax2.YLim(2))]);
        xlabel('Bump Magnitude (N)')
        ylabel('RT (s)');
        formatForLee(gcf);
        ax1.YLim(2) = max(ax1.YLim(2),ax2.YLim(2));
        set(gca,'fontsize',opts.FONT_SIZE);

    end
%% plot mean rt as a function of stim code
    if(opts.PLOT_STIM)
        plots{end+1} = figure();
        plots{end}.Name = strcat(opts.FIGURE_PREFIX,'_stim_rtMeans');
        fitData.x = [];
        fitData.y = [];
        for cueIdx = 1:numel(cueInfo)
            if(cueInfo(cueIdx).stimCode ~= -1 && ~isempty(cueInfo(cueIdx).rt) && cueInfo(cueIdx).bumpMag == 0)
                plot(opts.STIM_PARAMS(cueInfo(cueIdx).stimCode+1),mean(cueInfo(cueIdx).rt),'k.','markersize',opts.MARKER_SIZE)
                hold on
                plot(opts.STIM_PARAMS(cueInfo(cueIdx).stimCode+1) + [0,0],mean(cueInfo(cueIdx).rt) + [-std(cueInfo(cueIdx).rt), std(cueInfo(cueIdx).rt)],'k')
                fitData.x(end+1,1) = opts.STIM_PARAMS(cueInfo(cueIdx).stimCode+1);
                fitData.y(end+1,1) = mean(cueInfo(cueIdx).rt);
            end
        end

        % if fit, fit with a decaying exponential
        f_means_stim = [];
        if(opts.FIT)
            [f_means_stim.fitObj,f_means_stim.gof] = fit(fitData.x,fitData.y,'a*exp(b*x)+c','startPoint',[0,0,0.15]);

            xData = linspace(min(fitData.x*0.9),max(fitData.x)*1.1,100);
            yData = f_means_stim.fitObj.a*exp(f_means_stim.fitObj.b*xData)+f_means_stim.fitObj.c;
            plot(xData,yData,'k--','linewidth',opts.LINE_WIDTH);
        end

        ax1 = gca;
        ax1.XLim(1) = 0;
        ax1.YLim(1) = 0;
        xlabel(opts.STIM_LABEL);
        ylabel('RT (s)');
        formatForLee(gcf);
        set(gca,'fontsize',opts.FONT_SIZE);

        % plot all rt's for bumps
        plots{end+1} = figure();
        plots{end}.Name = strcat(opts.FIGURE_PREFIX,'_stim_rtAllDots');
        fitData.x = [];
        fitData.y = [];
        for cueIdx = 1:numel(cueInfo)
            if(cueInfo(cueIdx).stimCode ~= -1 && ~isempty(cueInfo(cueIdx).rt) && cueInfo(cueIdx).bumpMag == 0)
                plot(opts.STIM_PARAMS(cueInfo(cueIdx).stimCode+1),cueInfo(cueIdx).rt,'k.','markersize',opts.MARKER_SIZE)
                hold on
                fitData.x = [fitData.x,opts.STIM_PARAMS(cueInfo(cueIdx).stimCode+1)+zeros(size(cueInfo(cueIdx).rt))];
                fitData.y = [fitData.y,cueInfo(cueIdx).rt];
            end
        end

        
        if(opts.FIT)
            [f_all_stim.fitObj,f_all_stim.gof] = fit(fitData.x',fitData.y','a*exp(b*x)+c','startPoint',[0,0,0.15]);

            xData = linspace(min(fitData.x*0.9),max(fitData.x)*1.1,100);
            yData = f_all_stim.fitObj.a*exp(f_all_stim.fitObj.b*xData)+f_all_stim.fitObj.c;
            plot(xData,yData,'k--','linewidth',opts.LINE_WIDTH);
        end

        ax2 = gca;
        ax2.XLim(1) = 0;
        ylim([0,max(ax1.YLim(2),ax2.YLim(2))]);
        xlabel(opts.STIM_LABEL)
        ylabel('RT (s)');
        formatForLee(gcf);
        set(gca,'fontsize',opts.FONT_SIZE);
        ax1.YLim(2) = max(ax1.YLim(2),ax2.YLim(2));
    end

    %% plot psychometric curve for bump data
    
    if(opts.PLOT_BUMP)
        plots{end+1} = figure();
        plots{end}.Name = strcat(opts.FIGURE_PREFIX,'_bump_detectionCurve');
        hold on
        fitData.x = [];
        fitData.y = [];
        for cueIdx = 1:numel(cueInfo)
            if(cueInfo(cueIdx).bumpMag ~= 0 && ~isempty(cueInfo(cueIdx).rt))
                plot(cueInfo(cueIdx).bumpMag,cueInfo(cueIdx).percent_respond,'k.','markersize',opts.MARKER_SIZE)
                fitData.x(end+1,1) = cueInfo(cueIdx).bumpMag;
                fitData.y(end+1,1) = cueInfo(cueIdx).percent_respond;
            end
        end

        % fit psychometric curve

        if(opts.FIT)
            xData = linspace(0,max(fitData.x*1.1),100);

            if(opts.USE_ML_FIT)
                mlFit_fun = @(params)(sum((fitData.y-(params(1)+params(2)*erf(params(3)*(fitData.x-params(4))))).^2));
                min_options = optimset('MaxIter',10000,'MaxFunEvals',10000);
                g_bump = fminsearch(mlFit_fun, [.45 .4 .05 1],min_options);

                yData = g_bump(1) + g_bump(2)*erf(g_bump(3)*(xData-g_bump(4))); 
            else
                s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', [.5 .5 .1 90], 'Lower', [0 0 0 0], 'Upper', [1 1 10 1.5],...
                    'MaxFunEvals',10000,'MaxIter',1000,'TolFun',10E-8,'TolX',10E-8);
                ft = fittype('a+b*(erf(c*(x-d)))','options',s);
                [g_bump.fitObj,g_bump.gof] = fit(fitData.x', fitData.y', ft);

                yData = feval(g_bump.fitObj,xData);
            end

            plot(xData,yData,'k--','linewidth',opts.LINE_WIDTH);
        end

        ax = gca;
        ax.XLim(1) = 0;
        ax.YLim = [0,1];
        xlabel('Bump Magnitude (N)');
        ylabel('Proportion detected');
        formatForLee(gcf);
        set(gca,'fontsize',opts.FONT_SIZE);

    end
    %% plot psychometric curve for stim data
    
    if(opts.PLOT_STIM)
        plots{end+1} = figure();
        plots{end}.Name = strcat(opts.FIGURE_PREFIX,'_stim_detectionCurve');
        hold on
        fitData.x = [];
        fitData.y = [];
        for cueIdx = 1:numel(cueInfo)
            if(cueInfo(cueIdx).stimCode ~= -1 && cueInfo(cueIdx).bumpMag == 0)
                plot(opts.STIM_PARAMS(cueInfo(cueIdx).stimCode+1),cueInfo(cueIdx).percent_respond,'k.','markersize',opts.MARKER_SIZE)
                fitData.x(end+1,1) = opts.STIM_PARAMS(cueInfo(cueIdx).stimCode+1);
                fitData.y(end+1,1) = cueInfo(cueIdx).percent_respond;
            end
        end

        % fit psychometric curve
        if(opts.FIT)
            xData = linspace(0,max(fitData.x*1.1),100);

            if(opts.USE_ML_FIT)
                mlFit_fun = @(params)(sum((fitData.y-(params(1)+params(2)*erf(params(3)*(fitData.x-params(4))))).^2));
                min_options = optimset('MaxIter',10000,'MaxFunEvals',10000);
                g_stim = fminsearch(mlFit_fun, [.45 .4 .05 1],min_options);

                yData = g_stim(1) + g_stim(2)*erf(g_stim(3)*(xData-g_stim(4))); 
            else
                s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', [.5 .5 .1 90], 'Lower', [0 0 0 0], 'Upper', [1 1 10 1.5],...
                    'MaxFunEvals',10000,'MaxIter',1000,'TolFun',10E-8,'TolX',10E-8);
                ft = fittype('a+b*(erf(c*(x-d)))','options',s);
                [g_stim.fitObj,g_stim.gof] = fit(fitData.x', fitData.y', ft);

                yData = feval(g_stim.fitObj,xData);
            end

            plot(xData,yData,'k--','linewidth',opts.LINE_WIDTH);
        end

        ax = gca;
        ax.XLim(1) = 0;
        ax.YLim = [0,1];
        xlabel(opts.STIM_LABEL);
        ylabel('Proportion detected');
        formatForLee(gcf);
        set(gca,'fontsize',opts.FONT_SIZE);

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
    outputData.rt_fit_means_stim = f_means_stim;
    outputData.rt_fit_all_stim = f_all_stim;
    outputData.rt_fit_means_bump = f_means_bump;
    outputData.rt_fit_all_bump = f_all_bump;
    outputData.psychometric_fit_stim = g_stim;
    outputData.psychometric_fit_bump = g_bump;
    outputData.cueInfo = cueInfo;
    outputData.stim_learn_fit = learn_mdl_stim;
    outputData.bump_learn_fit = learn_mdl_bump;
end


function [opts] = configureOpts(optsInput)

    opts = [];
    
    opts.MIN_BIN = 0.1;
    opts.MAX_BIN = 0.6;
    opts.BIN_SIZE = 0.02;
    
    opts.LINE_WIDTH = 1.5;
    opts.MARKER_SIZE = 12;
    opts.FIT = 1;
    opts.USE_ML_FIT = 1;
    opts.FONT_SIZE = 16;
    
    opts.COLORS = {'k','r',[0 0.6 0],'b','m',[0.5,0.5,0.5],[0.4,0.4,0.2]};

    opts.SAVE_FIGURES = 0;
    opts.FIGURE_PREFIX = '';
    opts.FOLDER_PATH = '';
    
    opts.BUMP_MAGS = [];
    opts.STIM_CODES = [];
    opts.STIM_PARAMS = [];
    
    opts.PLOT_BUMP = 1;
    opts.PLOT_STIM = 1;
    
    opts.NUM_SEQUENCE_FAILS = 3;
    opts.STIM_LABEL = 'Stim amp (\muA)';
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