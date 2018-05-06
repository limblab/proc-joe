function [outputData,plots] = plotReactionTimeDataTD(td_reward,td_all,opts)

    %% configure opts
    opts = configureOpts(opts);
    
    
    %% get indexes for each cue type
    bumpList = unique([td_reward.bumpMagnitude]);
    bumpList = bumpList(~isnan(bumpList)); % bump mag of 0 means no bump
    
    stimCodeList = unique([td_reward.stimCode]);
    stimCodeList = stimCodeList(~isnan(stimCodeList)); 
    stimCodeList(end+1) = -1; % representing no stim
    
    cueIdx = 1;
    cueInfo = [];
    
    
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
            
            % store cue info
            cueInfo(cueIdx).bumpMag = bumpList(b);
            cueInfo(cueIdx).stimCode = stimCodeList(s);
            
            
            cueIdx = cueIdx + 1;

        end
    end
    

    %% make a histogram
    % get bin counts
    bE = opts.MIN_BIN:opts.BIN_SIZE:opts.MAX_BIN;
    %  each cue
    for cueIdx = 1:numel(cueInfo)
        bC(cueIdx,:) = histcounts(cueInfo(cueIdx).rt,bE)/numel(cueInfo(cueIdx).rt);
        cueIdx = cueIdx + 1;
    end
    % make plot
    plots = figure();
    hold on
    for cueIdx = 1:size(bC,1)
        plot(bE(1:end-1)+mode(diff(bE))/2,bC(cueIdx,:),'-','color',opts.COLORS{cueIdx},'linewidth',opts.LINE_WIDTH)
    end
    
    %% plot mean rt as a function of bump magnitude
    figure();
    
    fitData.x = [];
    fitData.y = [];
    for cueIdx = 1:numel(cueInfo)
        if(cueInfo(cueIdx).bumpMag ~= 0)
            plot(cueInfo(cueIdx).bumpMag,mean(cueInfo(cueIdx).rt),'k.','markersize',opts.MARKER_SIZE)
            hold on
            plot(cueInfo(cueIdx).bumpMag + [0,0],mean(cueInfo(cueIdx).rt) + [-std(cueInfo(cueIdx).rt), std(cueInfo(cueIdx).rt)],'k')
            fitData.x(end+1,1) = cueInfo(cueIdx).bumpMag;
            fitData.y(end+1,1) = mean(cueInfo(cueIdx).rt);
        end
    end
    
    % if fit, fit with a decaying exponential
    if(opts.FIT)
        f = [];
        [f.fitObj,f.gof] = fit(fitData.x,fitData.y,'a*exp(b*x)+c','startPoint',[0,0,0.15]);
    
        xData = linspace(min(fitData.x*0.9),max(fitData.x)*1.1,100);
        yData = f.fitObj.a*exp(f.fitObj.b*xData)+f.fitObj.c;
        plot(xData,yData,'k--','linewidth',opts.LINE_WIDTH);
    end
    
    xlim([0,max(xData)]);
    ylim([0,max(yData)*1.1]);
    
    % plot all rt's
    figure();
    hold on
    for cueIdx = 1:numel(cueInfo)
        if(cueInfo(cueIdx).bumpMag ~= 0)
            plot(cueInfo(cueIdx).bumpMag,cueInfo(cueIdx).rt,'k.','markersize',opts.MARKER_SIZE)
        end
    end
    
    ax = gca;
    xlim([0,max(fitData.x)]);
    ylim([0,ax.YLim(2)]);
    
    %% plot rt as a function of stim code
%     
%     

        

    %% plot psychometric curve for bump data
    figure();
    hold on
    fitData.x = [];
    fitData.y = [];
    for cueIdx = 1:numel(cueInfo)
        if(cueInfo(cueIdx).bumpMag ~= 0)
            plot(cueInfo(cueIdx).bumpMag,cueInfo(cueIdx).percent_respond,'k.','markersize',opts.MARKER_SIZE)
            fitData.x(end+1,1) = cueInfo(cueIdx).bumpMag;
            fitData.y(end+1,1) = cueInfo(cueIdx).percent_respond;
        end
    end
    
    % fit psychometric curve
    
    g = [];
    xData = linspace(0,max(fitData.x*1.1),100);

    if(opts.USE_ML_FIT)
        mlFit_fun = @(params)(sum((fitData.y-(params(1)+params(2)*erf(params(3)*(fitData.x-params(4))))).^2));
        min_options = optimset('MaxIter',10000,'MaxFunEvals',10000);
        g = fminsearch(mlFit_fun, [.45 .4 .05 1],min_options);

        yData = g(1) + g(2)*erf(g(3)*(xData-g(4))); 
    else
        s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', [.5 .5 .1 90], 'Lower', [0 0 0 0], 'Upper', [1 1 10 1.5],...
            'MaxFunEvals',10000,'MaxIter',1000,'TolFun',10E-8,'TolX',10E-8);
        ft = fittype('a+b*(erf(c*(x-d)))','options',s);
        [g.fitObj,g.gof] = fit(fitData.x', fitData.y', ft);

        yData = feval(g.fitObj,xData);
    end
    
    plot(xData,yData,'k--','linewidth',opts.LINE_WIDTH);
    
    xlim([0,max(fitData.x)*1.1]);
    ylim([0,1])
    
    %% setup output data
    outputData.fit = f.fitObj;
    outputData.gof = f.gof;
    
end


function [opts] = configureOpts(optsInput)

    opts = [];
    
    opts.MIN_BIN = 0.1;
    opts.MAX_BIN = 0.6;
    opts.BIN_SIZE = 0.025;
    
    opts.LINE_WIDTH = 1.5;
    opts.MARKER_SIZE = 12;
    opts.FIT = 1;
    opts.USE_ML_FIT = 1;
    
    opts.COLORS = {'r',[0 0.5 0],'b','k','m',[0.5,0.5,0.2],[0.3,0.3,0.3]};

    %% check if in optsSave and optsSaveInput, overwrite if so
    try
        inputFieldnames = fieldnames(optsInput);
        for fn = 1:numel(inputFieldnames)
           if(isfield(opts,inputFieldnames{fn}))
               opts.(inputFieldnames{fn}) = optsInput.(inputFieldnames{fn});
           end
        end
    catch
        % do nothing, [] was inputted which means use default setting
    end
    

end