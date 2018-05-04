function [outputData,plots] = plotReactionTimeDataTD(td,opts)

    %% configure opts
    opts = configureOpts(opts);
    
    
    %% get indexes for each cue type
    bumpList = unique([td.bumpMagnitude]);
    bumpList = bumpList(~isnan(bumpList)); % bump mag of 0 means no bump
    
    stimCodeList = unique([td.stimCode]);
    stimCodeList = stimCodeList(~isnan(stimCodeList)); 
    stimCodeList(end+1) = -1; % representing no stim
    
    cueIdx = 1;
    cueInfo = [];
    
    
    %% for each cue, store indexes in td, get reaction time
    for b = 1:numel(bumpList)
        for s = 1:numel(stimCodeList)
            % get idx in td
            if(stimCodeList(s) == -1) % handle no stim, just a bump case
                cueInfo(cueIdx).td_idx = find(isEqual(bumpList(b),[td.bumpMagnitude]) & ~[td.isStimTrial]);
            else % handle bump and stim (note, bump == 0 means no bump)
                cueInfo(cueIdx).td_idx = find(isEqual(bumpList(b),[td.bumpMagnitude]) & isEqual(stimCodeList(s),[td.stimCode]) & [td.isStimTrial]);
            end
            
            % find reaction time
            cueInfo(cueIdx).rt = mode([td.bin_size])*[(td(cueInfo(cueIdx).td_idx).idx_movement_on)] -  [td(cueInfo(cueIdx).td_idx).goCueTime];
            
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
    hold on
    
    fitData.x = [];
    fitData.y = [];
    for cueIdx = 1:numel(cueInfo)
        if(cueInfo(cueIdx).bumpMag ~= 0)
            plot(cueInfo(cueIdx).bumpMag,mean(cueInfo(cueIdx).rt),'k.','markersize',opts.MARKER_SIZE)
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
    %% plot as a function of stim code
%     
%     

        
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