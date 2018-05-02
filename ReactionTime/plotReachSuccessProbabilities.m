function [g] = plotReachSuccessProbabilities(behaviorData,opts)

    % plot and fit pschyometric functions to the behavior data for bump and
    % stim cues
    
    %% configure opts
    opts = configureOpts(opts);
    outputData = [];
    
    %% plot bump behavior data
    bumpIndices = find(~isnan(behaviorData.bumpMag));
    
    figure();
    plot(behaviorData.bumpMag(bumpIndices),behaviorData.reachSuccess(bumpIndices),'.','color','k','markersize',opts.MARKER_SIZE);
    hold on
    ylim([0,1]);
    
    %% fit with a psychometric curve
    g = [];
    xData = linspace(min(behaviorData.bumpMag(bumpIndices)*0.9),max(behaviorData.bumpMag(bumpIndices)*1.1),100);

    if(opts.USE_ML_FIT)
        mlFit_fun = @(params)(sum((behaviorData.reachSuccess(bumpIndices)-(params(1)+params(2)*erf(params(3)*(behaviorData.bumpMag(bumpIndices)-params(4))))).^2));
        min_options = optimset('MaxIter',10000);
        g = fminsearch(mlFit_fun, [.45 .4 .05 1],min_options);

        yData = g(1) + g(2)*erf(g(3)*(xData-g(4))); 
    else
        s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', [.5 .5 .1 90], 'Lower', [0 0 0 0], 'Upper', [1 1 10 1.5],...
            'MaxFunEvals',10000,'MaxIter',1000,'TolFun',10E-8,'TolX',10E-8);
        ft = fittype('a+b*(erf(c*(x-d)))','options',s);
        [g.fitObj,g.gof] = fit(behaviorData.bumpMag(bumpIndices)', behaviorData.reachSuccess(bumpIndices)', ft);

        yData = feval(g.fitObj,xData);
    end
    
    plot(xData,yData,'k--','linewidth',opts.LINE_WIDTH);
end

function [opts] = configureOpts(optsInput)

    opts = [];
    
    opts.LINE_WIDTH = 1.5;
    
    opts.COLORS = {'r',[0 0.5 0],'b','k','m',[0.5,0.5,0.2]};
    opts.MARKER_SIZE = 20;

    opts.USE_ML_FIT = 0;
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