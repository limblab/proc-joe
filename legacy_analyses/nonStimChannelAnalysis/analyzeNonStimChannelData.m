%% set folder names 
    inputData.folderpath = 'C:\Users\joh8881\Desktop\Han_20190923_trains_noAmp\';
    inputData.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';

    % inputData.mapFileName = 'mapFileR:\limblab-archive\Retired Animal Logs\Monkeys\Chips_12H1\map_files\left S1\SN 6251-001455.cmp';
%     inputData.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp';


    folderpath = inputData.folderpath; % rest of code uses folderpath currently...may have switched this, not 100% certain

    inputData.task='taskCObump';
    inputData.ranBy='ranByJoseph'; 
    inputData.array1='arrayLeftS1'; 
    inputData.monkey='monkeyHan';
    inputData.labnum = 6;

    pwd=cd;
    cd(inputData.folderpath)
    fileList = dirSorted('*spikesExtracted.nev*');
    stimInfoFileList = dirSorted('*stimInfo*');


%% extract relevant data for all units -- recommend saving arrayData after this step
    tic

    optsExtract.STIMULI_RESPONSE = 'all';
    optsExtract.STIMULATIONS_PER_TRAIN = 11;
    optsExtract.STIMULATION_BATCH_SIZE = 1000;

    optsExtract.USE_STIM_CODE = 0;
    optsExtract.STIM_ELECTRODE = {};
    optsExtract.CHAN_LIST = [];

    optsExtract.PRE_TIME = 150/1000; % made negative in the function
    optsExtract.POST_TIME = 500/1000;

    optsExtract.BIN_SIZE = 5/1000;
    optsExtract.TIME_AFTER_STIMULATION_WAVEFORMS = 10/1000;
    optsExtract.USE_ANALOG_FOR_STIM_TIMES = 1; % this uses the analog sync line to get stim times, not sure why you would want to do anything else
    optsExtract.GET_KIN = 1;
    optsExtract.GET_FORCE = 0;
    
    arrayData = extractDataAroundStimulations(inputData,fileList,stimInfoFileList,optsExtract);

    toc
    
    
%% plot raster, and PSTH for the given unit above
for arrIdx = 4%1:numel(arrayData)
% arrIdx = 1;
    % plot raster, and PSTH for the given unit above

%     optsPlotFunc.BIN_SIZE = optsExtract.BIN_SIZE;
    optsPlotFunc.MAKE_SUBPLOTS = 1;
    optsPlotFunc.PLOT_TITLE = 0;    

    optsPlotFunc.BIN_SIZE = mode(diff(arrayData{1}.binEdges{1,1}));
    optsPlotFunc.FIGURE_SAVE = 0;
    optsPlotFunc.FIGURE_DIR = inputData.folderpath;
    optsPlotFunc.FIGURE_PREFIX = 'Han_20190923';

    optsPlotFunc.PRE_TIME = 5/1000;
    optsPlotFunc.POST_TIME = 10/1000;
    optsPlotFunc.SORT_DATA = '';

    optsPlotFunc.PLOT_AFTER_STIMULATION_END = 1;
    optsPlotFunc.STIMULATION_LENGTH = 0.4+0.53;
    
    rasterPlots = plotRasterStim(arrayData{arrIdx},arrayData{arrIdx}.NN,optsPlotFunc);

    optsPlotFunc.PLOT_ALL_ONE_FIGURE = 0;
    optsPlotFunc.PLOT_LINE = 1;
    optsPlotFunc.PLOT_ALL_WAVES_ONE_FIGURE = 0;
% %     
%     PSTHPlots = plotPSTHStim(arrayData{arrIdx},1,optsPlotFunc);

end


%% get response amplitude for each condition/neuron
    
    baseline_list = [];
    for unit = 1:numel(arrayData)
        arrayData{unit} = getResponseAmplitude(arrayData{unit},unit,0.2,2); % bin size (ms), window width (bins on either side of peak)
        baseline_list(unit) = mean(mean(arrayData{unit}.mean_baseline));
    end

%% heatmap across whole array

%     inputData.folderpath = 'C:\Users\joh8881\Desktop\Han_20190930_trains_noAmp\';
    inputData.folderpath = 'E:\Data\Joseph\Duncan_stim_data\Duncan_20191026_trains_noAmp\';

    inputData.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
%     inputData.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp';
%     inputData.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\right S1 20180919\SN 6251-001804.cmp';
%     inputData.folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\StimRecData\Han\Multielec\Han_20190406_2elecs\';
%    inputData.mapFileName = 'mapFileZ:\Basic_Sciences\Phys\L_MillerLab\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
%    inputData.mapFileName = 'mapFileZ:\Basic_Sciences\Phys\L_MillerLab\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp';
    

    opts.STIM_ELECTRODE_PLOT = [1:size(arrayData{1}.binEdges,1)];
    opts.WAVEFORM_TYPES_PLOT = 4%[1:size(arrayData{1}.binEdges,2)];
    opts.MAX_RATIO = 1;
    opts.MIN_RATIO = -1;
    opts.LOG_SCALE = 0;
    opts.LOG_PARAM = 9;
    
    [stimHeatmapHandles] = plotHeatmaps(arrayData,inputData.mapFileName(8:end),opts);
    

%% Plot evoked response vs distance (400um between channels)
    
    inputData.use_responsive_only = 0;
    inputData.color_idx = [4,3,2,5];
    inputData.use_boxplots = 0;
    
    inputData.distance_bin_size = 300;
    inputData.dist_offset = [-75,-25,25,75];
    inputData.linewidth = 1.75;
    inputData.box_width = 30;
    inputData.outlier_marker = '.';
    inputData.outlier_marker_size = 12;
    
    
    amp_data = plotResponseAmpVsDistance(arrayData,inputData);

    
%% Plot latency vs distance (400um between channels)
    inputData.color_idx = [4,3,2,5];
    inputData.use_boxplots = 0;
    
    inputData.distance_bin_size = 600;
    inputData.dist_offset = [-75,-25,25,75];
    inputData.linewidth = 1.75;
    inputData.box_width = 30;
    inputData.outlier_marker = '.';
    inputData.outlier_marker_size = 12;
    
    latency_data = plotLatencyVsDistance(arrayData,inputData);

 
%% plot PD difference heatmap and stim heatmap for each condition
% PD differences compared to center channel
%     [pdHeatmapData,pdAlphaData] = getHeatmapDataPD(td_all,pd_all,mapData);
%    inputData.mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
   inputData.mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp';

    for i = 1:numel(heatmapData)
        [combinedHeatmapHandles] = plotPDandStimHeatmaps(heatmapData{i},pdHeatmapData,inputData.mapFileName);
    end
    
%% compare stim response and neural response from the task for each condition
    testAngles = -180:1:180;
    angleDifferences = zeros(numel(heatmapData),1);
    bestAngles = zeros(numel(heatmapData),1);
    compareData = {};
    for i = 1:numel(heatmapData)
        compareData{i} = compareHeatmaps(heatmapData{i},pdHeatmapData,testAngles,inputData.mapFileName);
        figure();
        plot(testAngles,compareData{i}.metricAllAngles);
        hold on
        plot(compareData{i}.centerChanPD + [0,0], [-0.2,0.2],'r--')
        angleDifferences(i) = angleDiff(compareData{i}.centerChanPD,compareData{i}.bestAngle,0,0);
        bestAngles(i) = compareData{i}.bestAngle;
        xlabel('Angle');
        ylabel('Metric');
    end
  
    figure();
    subplot(2,1,1)
    histogram(angleDifferences,18)
    subplot(2,1,2)
    histogram(bestAngles,[-180:20:180])
    xlim([-180,180])
    
%% bootstrap stim channels to get a uniform distribution of PDs, then compute metrics
    testAngles = -180:1:180;
    numTests = 200;
    % pdHeatmapData, alphaData
    pdHeatmapWeight = [];
    row = []; col = [];
    pdAll = [];
    angleDifferences = zeros(numTests,numel(heatmapData));
    bestAngles = zeros(numTests,numel(heatmapData));
    pdHeatmapFlatData = [];
    
    for rowIdx = 1:size(pdHeatmapData,1)
        for colIdx = 1:size(pdHeatmapData,2)
            if(alphaData(rowIdx,colIdx) == 1)
                weight_idx = find(pdHeatmapData(rowIdx,colIdx) > PDBinEdges,1,'last');
                pdHeatmapWeight(end+1,1) = 1/(PDBinCount(weight_idx));
                pdHeatmapFlatData(end+1,1) = pdHeatmapData(rowIdx,colIdx);
                row(end+1,1) = rowIdx; col(end+1,1) = colIdx;
            end
        end
    end
    
    for t = 1:numTests
        [~,pdIdx] = datasample(pdHeatmapWeight,numel(pdHeatmapWeight),'Replace',true,'Weights',pdHeatmapWeight);
        % make pdHeatmapData based on datasample
        pdHeatmapDataSample = nan(size(pdHeatmapData));
        for i = 1:numel(pdIdx)
            if(isnan(pdHeatmapDataSample(row(pdIdx(i)),col(pdIdx(i)))))
                pdHeatmapDataSample(row(pdIdx(i)),col(pdIdx(i))) = pdHeatmapFlatData(pdIdx(i));
            else
                pdHeatmapDataSample(row(pdIdx(i)),col(pdIdx(i))) = pdHeatmapDataSample(row(pdIdx(i)),col(pdIdx(i))) + pdHeatmapFlatData(pdIdx(i));
            end
        end
        
        for i = 1:numel(heatmapData)            
            compareData{i} = compareHeatmaps(heatmapData{i},pdHeatmapDataSample,testAngles,inputData.mapFileName);
            angleDifferences(t,i) = angleDiff(compareData{i}.centerChanPD,compareData{i}.bestAngle,0,0);
            bestAngles(t,i) = compareData{i}.bestAngle;
        end
        pdAll = [pdAll; pdHeatmapFlatData(pdIdx)];
    end
    
    figure();
    subplot(2,1,1)
    histogram(angleDifferences,18)
    subplot(2,1,2)
    histogram(pdAll,18);
%% shuffle stim data and see what effect that has

    testAngles = -180:2:180;
    numShuffles = 500;
    numPlot = 25;
    shuffledCompareData = cell(numShuffles,1);
    
    f=figure(); f.Position = [269 529 1019 420];
    subplot(1,2,1)
    bestAngles = zeros(numShuffles,1);
    for i = 1:numShuffles
        tempHeatmapData = heatmapData{1};%heatmapData{ceil(rand()*numel(heatmapData))};
        tempHeatmapData.dataRatioScaled = shuffleMatrix(tempHeatmapData.dataRatioScaled);
        
        shuffledCompareData{i} = compareHeatmaps(tempHeatmapData,pdHeatmapData,testAngles,inputData.mapFileName);
        if(i < numPlot)
            plot(testAngles,shuffledCompareData{i}.metricAllAngles,'r','linewidth',0.5);
        end
        hold on
        
        bestAngles(i) = shuffledCompareData{i}.bestAngle;
    end
    subplot(1,2,2)
    histogram(bestAngles,18);
 