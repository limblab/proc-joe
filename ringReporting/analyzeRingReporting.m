%% set file name and load file into cds

    folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\RingReporting\Duncan\Duncan_20181022\';
    mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Chips_12H1\map_files\left S1\SN 6251-001455.cmp';

    inputData.array = 'arrayNone';
    inputData.monkey = 'monkeyDuncan';
    inputData.ranBy = 'ranByJoe';
    inputData.lab = 6;
    inputData.mapFile = strcat('mapFile',mapFileName);
    inputData.task = 'taskRR';

    pwd=cd;
    cd(folderpath)
    fileList = dir('*nev*');

    cds = commonDataStructure();
    cds.file2cds(strcat(folderpath,fileList(1).name),inputData.array,inputData.monkey,inputData.ranBy,...
        inputData.lab,inputData.mapFile,inputData.task,'recoverPreSync','ignoreJumps','ignoreFilecat');
    cd(pwd);

%%
    opts.FIGURE_SAVE = 0;
    opts.FIGURE_PREFIX = 'Duncan_20181022_RR_training';
    opts.FIGURE_DIR = folderpath;

    opts.NUM_BINS_DIR = 8;
    opts.MAKE_FIGURES = 1;
    opts.PLOT_POLAR = 1;
    opts.MAX_TRIALS_PLOT = 2000;
    opts.CIRCLE_RADIUS = 5;
    opts.CIRCLE_DEPTH = 2;

    opts.DISTRIBUTION_BIN_SIZE = 5;
    opts.BUMP_MAGS = [1.25];
    opts.BUMP_NONLINEARITY = 0;

    behaviorData = processBehaviorRingReporting(cds,opts);

%% look at force onset
    params.event_list = {'bumpTime';'bumpDir';'otHoldTime'};
    params.trial_results = {'R','F'};
    td = parseFileByTrial(cds,params);
    td = removeBadTrials(td);

%%
    offset = [6,12];
    middle = [-3,33];
    bump_dir = [];
    acc_dir = [];
    reach_dir = [];
    for t = 1:numel(td)
        bump_dir(t) = td(t).target_direction*180/pi;
        pos_start = td(t).pos(td(t).idx_bumpTime+offset(1),:);
        pos_end = td(t).pos(td(t).idx_bumpTime+offset(2),:);
        acc_dir(t) = atan2(pos_end(2)-pos_start(2),pos_end(1)-pos_start(1))*180/pi;
        reach_dir(t) = atan2(td(t).pos(td(t).idx_otHoldTime,2)+middle(2),td(t).pos(td(t).idx_otHoldTime,1)+middle(1))*180/pi;
    end
    figure
    remove_idx = abs(bump_dir-acc_dir)*pi/180 > 1;
    bump_dir(remove_idx) = [];
    acc_dir(remove_idx) = [];
    
    plot(bump_dir*pi/180,acc_dir*pi/180,'.','markersize',12)
    hold on
%     plot(acc_dir,reach_dir,'.','markersize',12)
    plot([-pi,pi],[-pi,pi],'k--','linewidth',2)
    hold off
    formatForLee(gcf)
    ylabel('Movement direction')
    xlabel('Bump direction')
    xlim([-pi,pi])
    ylim([-pi,pi])
    set(gca,'fontsize',14)
%     f = gcf;
%     f.Name = 'Duncan_20180904_reachVsAcc_120ms';
%     saveFiguresLIB(f,folderpath,f.Name);


% fit the data with a function 
    [f,gof] = fit(acc_dir'*pi/180,bump_dir'*pi/180,'poly4');
    
    test = [-pi:0.01:pi];
   
    pred = feval(f,test);
    hold on
    plot(pred,test,'r','linewidth',2)
    