%% set file name and load file into cds

    folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\RingReporting\Duncan\Duncan_20190411\';
    mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp';

    inputData.array = 'arrayLeftS1';
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
    opts.FIGURE_PREFIX = 'Duncan_20190410_RR_training';
    opts.FIGURE_DIR = folderpath;

    opts.NUM_BINS_DIR = 8;
    opts.MAKE_FIGURES = 1;
    opts.PLOT_POLAR = 1;
    opts.MAX_TRIALS_PLOT = 2000;
    opts.CIRCLE_RADIUS = 5.75;
    opts.CIRCLE_DEPTH = 2;

    opts.DISTRIBUTION_BIN_SIZE = 5;
    opts.BUMP_MAGS = [];
    opts.BUMP_NONLINEARITY = 0;

    behaviorData = processBehaviorRingReporting(cds,opts);

%% look at force onset
    params.event_list = {'bumpTime';'bumpDir';'otHoldTime'};
    params.trial_results = {'R','F'};
    td = parseFileByTrial(cds,params);
    td = removeBadTrials(td);

%%
    offset = [0,30];
    middle = [-3,33];
    tgt_dir = []; bump_dir = [];
    acc_dir = [];
    reach_dir = [];
    for t = 1:numel(td)
        bump_dir(t) = td(t).bumpDir;
        tgt_dir(t) = td(t).target_direction*180/pi;
        pos_start = td(t).pos(td(t).idx_bumpTime+offset(1),:);
        pos_end = td(t).pos(td(t).idx_bumpTime+offset(2),:);
        acc_dir(t) = atan2(pos_end(2)-pos_start(2),pos_end(1)-pos_start(1))*180/pi;
        reach_dir(t) = atan2(td(t).pos(td(t).idx_otHoldTime,2)+middle(2),td(t).pos(td(t).idx_otHoldTime,1)+middle(1))*180/pi;
        dist_dir(t) = sqrt((td(t).pos(td(t).idx_otHoldTime,2)+middle(2)).^2 + (td(t).pos(td(t).idx_otHoldTime,1)+middle(1)).^2);
    end
    remove_idx = abs(tgt_dir-acc_dir)*pi/180 > 1;
    tgt_dir(remove_idx) = []; bump_dir(remove_idx) = [];
    acc_dir(remove_idx) = [];
    dist_dir(remove_idx) = [];
    reach_dir(remove_idx) = [];
%
    figure

    idx_change = bump_dir*pi/180 < -2.8;
    acc_dir(idx_change) = acc_dir(idx_change)-360;
    
    plot(bump_dir*pi/180,acc_dir*pi/180,'.','markersize',12)
    hold on
%     plot(acc_dir,reach_dir,'.','markersize',12)
    plot([-pi,pi],[-pi,pi],'k--','linewidth',2)
    hold off
    formatForLee(gcf)
    ylabel('acc dir')
    xlabel('bump dir')
    xlim([-pi,pi])
    ylim([-pi,pi])
    set(gca,'fontsize',14)
    
% % fit the data with a function 
    [f,gof] = fit(acc_dir'*pi/180,bump_dir'*pi/180,'poly5');
    
    acc_test = [-pi:0.01:pi];
   
    bump_pred = feval(f,acc_test);
    hold on
    plot(bump_pred,acc_test,'r','linewidth',2)
%     
%     
%     figure
% 
%     plot(acc_dir*pi/180,reach_dir*pi/180,'.','markersize',12)
%     hold on
% %     plot(acc_dir,reach_dir,'.','markersize',12)
%     plot([-pi,pi],[-pi,pi],'k--','linewidth',2)
%     hold off
%     formatForLee(gcf)
%     ylabel('reach direction')
%     xlabel('acc dir')
%     xlim([-pi,pi])
%     ylim([-pi,pi])
%     set(gca,'fontsize',14)
    
%%

    norm_dist_dir = dist_dir;

    figure();
    
    plot(tgt_dir*pi/180,norm_dist_dir,'.','markersize',12)
    hold on
    formatForLee(gcf)
    ylabel('Distance')
    xlabel('Bump direction')
    xlim([-pi,pi])
    set(gca,'fontsize',14)
%     f = gcf;
%     f.Name = 'Duncan_20180904_reachVsAcc_120ms';
%     saveFiguresLIB(f,folderpath,f.Name);


% fit the data with a function 
    [f_dist,gof] = fit(tgt_dir'*pi/180,norm_dist_dir','a*cos(b*x+c)+d');
    
    test = [-pi:0.01:pi];
   
    pred = feval(f_dist,test);
    hold on
    plot(test,pred,'r','linewidth',2)