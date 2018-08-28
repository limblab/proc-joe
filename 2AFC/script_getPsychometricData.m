%% determine filename and input data
    input_data.folderpath = 'C:\Users\jts3256\Desktop\Han_20180803_BD\';
    input_data.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';

    input_data.task='taskBD';
    input_data.ranBy='ranByJoseph'; 
    input_data.array1='arrayLeftS1'; 
    input_data.monkey='monkeyHan';
    input_data.labnum = 6;
    
    pwd=cd;
    cd(input_data.folderpath)
    fileList = dir('*.nev*');
    cd(pwd)
%% load in cds and extract data
    td_all = [];

    for fileNumber = 1:numel(fileList)  
        cds = commonDataStructure();
        cds.file2cds([input_data.folderpath fileList(fileNumber).name],input_data.task,input_data.ranBy,...
            input_data.monkey,input_data.labnum,input_data.array1,input_data.mapFileName,'recoverPreSync');
        cd(pwd);

    % convert cds to trial data
        params.event_list = {'tgtDir';'isPrimaryTgt';...
            'bumpTime';'bumpDir';'bumpMagnitude';'isStimTrial';'stimCode'};
        params.trial_results = {'R','F'};
        params.extra_time = [1,2];
        params.exclude_units = [0:255];
        td_temp = parseFileByTrial(cds,params);

        td_all = [td_all, td_temp];
    end
%% get psychometric curve data
% note: all bumpDir's are relative to tgtDir. Is primary target determines
% if the tgt in tgtDir or the opposite one is the primary (1 = tgtDir, 0 =
% other target).
    
% logic for making these curves: For each bump direction, count the total
% rewards and trials to get a percent correct. Then, do 1-that percent for
% any bump dir > 90 as rewards here represent the opposite target
    input_data.max_trial_time = 50000;
    input_data.min_trial_time = 0;
    
    input_data.num_bootstrap = 0;
        
    psych_data = getPsychometricCurveData(td_all,input_data);
    
%% plot psych data
    
    input_data.colors = {'k','r','b',[0,0.5,0]};
    f = plotPsychometricCurve(psych_data(1:end),input_data);

% legend
    ax = gca;
    ax.Children;
    uistack(ax.Children(3),'top');
    uistack(ax.Children(5),'top');
    l=legend('bump','0-deg','180-deg');
    set(l,'box','off');
    
%% save
    f = gcf;
    f.Name = 'Han_20180805_BD_psychometricCurve';
    saveFiguresLIB(f,input_data.folderpath,f.Name);
    
    
%% draw circle with line in 0 deg direction
    f = figure;
    theta = 90+[0,0;180,180];
    
    rho = [0,1];
    polarplot(theta(1,:)*pi/180,rho,'color','r','linewidth',3);
    hold on
    polarplot(theta(2,:)*pi/180,rho,'color','b','linewidth',3);
    ax = f.Children;
    ax.RTickLabel = {};
    ax.ThetaTick = [theta(1,1),theta(2,1)];
    ax.ThetaTickLabel = {'0','180'};
    ax.ThetaMinorGrid = 'on';
    ax.LineWidth = 2;

    f.Name = 'Han_20180805_BD_polarDir_90deg270deg';
    saveFiguresLIB(f,input_data.folderpath,f.Name);
    
           