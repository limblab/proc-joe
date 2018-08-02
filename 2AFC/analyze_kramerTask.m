%% determine filename and input data
    inputData.folderpath = 'C:\Users\jts3256\Desktop\Han_20180802_BD\';
    inputData.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';

    inputData.task='taskBD';
    inputData.ranBy='ranByJoseph'; 
    inputData.array1='arrayLeftS1'; 
    inputData.monkey='monkeyHan';
    inputData.labnum = 6;
    
    pwd=cd;
    cd(inputData.folderpath)
    fileList = dir('*.nev*');
    cd(pwd)
%% load in cds and extract data
    td_all = [];
    num_trials = 0;
        
    cds = commonDataStructure();
    cds.file2cds([inputData.folderpath fileList(1).name],inputData.task,inputData.ranBy,...
        inputData.monkey,inputData.labnum,inputData.array1,inputData.mapFileName,'recoverPreSync');
    cd(pwd);

% convert cds to trial data
    params.event_list = {'tgtDir';'isPrimaryTgt';...
        'bumpTime';'bumpDir';'bumpMagnitude';'isStimTrial';'stimCode'};
    params.trial_results = {'R','F'};
    params.extra_time = [1,2];
    params.exclude_units = [0:255];
    td_all = parseFileByTrial(cds,params);

%% plot psychometric curves
% note: all bumpDir's are relative to tgtDir. Is primary target determines
% if the tgt in tgtDir or the opposite one is the primary (1 = tgtDir, 0 =
% other target).

% logic for making these curves: For each bump direction, count the total
% rewards and trials to get a percent correct. Then, do 1-that percent for
% any bump dir > 90 as rewards here represent the opposite target

    bump_trial_mask = [td_all.isStimTrial] == 0 & [td_all.bumpMagnitude] > 0;
    td_bump = td_all(bump_trial_mask);
    bump_dirs = unique([td_bump.bumpDir]);
    bump_dirs = bump_dirs(~isnan(bump_dirs));
    bump_dirs = bump_dirs(bump_dirs <= 180);
    bump_correct = zeros(size(bump_dirs));
    bump_total = zeros(size(bump_dirs));
    
    for t = 1:numel(td_bump)
        if(td_bump(t).bumpDir > 180)
            bump_idx = find(bump_dirs == -1*(td_bump(t).bumpDir-360));
        else
            bump_idx = find(bump_dirs == td_bump(t).bumpDir);
        end
        
        if(td_bump(t).result == 'R')
            bump_correct(bump_idx) = bump_correct(bump_idx) + 1;
        end
        bump_total(bump_idx) = bump_total(bump_idx) + 1;
    end

    bump_percent_correct = bump_correct./bump_total;
    
    % 1-bump_percent_correct for dirs > 90 to make these figures
    bump_percent_correct(bump_dirs > 90) = 1 - bump_percent_correct(bump_dirs > 90);
    
    % make psychometric curve
    figure();
    plot(bump_dirs,bump_percent_correct,'k.','markersize',14)
    xlabel('Bump Direction')
    ylabel('Percent to 0-deg target')
    xlim([0,180])
    ylim([0,1])
    formatForLee(gcf)
    set(gca,'fontsize',14)
    
% fit psychometric curve
    s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', [.5 -0.5 .1 100],...
        'MaxFunEvals',10000,'MaxIter',1000,'TolFun',10E-8,'TolX',10E-8,...
        'Weights',bump_total/sum(bump_total));
    ft = fittype('a+b*(erf(c*(x-d)))','options',s);

    [g_mdl.fitObj,g_mdl.gof] = fit(bump_dirs', bump_percent_correct', ft);

    x_fit = [0:180];
    y_fit = feval(g_mdl.fitObj,x_fit);
%       
    hold on
    plot(x_fit,y_fit,'color','k','linewidth',1.5)
% plot for stim trials
    stim_codes = unique([td_all.stimCode]);
    stim_codes = stim_codes(~isnan(stim_codes));
    colors = {'r','b'};
    
    for st = 1:numel(stim_codes)
        stim_trial_mask = [td_all.isStimTrial] == 1 & [td_all.stimCode] == stim_codes(st) & [td_all.bumpMagnitude] > 0;
        td_stim = td_all(stim_trial_mask);
        bump_dirs = unique([td_stim.bumpDir]);
        bump_dirs = bump_dirs(~isnan(bump_dirs));
        bump_dirs = bump_dirs(bump_dirs <= 180);
        bump_correct = zeros(size(bump_dirs));
        bump_total = zeros(size(bump_dirs));

        for t = 1:numel(td_stim)
            if(td_stim(t).bumpDir > 180)
                bump_idx = find(bump_dirs == -1*(td_stim(t).bumpDir-360));
            else
                bump_idx = find(bump_dirs == td_stim(t).bumpDir);
            end
            if(td_stim(t).result == 'R')
                bump_correct(bump_idx) = bump_correct(bump_idx) + 1;
            end
            bump_total(bump_idx) = bump_total(bump_idx) + 1;
        end

        bump_percent_correct = bump_correct./bump_total;

        % 1-bump_percent_correct for dirs > 90 to make these figures
        bump_percent_correct(bump_dirs > 90) = 1 - bump_percent_correct(bump_dirs > 90);

        % make psychometric curve
%         figure();
        hold on
        plot(bump_dirs,bump_percent_correct,'.','color',colors{st},'markersize',14)
%         xlabel('Bump Direction')
%         ylabel('Percent to 0-deg target')
%         xlim([0,180])
%         ylim([0,1])
%         formatForLee(gcf)
%         set(gca,'fontsize',14)

        % fit psychometric curve
        s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', [.5 -0.5 .1 100],...
            'MaxFunEvals',10000,'MaxIter',1000,'TolFun',10E-8,'TolX',10E-8,...
            'Weights',bump_total/sum(bump_total));
        ft = fittype('a+b*(erf(c*(x-d)))','options',s);
        
        [g_mdl.fitObj,g_mdl.gof] = fit(bump_dirs', bump_percent_correct', ft);
        
        x_fit = [0:180];
        y_fit = feval(g_mdl.fitObj,x_fit);
%         
        plot(x_fit,y_fit,'color',colors{st},'linewidth',1.5)

    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
           