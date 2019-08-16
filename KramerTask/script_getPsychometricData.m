%% determine filename and input data
    input_data.folderpath = 'C:\Users\jts3256\Desktop\Duncan_BD_data\';
%     input_data.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
    input_data.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\right S1 20180919\SN 6251-001804.cmp';


    input_data.task='taskBD';
    input_data.ranBy='ranByJoseph'; 
    input_data.array1='arrayLeftS1'; 
    input_data.monkey='monkeyDuncan';
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
        params.exclude_units = [1:255];
        td_temp = parseFileByTrial(cds,params);
        if(isfield(td_temp,'LeftS1_unit_guide'))
            td_temp = rmfield(td_temp,'LeftS1_unit_guide');
            td_temp = rmfield(td_temp,'LeftS1_spikes');
        end
        if(isfield(td_temp,'force'))
            td_temp = rmfield(td_temp,'force');
        end
        td_all = [td_all, td_temp];
    end
    
    
%% get psychometric curve data
% note: all bumpDir's are relative to tgtDir. Is primary target determines
% if the tgt in tgtDir or the opposite one is the primary (1 = tgtDir, 0 =
% other target).
    
% logic for making these curves: For each bump direction, count the total
% rewards and trials to get a percent correct. Then, do 1-that percent for
% any bump dir > 90 as rewards here represent the opposite target

    input_data.max_trial_time = 10;
            
    input_data.num_bootstrap = 0;
    
    psych_data = getPsychometricCurveData(td_all,input_data);
    
%% plot psych data
    
    input_data.colors = {'k','r','b',[0,0.5,0],'m','g',[0.5,0.5,0.5]};
    input_data.psych_data_idx_list = [];
    input_data.plot_bootstrap = 0;
    
    input_data.date = '20190715';
    
    for i = 1:size(psych_data,2)
        input_data.axis = i;
        plotPsychometricCurve(psych_data,input_data);
    
    % legend
        ax = gca;
        ax.Children;
    %     uistack(ax.Children(3),'top');
    %     uistack(ax.Children(5),'top');
    %     l=legend('bump','0-deg','180-deg');
    %     set(l,'box','off');
    end

%% downsample TD
    td_all_use = binTD(td_all,50);

    
%% look at kinematics during bump for stim and non stim trials
    window = [0,400];
    input_data.colors = {'k','r','b',[0,0.5,0],'m','g',[0.5,0.5,0.5]};
    f=figure(); hold on;
    for trial = 1:numel(td_all)
        if(~isnan(td_all(trial).idx_bumpTime) && td_all(trial).bumpDir == 85 && td_all(trial).tgtDir == 0)
            if(~td_all(trial).isStimTrial)
                c = input_data.colors{1};
            else % it is a stim trial
                c = input_data.colors{td_all(trial).stimCode+2}; % 0 indexed, add another 1 to skip the bump color.
            end
            offset = td_all(trial).pos(td_all(trial).idx_bumpTime,:);
            
            plot(td_all(trial).pos(td_all(trial).idx_bumpTime+window(1):td_all(trial).idx_bumpTime+window(2),1)-offset(1),...
                td_all(trial).pos(td_all(trial).idx_bumpTime+window(1):td_all(trial).idx_bumpTime+window(2),2)-offset(2),...
                'color',c,'linewidth',1.5)

        end
        

    end
    % set x and y limits to be the same
    limits = [min(f.Children.XLim(1),f.Children.YLim(1)),max(f.Children.XLim(2), f.Children.YLim(2))];
    f.Children.XLim = limits;
    f.Children.YLim = limits;
    
    
%% Friedman test
    
    
    
    
    
%% choice direction stuff
%% start with simple neurometric curves -- plot firing rate as function of bump direction for each axis
    input_data.make_plot = 0;
    input_data.window = [0,0.4]; % s
    neurometric_data = makeNeurometricCurves(td_all_use,psych_data,input_data);
    
%% perform analysis done by (Ingaki, 2019) -- coding direction
% find n-dimensional vector (n = num neurons) where each entry is the
% average difference between spikes during correct left and right trials

    input_data.sample_rate = 0.5;
    input_data.bump_range = 90 + 50*[-1,1];
    for axis = 1:numel(neurometric_data)
        cd_data{axis} = getCodingDirection(neurometric_data{axis},input_data);
    end
    
%% predict decision and visualize projection onto cd over time
    input_data.bump_range = 90 + 90*[-1,1];

    for axis = 1:numel(neurometric_data)
        cd_pred_data{axis} = getCodingDirectionPredictions(neurometric_data{axis},cd_data{axis},input_data);
    end
    



    

    