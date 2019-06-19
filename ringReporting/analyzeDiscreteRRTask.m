%% set file name and load file into cds

    input_data.folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\KramerTask\Han_20190521_BDmanyTgt_training\';
%     mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp';
    mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
    
    
    input_data.date = '20190520';
    input_data.array = 'arrayLeftS1';
    input_data.monkey = 'monkeyHan';
    input_data.ranBy = 'ranByJoe';
    input_data.lab = 6;
    input_data.mapFile = strcat('mapFile',mapFileName);
    input_data.task = 'taskBD';

    pwd=cd;
    cd(input_data.folderpath)
    fileList = dir('*nev*');
    input_data.fileNum = 1;
  
    input_data.center_y = -33;
    input_data.center_x = 3;
    input_data.num_bins = 8;
%% convert file to cds and td
    td_all = [];
    for f = 1:numel(fileList)
        cds = commonDataStructure();
        cds.file2cds(strcat(input_data.folderpath,fileList(f).name),input_data.array,input_data.monkey,input_data.ranBy,...
            input_data.lab,input_data.mapFile,input_data.task,'recoverPreSync','ignoreJumps','ignoreFilecat');
        cd(pwd);
    
        params.event_list = {'goCueTime';'tgtDir';'bumpDir';'bumpTime';...
            'bumpMagnitude';'stimCode';'showOuterTarget';'numTargets';'angleTolerance';...
            'correctAngle'};
        params.trial_results = {'R','F'};
        params.extra_time = [1,2];
        params.include_ts = 0;
        params.exclude_units = [1:255];
        td_temp = parseFileByTrial(cds,params);
        
        td_all = [td_all,td_temp];
    end
    
    input_data.tgt_width = mode([cds.trials.angleTolerance]);
    


%% get bump trials, plot tgt_dir vs reach_dir. 
% plot percentage correct in bins around a polar plot

    td_bump = td_all(~isnan([td_all.bumpDir]));
    
    % get pos at idx_endTime
    reach_angles = zeros(numel(td_bump),1);
    correct_angle = zeros(numel(td_bump),1);
    is_rewarded = zeros(numel(td_bump),1);
    for t = 1:numel(td_bump)
        correct_angle(t) = td_bump(t).correctAngle;
        
        reach_angles(t) = getReachAngle(td_bump,t,[input_data.center_x,input_data.center_y]);
%         is_rewarded(t) = td_bump(t).result == 'R';
        is_rewarded(t) = abs(reach_angles(t) - correct_angle(t)) < td_bump(t).angleTolerance;
    end
    f=figure(); hold on
    f.Name = [input_data.monkey(7:end),'_',input_data.date,'_bump_reachVsTgt'];
    
    % unity black line
    plot([-180,180],[-180,180],'k--','linewidth',2);
    % error tolerance lines (4 to compensate for wrap arounds)
    plot([-180,180],[-180,180]+input_data.tgt_width/2,'r--','linewidth',2);
        plot([-180,180],[-180,180]+360-input_data.tgt_width/2,'r--','linewidth',2);
        
    plot([-180,180],[-180,180]-input_data.tgt_width/2*180/180,'r--','linewidth',2);
        plot([-180,180],[-180,180]-360+input_data.tgt_width/2,'r--','linewidth',2);
      
    % plot reach vs tgt
    plot(correct_angle*180/pi,reach_angles*180/pi,'.','markersize',16)
    % format
    xlim([-180,180]); ylim([-180,180]);
    formatForLee(gcf); set(gca,'fontsize',14)
    xlabel('Tgt dir'); ylabel('Reach dir');
    
    % get percent correct in bins around a circle
    bin_edges = [-pi:2*pi/input_data.num_bins:pi];
    percent_correct = zeros(numel(bin_edges)-1,1);
    for i = 1:numel(bin_edges)-1
        trial_mask = correct_angle >= bin_edges(i) & correct_angle < bin_edges(i+1);
        percent_correct(i) = sum(is_rewarded(trial_mask))/sum(trial_mask);
    end
    
    f=figure();
    f.Name = [input_data.monkey(7:end),'_',input_data.date,'_bump_percentCorrect'];
    polarplot([bin_edges(1:end-1),bin_edges(1)]+mode(diff(bin_edges))/2,[percent_correct;percent_correct(1)],...
        '--.','markersize',20,'linewidth',2)
    f.Children(1).RLim = [0,1]; % radius limits from 0 to 1
    ax = gca;
	ax.ThetaTickLabel = {'0','30','60','90','120','150','180','-150','-120','-90','-60','-30'};
    
 