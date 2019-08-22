%% set file name and load file into cds

    input_data.folderpath = 'C:\Users\jts3256\Desktop\Duncan_stim_evoked_move\Duncan_20190814_RR_50uA_headstage1\';
    input_data.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp';
%     input_data.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';

    input_data.date = '20190814';
    input_data.array = 'arrayLeftS1';
    input_data.monkey = 'monkeyDuncan';
    input_data.ranBy = 'ranByJoe';
    input_data.lab = 6;
    input_data.task = 'taskRR';

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
            input_data.lab,input_data.mapFileName,input_data.task,'recoverPreSync','ignoreJumps','ignoreFilecat');
        cd(pwd);
    
        if(strcmpi(input_data.task,'RR'))
            params.event_list = {'otHoldTime';'goCueTime';'tgtDir';'bumpDir';'bumpTime';'bumpMagnitude';'stimTime';'stimCode';'catchTrial';'showOuterTarget'};
        else
            params.event_list = {'goCueTime';'bumpDir';'bumpTime';'bumpMagnitude';'stimTime';'stimCode'};
        end
        params.trial_results = {'R','F','A','I'};
        params.extra_time = [1,2];
        params.include_ts = 0;
        params.exclude_units = [255];
        td_temp = parseFileByTrial(cds,params);
%         td_temp = td_temp(~isnan([td_temp.catchTrial]));
        
        td_all = [td_all,td_temp];
    end
    
%     input_data.tgt_width = mode([cds.trials.tgtWidth]);
    
%% use sync to get stim times:
    aIdx = 3; syncName = 'ainp16';
    stimOnIdx=find(diff(cds.analog{aIdx}.(syncName)-mean(cds.analog{aIdx}.(syncName))>3)>.5);
    for j=1:numel(stimOnIdx)
        if j<numel(stimOnIdx)
            next=stimOnIdx(j+1);
        else
            next=numel(cds.analog{aIdx}.(syncName));
        end
    end

    stimOnTimes = cds.analog{aIdx}.t(stimOnIdx);
    timeDiffs = [];
    td_all_adj = td_all;
    for tr = 1:size(cds.trials,1)
        stimOnIdx = find(stimOnTimes > cds.trials.stimTime(tr),1,'first');
        timeDiff = stimOnTimes(stimOnIdx) - cds.trials.stimTime(tr);
        
        tdAllIdx = find([td_all_adj.trial_id] == cds.trials.number(tr));
        if(~isempty(tdAllIdx) && ~isempty(stimOnIdx))
            td_all_adj(tdAllIdx).idx_stimTime = td_all_adj(tdAllIdx).idx_stimTime + floor(timeDiff*1000);
            timeDiffs(end+1,1) = floor(timeDiff*1000);
        end
    end
    
%% load in pattern file to set predicted direction and have relevant info

%% get bump trials, plot tgt_dir vs reach_dir. 
% plot percentage correct in bins around a polar plot

    td_bump = td_all(~isnan([td_all.bumpDir]) & ~[td_all.catchTrial] & isnan([td_all.stimCode]));
    
    % get pos at idx_endTime
    reach_angles_bump = zeros(numel(td_bump),1);
    tgt_angles_bump = zeros(numel(td_bump),1);
    is_rewarded = zeros(numel(td_bump),1);
    for t = 1:numel(td_bump)
        tgt_angles_bump(t) = td_bump(t).target_direction;
        
        reach_angles_bump(t) = getReachAngle(td_bump,t,[input_data.center_x,input_data.center_y]);
        is_rewarded(t) = td_bump(t).result == 'R';
%         is_rewarded(t) = abs(reach_angles(t) - tgt_angles(t))*180/pi < input_data.tgt_width/2;
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
    plot(tgt_angles_bump*180/pi,reach_angles_bump*180/pi,'.','markersize',16)
    % format
    xlim([-180,180]); ylim([-180,180]);
    formatForLee(gcf); set(gca,'fontsize',14)
    xlabel('Tgt dir'); ylabel('Reach dir');
    
    % get percent correct in bins around a circle
    bin_edges = [-pi:2*pi/input_data.num_bins:pi];
    percent_correct = zeros(numel(bin_edges)-1,1);
    for i = 1:numel(bin_edges)-1
        trial_mask = tgt_angles_bump >= bin_edges(i) & tgt_angles_bump < bin_edges(i+1);
        percent_correct(i) = sum(is_rewarded(trial_mask))/sum(trial_mask);
    end
    
    f=figure();
    f.Name = [input_data.monkey(7:end),'_',input_data.date,'_bump_percentCorrect'];
    polarplot([bin_edges(1:end-1),bin_edges(1)]+mode(diff(bin_edges))/2,[percent_correct;percent_correct(1)],...
        '--.','markersize',20,'linewidth',2)
    f.Children(1).RLim = [0,1]; % radius limits from 0 to 1
    ax = gca;
	ax.ThetaTickLabel = {'0','30','60','90','120','150','180','-150','-120','-90','-60','-30'};
    
%% get catch trials, plot end point reaches around a circle and x's in middle if he stayed in middle
    td_catch = td_all([td_all.catchTrial]==1 & ~isnan([td_all.idx_otHoldTime]));
    reach_angles = zeros(numel(td_catch),1);
    
    for t = 1:numel(td_catch)
        reach_angles(t) = getReachAngle(td_catch,t,[input_data.center_x,input_data.center_y])*180/pi; % output is in radians.
    end
    
    % histogram instead of polar plot
    f=figure();
    f.Name = [input_data.monkey(7:end),'_',input_data.date,'_catchTrials'];
    bin_edges = [-180:360/input_data.num_bins:180];
    num_reaches = histcounts(reach_angles,bin_edges)/numel(reach_angles);
    histogram('BinEdges',bin_edges,'BinCounts',num_reaches)
    xlabel('Angle')
    ylabel('Proportion of reaches');
    formatForLee(gcf)
    set(gca,'fontsize',14)
    
%% for stim trials, plot pred_dir vs reach_dir total, also reach dir around circle for each code
    td_stim = td_all(~isnan([td_all.stimCode]) & ~[td_all.catchTrial]);
    
    % get pos at idx_endTime
    stim_code = zeros(numel(td_stim),1);
    reach_angles = zeros(numel(td_stim),1);
    pred_angles = zeros(numel(td_stim),1);
    is_rewarded = zeros(numel(td_stim),1);
    is_biomimetic = zeros(numel(td_stim),1);
    for t = 1:numel(td_stim)
        stim_code(t) = td_stim(t).stimCode+1; % td.stim code is 0-15, we want 1-16 so add 1
        if(exist('pattern_data') && stim_code(t) <= numel(pattern_data.pred_dirs))
            pred_angles(t) = pattern_data.pred_dirs(stim_code(t)); % in degrees
            is_biomimetic(t) = pattern_data.is_biomimetic(stim_code(t));
        elseif(td_stim(t).bumpMagnitude > 0 && stim_code(t) == 2)
            pred_angles(t) = td_stim(t).tgtDir;
            is_biomimetic(t) = 0;
        else
            pred_angles(t) = nan;
            is_biomimetic(t) = -1;
        end
        reach_angles(t) = getReachAngle(td_stim,t,[input_data.center_x,input_data.center_y])*180/pi; % output is in radians.
        is_rewarded(t) = td_stim(t).result == 'R';
        
    end
    f=figure(); hold on
    f.Name = [input_data.monkey(7:end),'_',input_data.date,'_stim_reachVsPred'];
    % unity black line
    plot([-180,180],[-180,180],'k--','linewidth',2);
    % error tolerance lines (4 to compensate for wrap arounds)
    plot([-180,180],[-180,180]+input_data.tgt_width/2,'r--','linewidth',2);
        plot([-180,180],[-180,180]+360-input_data.tgt_width/2,'r--','linewidth',2);
        
    plot([-180,180],[-180,180]-input_data.tgt_width/2*180/180,'r--','linewidth',2);
        plot([-180,180],[-180,180]-360+input_data.tgt_width/2,'r--','linewidth',2);
      
        
    % plot bump data 
        plot(tgt_angles_bump*180/pi,reach_angles_bump*180/pi,'.','markersize',12)
    % plot reach vs tgt for bio patterns as dots
    plot(pred_angles(is_biomimetic==1),reach_angles(is_biomimetic==1),'.','color',[0,0,0.8],'markersize',20,'linewidth',2)
    % plot reach vs tgt for nonbio patterns as x's
    plot(pred_angles(~is_biomimetic),reach_angles(~is_biomimetic),'x','color',[0,0.5,0],'markersize',20,'linewidth',2)
    


    % format
    xlim([-180,180]); ylim([-180,180]);
    formatForLee(gcf); set(gca,'fontsize',14)
    xlabel('Tgt dir'); ylabel('Reach dir');
    
    % for each stim code (pattern,elec,etc.), plot a histogram with end
    % point reaches
    for i = 1:numel(unique(stim_code))
        f=figure();
        f.Name = [input_data.monkey(7:end),'_',input_data.date,...
            '_stimReaches_patternNum',num2str(i),'_isbio',num2str(is_biomimetic(find(stim_code==i,1,'first')))];
        
        bin_edges = [-180:360/input_data.num_bins:180];
        num_reaches = histcounts(reach_angles(stim_code==i),bin_edges);
        histogram('BinEdges',bin_edges,'BinCounts',num_reaches)
        
        if(exist('pattern_data'))
            hold on
            plot(pattern_data.pred_dirs(i)+[0,0],[0,max(num_reaches)],'r-','linewidth',2)
        end

        xlabel('Angle')
        ylabel('Number of reaches');
        formatForLee(gcf)
        set(gca,'fontsize',14)
    end
    
%% time to target for bumps, stim and catch trials
    spread = 0.25;

    time_to_target_stim = ([td_stim.idx_otHoldTime] - [td_stim.idx_goCueTime])*td_stim(1).bin_size;
    time_to_target_bump = ([td_bump.idx_otHoldTime] - [td_bump.idx_goCueTime])*td_bump(1).bin_size;
    time_to_target_catch = ([td_catch.idx_otHoldTime] - [td_catch.idx_goCueTime])*td_catch(1).bin_size;

    
    f = figure();
    f.Name = [input_data.monkey(7:end),'_',input_data.date,'_timeToTarget'];
    
    % bump
    counter = 1;
    plot(rand(size(time_to_target_bump))*spread-spread/2 + counter,time_to_target_bump,'.','markersize',12,'color',getColorFromList(1,0));
    hold on
    counter = counter + 1;
    % catch
    plot(rand(size(time_to_target_catch))*spread-spread/2 + counter,time_to_target_catch,'.','markersize',12,'color',getColorFromList(1,1));
    counter = counter + 1;
    % stim
    for stim_idx = 1:numel(unique(stim_code))
        stim_mask = [td_stim.stimCode]+1 == stim_idx;
        plot(rand(size(time_to_target_stim(stim_mask)))*spread-spread/2 + counter,time_to_target_stim(stim_mask),'.','markersize',12,'color',getColorFromList(1,2));
        counter = counter+1;
    end
    
    l = legend('bump','catch','stim'); set(l,'box','off');
    
    xlabel('Condition')
    ylabel('Time to target (s)');
    formatForLee(gcf)
    set(gca,'fontsize',14)
    
    
    
    
    
%     
%     
    %% plot reach kinematics,
    num_plot = 10;
    counter = 0;
    figure();
    for t = 1:numel(td_stim)
        if(counter < num_plot && td_stim(t).stimCode == 0)
            counter = counter + 1;
            plot(td_stim(t).pos(td_stim(t).idx_stimTime-50:td_stim(t).idx_otHoldTime,1),...
                td_stim(t).pos(td_stim(t).idx_stimTime-50:td_stim(t).idx_otHoldTime,2),'b','linewidth',2);
            hold on
            plot(td_stim(t).pos(td_stim(t).idx_stimTime,1),td_stim(t).pos(td_stim(t).idx_stimTime,2),'k.','markersize',20);
            plot(td_stim(t).pos(td_stim(t).idx_goCueTime,1),td_stim(t).pos(td_stim(t).idx_goCueTime,2),'r.','markersize',20);
        end
    end
    
%% plot velocity
    td_stim = td_all(~isnan([td_all.stimCode]));

    num_plot = 1000;
    counter = 0;
    figure();
    peak_speeds = cell(numel(unique([td_stim.stimCode])),1);
    latencies = cell(numel(unique([td_stim.stimCode])),1);
    for t = 1:numel(td_stim)

        [~,peak_speed_idx] = max(abs(td_stim(t).vel(td_stim(t).idx_stimTime+100:td_stim(t).idx_stimTime+300,2)));
        peak_speed_idx = peak_speed_idx + 100;
        peak_speeds{td_stim(t).stimCode+1}(end+1,1) = td_stim(t).vel(peak_speed_idx+td_stim(t).idx_stimTime-1,2) - td_stim(t).vel(td_stim(t).idx_stimTime,2);
    
        [~,peak_acc_idx] = max((td_stim(t).acc(td_stim(t).idx_stimTime+100:td_stim(t).idx_stimTime+300,2)));
        peak_acc_idx = peak_acc_idx + 100;
        peak_acc =  td_stim(t).acc(peak_acc_idx+td_stim(t).idx_stimTime-1,2);
        
        if(peak_speeds{td_stim(t).stimCode+1}(end) > 3)
            threshold = abs(peak_acc*0.5);
            peak_mask = 1:1:numel(td_stim(t).acc(td_stim(t).idx_stimTime:td_stim(t).idx_stimTime+300,2));
            peak_mask = peak_mask < peak_speed_idx;
            peak_mask = peak_mask';
            above_threshold = find(td_stim(t).acc(td_stim(t).idx_stimTime:td_stim(t).idx_stimTime+300,2) > threshold & peak_mask);
            if(isempty(above_threshold))
                latency = nan;
            else
                latency = above_threshold(1);
            end
            latencies{td_stim(t).stimCode+1}(end+1,1) = latency;
        else
            latency = nan;
            latencies{td_stim(t).stimCode+1}(end+1,1) = latency;
        end
        
        if(counter < num_plot && td_stim(t).stimCode == 8)
            plot(td_stim(t).acc(td_stim(t).idx_stimTime-20:td_stim(t).idx_stimTime+500,2),'b');
            hold on
            plot(td_stim(t).idx_stimTime-td_stim(t).idx_stimTime+20,td_stim(t).acc(td_stim(t).idx_stimTime,2),'k.','markersize',20);
            if(~isnan(latency))
                plot(latency+20,td_stim(t).acc(td_stim(t).idx_stimTime + latency,2),'r.','markersize',20);
            end
        end
    end
%% plot latency and speed of movement during stimulation
    figure();
    subplot(1,2,1)
    for i = 1:numel(peak_speeds)
        plot(i,peak_speeds{i},'.','markersize',12,'color','k')
        hold on
    end
    
    xlim([0,11])
    ylabel('y-vel (cm/s)')
    formatForLee(gcf)
    set(gca,'fontsize',16)
    set(gca,'XTick',[1:numel(peak_speeds)])
    set(gca,'XMinorTick','off')
    
    subplot(1,2,2)
    for i = 1:numel(peak_speeds)
        plot(i,latencies{i},'.','markersize',12,'color','k')
        hold on
    end
    
    xlim([0,11])
    ylim([0,200])
    ylabel('Latency (ms)')
    xlabel('channel')
    formatForLee(gcf)
    set(gca,'fontsize',16)
    set(gca,'XTick',[1:numel(peak_speeds)])
    set(gca,'XMinorTick','off')
    
    
    
    
%     
% %%
%     opts.FIGURE_SAVE = 0;
%     opts.FIGURE_PREFIX = 'Duncan_20190410_RR_training';
%     opts.FIGURE_DIR = folderpath;
% 
%     opts.NUM_BINS_DIR = 8;
%     opts.MAKE_FIGURES = 1;
%     opts.PLOT_POLAR = 1;
%     opts.MAX_TRIALS_PLOT = 2000;
%     opts.CIRCLE_RADIUS = 5.75;
%     opts.CIRCLE_DEPTH = 2;
% 
%     opts.DISTRIBUTION_BIN_SIZE = 5;
%     opts.BUMP_MAGS = [];
%     opts.BUMP_NONLINEARITY = 0;
% 
%     behaviorData = processBehaviorRingReporting(cds,opts);
% 
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
% %% look at force onset
%     params.event_list = {'bumpTime';'bumpDir';'otHoldTime'};
%     params.trial_results = {'R','F'};
%     td = parseFileByTrial(cds,params);
%     td = removeBadTrials(td);
% 
% %%
%     offset = [0,100];
%     middle = [-3,33];
%     tgt_dir = []; bump_dir = [];
%     acc_dir = [];
%     reach_dir = []; dist_dir = [];
%     for t = 1:300
%         bump_dir(t) = td(t).bumpDir;
%         tgt_dir(t) = td(t).target_direction*180/pi;
%         pos_start = td(t).pos(td(t).idx_bumpTime+offset(1),:);
%         pos_end = td(t).pos(td(t).idx_bumpTime+offset(2),:);
%         acc_dir(t) = atan2(pos_end(2)-pos_start(2),pos_end(1)-pos_start(1))*180/pi;
%         reach_dir(t) = atan2(td(t).pos(td(t).idx_otHoldTime,2)+middle(2),td(t).pos(td(t).idx_otHoldTime,1)+middle(1))*180/pi;
%         dist_dir(t) = sqrt((td(t).pos(td(t).idx_otHoldTime,2)+middle(2)).^2 + (td(t).pos(td(t).idx_otHoldTime,1)+middle(1)).^2);
%     end
%    
%     remove_idx = abs(tgt_dir-acc_dir)*pi/180 > 1;
%     tgt_dir(remove_idx) = []; bump_dir(remove_idx) = [];
%     acc_dir(remove_idx) = [];
%     dist_dir(remove_idx) = [];
%     reach_dir(remove_idx) = [];
% %
%     figure
%     
%     plot(tgt_dir*pi/180,acc_dir*pi/180,'.','markersize',12)
%     hold on
% %     plot(acc_dir,reach_dir,'.','markersize',12)
%     plot([-2*pi,2*pi],[-2*pi,2*pi],'k--','linewidth',2)
%     hold off
%     formatForLee(gcf)
% %     ylabel('acc dir')
% %     xlabel('bump dir')
%     xlim([-pi,pi])
%     ylim([-pi,pi])
%     set(gca,'fontsize',14)
%     
% %     
% % fit the data with a function 
% %     acc_dir(bump_dir < -170 & acc_dir > 150) = acc_dir(bump_dir < -170 & acc_dir > 150) - 180;
% % 
% %     [f,gof] = fit(acc_dir'*pi/180,bump_dir'*pi/180,'poly5');
% %     
% %     acc_test = [-pi:0.01:pi];
% %     bump_pred = feval(f,acc_test);
% %     hold on
% %     plot(bump_pred,acc_test,'r','linewidth',2)
% 
% %% do PCA on spikes during perturbation to see if nicely laid out
% 
%     idx_start = {'idx_bumpTime',0};
%     idx_end = {'idx_bumpTime',20};
%     
%     td_bump = trimTD(td_all,idx_start,idx_end);
% 
%     params_dim = [];
%     params_dim.algorithm = 'ppca';
%     params_dim.num_dims = 3;
%     [td_bump,info_out] = dimReduce(td_bump,params_dim)
%     
% 
%     %%
%     figure();
%     for tr = 1:numel(td_bump)
%         
%         if(~isnan(td_bump(tr).bumpDir))
%             color_to_use = (td_bump(tr).bumpDir + 180)/720 + [0,0,0];
%             plot(mean(td_bump(tr).LeftS1_ppca(:,2)),mean(td_bump(tr).LeftS1_ppca(:,3)),'.','color',color_to_use,'markersize',20)
%             hold on
%         end
%         
%     end
%     
%     
%     
%     
% %%
% 
%     norm_dist_dir = dist_dir;
% 
%     figure();
%     
%     plot(tgt_dir*pi/180,norm_dist_dir,'.','markersize',12)
%     hold on
%     formatForLee(gcf)
%     ylabel('Distance')
%     xlabel('Bump direction')
%     xlim([-pi,pi])
%     set(gca,'fontsize',14)
% %     f = gcf;
% %     f.Name = 'Duncan_20180904_reachVsAcc_120ms';
% %     saveFiguresLIB(f,folderpath,f.Name);
% 
% 
% % fit the data with a function 
%     [f_dist,gof] = fit(tgt_dir'*pi/180,norm_dist_dir','a*cos(b*x+c)+d');
%     
%     test = [-pi:0.01:pi];
%    
%     pred = feval(f_dist,test);
%     hold on
%     plot(test,pred,'r','linewidth',2)