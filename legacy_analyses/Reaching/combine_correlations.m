fpaths = {'D:\Lab\Data\DLC_videos\Crackle_20201203_rwFreeReach'; ...
        'D:\Lab\Data\DLC_videos\Han_20210623_rwFreeReach';...
        'D:\Lab\Data\DLC_videos\Rocket_20210723'};
    % monkey_names
    monkey_names = {'Crackle','Han','Rocket'};
    
    % get list of real PDs
    kin_corr_list = {};
    for i_file = 1:numel(fpaths)
        tuning_fname = dir([fpaths{i_file} filesep '*kinCorr*']);
        % load in tuning
        load([fpaths{i_file} filesep tuning_fname.name]);
        
        kin_corr_list{i_file} = kin_corr;
    end

    
    % plot correlations between 2D and 3D for each file
    % kin_corr is elbow-hand (vel-x, vel-y, vel-z, speed)
    
    
    markers = {'*','o','<'};
    marker_sizes = [12,6,8];
    
    figure(); hold on;
    for i=1:numel(kin_corr_list)
        kin_corr = kin_corr_list{i};
        
        % plot elbow-hand vel x
        plot(kin_corr.corr_2D(1),kin_corr.corr_3D(1),'marker',markers{i},'color',getColorFromList(1,2),...
            'markersize',marker_sizes(i),'linewidth',2)
        
        % plot elbow-hand vel y
        plot(kin_corr.corr_2D(2),kin_corr.corr_3D(2),'marker',markers{i},'color',getColorFromList(1,3),...
            'markersize',marker_sizes(i),'linewidth',2)
        
        % plot elbow-hand speed
        plot(kin_corr.corr_2D(4),kin_corr.corr_3D(4),'marker',markers{i},'color',getColorFromList(1,4),...
            'markersize',marker_sizes(i),'linewidth',2)
        
    end
    
    plot([0,1],[0,1],'k--','linewidth',1.5)
    xlabel('RT2D correlation');
    ylabel('RT3D correlation');
    
    formatForLee(gcf);
    set(gca,'fontsize',14);
    
    l=legend('Elbow-Hand Vel-X','Elbow-Hand Vel-Y', 'Elbow-Hand Speed');
    set(l,'box','off','location','best')
    
    