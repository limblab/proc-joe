%% make 3D PD distribution plot across monkeys, compare to vel and spindle pred
    fpaths = 'D:\Lab\Data\FreeReaching\summary_data';
    
    
    % get list of real PDs
    pd_list = [];
    tuning_fnames = dir([fpaths filesep '*S1_pd*']);
    for i_file = 1:numel(tuning_fnames)
        % load in tuning
        load([fpaths filesep tuning_fnames(i_file).name]);

        %2D PD, 3D PD
        pd_list = [pd_list; tuning_data.velPD_2D, tuning_data.velPD_zAng_2D, tuning_data.velPD_3D,tuning_data.velPD_zAng_3D];
    end

    % get hand vel distribution
    handvel_list = {};
    handvel_list{1} = []; handvel_list{2} = [];
    tuning_fnames = dir([fpaths filesep '*handvel_pd*']);
    for i_file = 1:numel(tuning_fnames)
        
        % load in tuning
        load([fpaths filesep tuning_fnames(i_file).name]);
        
        % avg 2 tuning is 3D task; planar PD, z angle
        handvel_list{1} = [handvel_list{1}; tuning_data.velPD_2D, tuning_data.velPD_zAng_2D]; 
        handvel_list{2} = [handvel_list{2}; tuning_data.velPD_3D,tuning_data.velPD_zAng_3D];
    end
    
    % get list of spindle PDs
    spindle_list = [];
    tuning_fnames = dir([fpaths filesep '*spindle_pd*']);
    for i_file = 1:numel(tuning_fnames)
        
        % load in tuning
        load([fpaths filesep tuning_fnames(i_file).name]);

        % avg 2 tuning is 3D task; planar PD, z angle
        spindle_list = [spindle_list; tuning_data.velPD_2D, tuning_data.velPD_zAng_2D, tuning_data.velPD_3D, tuning_data.velPD_zAng_3D];
    end
    
%% 3D representation of PD, scatter between z_ang and planar PD
    s1_pd_2D = [cos(pd_list(:,2)).*cos(pd_list(:,1)), cos(pd_list(:,2)).*sin(pd_list(:,1)), sin(pd_list(:,2))];
    s1_pd_3D = [cos(pd_list(:,4)).*cos(pd_list(:,3)), cos(pd_list(:,4)).*sin(pd_list(:,3)), sin(pd_list(:,4))];
    f=figure('Position',[680 558 853 420]);
    ax1=subplot(1,2,1); hold on;
    for i = 1:length(s1_pd_3D)
        plot3([0,s1_pd_3D(i,1)],[0,s1_pd_3D(i,2)],[0,s1_pd_3D(i,3)],'k');
    end
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    ax1.CameraPosition = [-9.3579 -13.2656 6.0376];
    xlim([-1,1]); ylim([-1,1]); zlim([-1,1]);
    axis square
    
    ax2=subplot(1,2,2);
    plot(180/pi*pd_list(:,3),180/pi*pd_list(:,4),'k.','markersize',12)
    
    xlim([-180,180]);
    ylim([-90,90]);
    formatForLee(gca);
    xlabel('Planar PD (deg)');
    ylabel('Z-ang (deg)');
    
%% rose plot of actual planar PD compared to hand vel and spindle predictions. Same for z_ang
    planar_nbins = 20; 
    figure('Position',[2294 267 1102 611]); % 3D task
    
    rlim_max= 0.3;
    subplot(2,3,1)
    polarhistogram(pd_list(:,3),planar_nbins,'Normalization','probability');
    rlim([0,rlim_max])
    subplot(2,3,2)
    polarhistogram(handvel_list{2}(:,1),planar_nbins,'Normalization','probability');
    rlim([0,rlim_max])
    subplot(2,3,3)
    polarhistogram(spindle_list(:,3),planar_nbins,'Normalization','probability');
    rlim([0,rlim_max])

    
    z_ang_bins = linspace(-90,90,20);
    ax1=subplot(2,3,4)
    histogram(180/pi*pd_list(:,4),z_ang_bins,'Normalization','probability');
    xlabel('Z-ang (deg)');
    ylabel('Number of neurons');
    formatForLee(gcf);
    ax2=subplot(2,3,5)
    histogram(180/pi*handvel_list{2}(:,2),z_ang_bins,'Normalization','probability');
    formatForLee(gcf);
    ax3=subplot(2,3,6)
    histogram(180/pi*spindle_list(:,4),z_ang_bins,'Normalization','probability');
    formatForLee(gcf);
    linkaxes([ax1,ax2,ax3],'xy');
    
%% compare 2D and 3D PDs
    
    

    figure('Position',[2027 456 872 420]);
    ax1=subplot(1,2,1); hold on;
    pd_2d = 180/pi*pd_list(:,1);
    pd_3d = 180/pi*pd_list(:,3);
    pd_3d_adj = pd_2d + angleDiff(pd_2d,pd_3d,0,1);
    
    plot(pd_2d,pd_3d_adj,'ko','linewidth',1.25)
    plot([-180,180],[-180,180],'k--');
    plot([-180,180],[-270,90],'k:');
    plot([-180,180],[-90,270],'k:');
    pbaspect([0.5,1,1])
    
    
    formatForLee(ax1)
    xlabel('RT2D planar PD (deg)');
    ylabel('RT3D planr PD (deg)');
    xlim([-180,180]);
    ylim([-360,360])
    
    ax1=subplot(1,2,2)
    hist_bins = -180:10:180;
    histogram(angleDiff(pd_2d,pd_3d,0,1),hist_bins,'Normalization','Probability')
    xlim([-180,180]);
    xlabel('Planar PD_2_D - Planar PD_3_D (deg)');
    ylabel('Proportion of neurons');
    formatForLee(ax1);