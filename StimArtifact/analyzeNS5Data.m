%% load in a ns5

    folderpath = 'C:\Users\jts3256\Desktop\Han_20190730_trains_differentIPIs\chan21\\';
    analog_pin_idx = 1;
    sync_idx = 2;
    artifact_data = {};
    pwd = cd;
    window = [-2,8]; % ms
    pulse_width = 200; % us
    interphase = 53; % us
    
    window_idx = window*30; % convert to data points
    
    cd(folderpath);
    file_list = dir('*.ns5');
    
    for file_num = 1%:numel(file_list)
        NS5 = openNSx([folderpath,file_list(file_num).name],'uV');
    
%         artifact_data{file_num} = NS5.Data(analog_pin_idx,:);
%         sync_line_data{file_num} = NS5.Data(sync_idx,:);
    end
    
    %% pick a file (idx) and plot anodic and cathodic data
    
    for file_num = 1:2
        disp(file_list(file_num).name);

        stim_on=find(diff(sync_line_data{file_num}-mean(sync_line_data{file_num})>3)>.5);

        anodic_idx = 2:2:10;%numel(stim_on);
        cathodic_idx = 1:2:10;%numel(stim_on);

        plot_data = [];
        for st = 1:numel(stim_on)
            plot_data(:,st) = artifact_data{file_num}((stim_on(st)+window_idx(1)):(stim_on(st)+window_idx(2)));
        end

        %
        x_data = [window_idx(1):window_idx(2)]'/30 - (2*pulse_width + interphase)/1000;
        f=figure();
        f.Name = file_list(file_num).name(1:end-10);
        
        subplot(2,1,1)
        plot(x_data,mean(plot_data(:,cathodic_idx)'));
        xlim([-2,4])
        ylim([-5000,5000])
        ylabel('Voltage (\muV)');
        formatForLee(gcf)

        subplot(2,1,2)
        plot(x_data,mean(acausalFilter(plot_data(:,cathodic_idx))'));
        xlim([-2,4])
        ylim([-500,500])
        xlabel('Time after stimulation offset (ms)');
        formatForLee(gcf)
        
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    