%% load in a ns5

    folderpath = 'E:\Data\Joseph\Han_stim_data\Han_20191210_dukeBoards_lowpassFilters\noaux\';
%     folderpath = 'C:\Users\jts3256\Desktop\Duncan_stim_data\getIPI\';

        
    cd(folderpath);
    file_list = dir('*.ns5');
    
    analog_pin_idx = 1;
    sync_idx = 2;
    artifact_data = {};
    sync_line_data = {};
    pwd = cd;
    window = [-4,50]; % ms
    pulse_width_1 = zeros(numel(file_list),1); % us
    pulse_width_2 = zeros(size(pulse_width_1)); % us
    interphase = 53; % us
    
    window_idx = window*30; % convert to data points

%     f= figure();
    condition_counter = 1;
    for file_num = 1:numel(file_list)
        disp(file_list(file_num).name);
        NS5 = openNSx([folderpath,file_list(file_num).name],'uV');
    
% 
        artifact_data{file_num} = NS5.Data(analog_pin_idx,:);
        sync_line_data{file_num} = NS5.Data(sync_idx,:);
        stim_on=find(diff(sync_line_data{file_num}-mean(sync_line_data{file_num})>3)>.5);
        
        pre_stim = artifact_data{file_num}(1:stim_on(1)-1200);
        during_stim = artifact_data{file_num}(stim_on(1)-1200:stim_on(end)+300);
        stim_on = stim_on - stim_on(1) + 1201;
        x_data = (0:1:(numel(pre_stim)-1))/30;
        f=figure();
        f.Name = file_list(file_num).name;
        plot(x_data,pre_stim);
        xlim([0,15])
        ylim([-200,200])
        xlabel('Time (ms)');
        ylabel('Voltage (uV)');
%         f=figure();
%         f.Name = file_list(file_num).name;
%         periodogram(pre_stim)
%         subplot(2,2,condition_counter)
%         condition_counter = condition_counter + 1;
%         
%         periodogram(pre_stim);
%         ylim([-40,100])
%         title('')
        
%         
%         % get pulse widths
        pw1_idx = strfind(file_list(file_num).name,'PW1');
        pw2_idx = strfind(file_list(file_num).name,'PW2');
        amp1_idx = strfind(file_list(file_num).name,'A1');
        amp2_idx = strfind(file_list(file_num).name,'A2');
        stim_idx = strfind(file_list(file_num).name,'stim');
        chan_idx = strfind(file_list(file_num).name,'chan');
        underscore_idx = strfind(file_list(file_num).name,'_');
        
        pulse_width_1(file_num) = str2num(file_list(file_num).name(pw1_idx+4:underscore_idx(find(underscore_idx > pw1_idx,1,'first'))-1));
        pulse_width_2(file_num) = str2num(file_list(file_num).name(pw2_idx+4:underscore_idx(find(underscore_idx > pw2_idx,1,'first'))-1));
        amp_1(file_num) = str2num(file_list(file_num).name(amp1_idx+3:underscore_idx(find(underscore_idx > amp1_idx,1,'first'))-1));
        amp_2(file_num) = str2num(file_list(file_num).name(amp2_idx+3:underscore_idx(find(underscore_idx > amp2_idx,1,'first'))-1));
        stim_chan(file_num) = str2num(file_list(file_num).name(chan_idx+4:stim_idx-1));
    end
    cd(pwd);
    
%% pick a file (idx) and plot anodic and cathodic data
    make_plot = 0;
    threshold = 0.1;
    peak_data = cell(numel(file_list),1);
    
    condition_counter = 1;
    
    for file_num = 1:numel(file_list)
        disp(file_list(file_num).name);
        
        stim_on=find(diff(sync_line_data{file_num}-mean(sync_line_data{file_num})>3)>.5);
            
        anodic_idx = 2:2:numel(stim_on);
        cathodic_idx = 1:2:numel(stim_on);

        plot_data = [];
        x_data = [window_idx(1):window_idx(2)]'/30 - (pulse_width_1(file_num) + pulse_width_2(file_num) + interphase)/1000;

        for st = 1:numel(stim_on)
            plot_data(:,st) = artifact_data{file_num}((stim_on(st)+window_idx(1)):(stim_on(st)+window_idx(2)));
        end
        
        filtered_cathodic_data = mean(acausalFilter(plot_data(:,cathodic_idx))');
        filtered_anodic_data = mean(acausalFilter(plot_data(:,anodic_idx))');
        
        max_cathodic_peak = max(filtered_cathodic_data(find(x_data > 4,1,'first'):find(x_data > 7,1,'first')));
        min_cathodic_peak = min(filtered_cathodic_data(find(x_data > 4,1,'first'):find(x_data > 7,1,'first')));
        
        max_anodic_peak = max(filtered_anodic_data(find(x_data > 4,1,'first'):find(x_data > 7,1,'first')));
        min_anodic_peak = min(filtered_cathodic_data(find(x_data > 4,1,'first'):find(x_data > 7,1,'first')));
        
        [pks,locs,~,~] = findpeaks(abs(filtered_cathodic_data),'MinPeakHeight',max_cathodic_peak*threshold);
        locs(pks < max_cathodic_peak*threshold) = [];
        locs(x_data(locs) < 0.4 | x_data(locs) > 4) = [];
        
        peak_data{file_num}.cathodic.t_post_stim = x_data(locs);
        peak_data{file_num}.cathodic.peak = filtered_cathodic_data(locs);
        peak_data{file_num}.cathodic.max_peak = max_cathodic_peak;
        peak_data{file_num}.cathodic.min_peak = min_cathodic_peak;

        
        [pks,locs,~,~] = findpeaks(abs(filtered_anodic_data),'MinPeakHeight',max_anodic_peak*threshold);
        locs(pks < max_anodic_peak*threshold) = [];
        locs(x_data(locs) < 0.4 | x_data(locs) > 4) = [];
        
        peak_data{file_num}.anodic.t_post_stim = x_data(locs);
        peak_data{file_num}.anodic.peak = filtered_anodic_data(locs);
        peak_data{file_num}.anodic.max_peak = max_anodic_peak;
        peak_data{file_num}.anodic.min_peak = min_anodic_peak;
        
        peak_data{file_num}.amp1 = amp_1(file_num);
        peak_data{file_num}.amp2 = amp_2(file_num);
        peak_data{file_num}.pw1 = pulse_width_1(file_num);
        peak_data{file_num}.pw2 = pulse_width_2(file_num);
        peak_data{file_num}.chan = stim_chan(file_num);
        
        if(make_plot)
            f=figure();
            f.Name = file_list(file_num).name(1:end-10);

            subplot(2,1,1)
            plot(x_data,(plot_data(:,1:2:6)'),'linewidth',1,'color',getColorFromList(1,0));
            hold on
            plot(x_data,(plot_data(:,2:2:6)'),'linewidth',1,'color',getColorFromList(1,1));
            xlim([-4,12])
            ylim([-5000,5000])
            ylabel('Voltage (\muV)');
            formatForLee(gcf)
            xlabel('Time after stimulation offset (ms)');
            set(gca,'fontsize',14)
            subplot(2,1,2)
            plot(x_data,(acausalFilter(plot_data(:,1:2:6))'),'color',getColorFromList(1,0));
            hold on
            plot(x_data,(acausalFilter(plot_data(:,2:2:6))'),'color',getColorFromList(1,1));
            xlim([-4,12])
            ylim([-1000,1000])

            xlabel('Time after stimulation offset (ms)');
            formatForLee(gcf)
        end
%         
        
    end
    
