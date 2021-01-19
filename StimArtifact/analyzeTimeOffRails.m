%% load in a ns5


    folderpath = 'C:\Users\jts3256\Desktop\Duncan_Han_blackrockAmpRecovery\';
        
    cd(folderpath);
    file_list = dir('*.ns5');
    
    sync_idx = 97;
    artifact_data = {};
    sync_line_data = {};
    pwd = cd;
    window = [-8,60]; % ms
    pulse_width_1 = zeros(numel(file_list),1); % us
    pulse_width_2 = zeros(size(pulse_width_1)); % us
    interphase = 53; % us
    
    window_idx = window*30; % convert to data points

    
    for file_num = 1:numel(file_list)
       
        % get pulse widths
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
        
        
        
        NS5 = openNSx([folderpath,file_list(file_num).name],'uV');
    
        artifact_data{file_num} = NS5.Data(stim_chan(file_num),:);
        sync_line_data{file_num} = NS5.Data(end,:);
    end
    cd(pwd);
    
%% pick a file (idx) and plot anodic and cathodic data
    time_off_rails_data = [];
    max_amp = 30;
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
        
        cathodic_data = acausalFilter(mean(plot_data(:,cathodic_idx)'));
        anodic_data = acausalFilter(mean(plot_data(:,anodic_idx)'));       

        temp_cathodic = find(abs(cathodic_data(find(x_data >= 0,1,'first')+5:end)) > max_amp,1,'last');
        temp_anodic = find(abs(anodic_data(find(x_data >= 0,1,'first')+5:end)) > max_amp,1,'last');
        
        time_off_rails_data(file_num,:) = [temp_cathodic,temp_anodic] +5;
        
        
%         f=figure();
%         f.Name = [file_list(file_num).name(1:end-10),'_raw'];
%         
%         plot(x_data,(plot_data(:,1:2:6)'),'color',getColorFromList(1,0),'linewidth',1.5);
%         hold on 
%         plot(x_data,(plot_data(:,2:2:6)'),'color',getColorFromList(1,1),'linewidth',1.5);
%         xlim([-4,12])
%         ylim([-9000,9000])
%         ylabel('Voltage (\muV)');
%         formatForLee(gcf)
%         xlabel('Time after stimulation offset (ms)');
%         set(gca,'fontsize',14)

%                 %
%             f=figure();
%             f.Name = file_list(file_num).name(1:end-10);
%     %         
%             plot(x_data,(cathodic_data'),'r','linewidth',1.5);
%             hold on
%             plot(x_data,(anodic_data'),'b','linewidth',1.5);
%             xlim([-4,10])
%             ylim([-100,100])
%             ylabel('Voltage (\muV)');
%             formatForLee(gcf)
%             xlabel('Time after stimulation offset (ms)');
%             set(gca,'fontsize',14)
%             hold on
%             try
%                 plot(x_data(temp_cathodic+find(x_data >= 0,1,'first')+5)+[0,0],[-9000,9000],'--','color',[0,0.5,0],'linewidth',4);
%                 hold on
%                 plot(x_data(temp_anodic+find(x_data >= 0,1,'first')+5)+[0,0],[-9000,9000],'k--','linewidth',4);
%             catch
%             end
%         pause;
%         close all
    end
    
    
%% get average response and std response over all channels tested and plot
    channels = unique(stim_chan);
    amps = unique(amp_1);

    time_off_rails_cathodic = nan(numel(amps),numel(channels));
    time_off_rails_anodic = nan(numel(amps),numel(channels));
    
    % for each file
    for r = 1:size(time_off_rails_data,1)
        % find amp idx and chan idx
        amp_idx = find(amps == amp_1(r));
        chan_idx = find(channels == stim_chan(r));
        if(pulse_width_1(r) == 200 && pulse_width_2(r) == 200) % make sure this is a balanced phase case
            time_off_rails_cathodic(amp_idx, chan_idx) = time_off_rails_data(r,1)/30;
            time_off_rails_anodic(amp_idx, chan_idx) = time_off_rails_data(r,2)/30;
        end


    end
        
    mean_time_cathodic = mean(time_off_rails_cathodic,2,'omitnan');
    mean_time_anodic = mean(time_off_rails_anodic,2,'omitnan');
    
    
    std_time_cathodic = std(time_off_rails_cathodic,0,2,'omitnan');
    std_time_anodic = std(time_off_rails_anodic,0,2,'omitnan');
    
    x_data = amps;
    
    f=figure();
%     boxplot(time_off_rails_cathodic');
    errorbar(amps,mean_time_cathodic,std_time_cathodic,'color',getColorFromList(1,0),'marker','.','markersize',20,'linewidth',1.5);
    hold on
    errorbar(amps,mean_time_anodic,std_time_anodic,'color',getColorFromList(1,1),'marker','.','markersize',20,'linewidth',1.5);
%     
%     xlabel('Amplitude (\muA)');
%     ylabel('Mean time off rails post stim (ms)');
%     formatForLee(gcf);
%     set(gca,'fontsize',14);
%     l=legend('Cathodic','Anodic');
%     set(l,'box','off','fontsize',12,'location','best')
    
    
    
    
    
    