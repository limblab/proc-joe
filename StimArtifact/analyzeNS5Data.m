%% load in a ns5


    folderpath = 'C:\Users\jts3256\Desktop\Han_stim_data\Han_20190821_dukeProjBox\chan39\';
        
    cd(folderpath);
    file_list = dir('*.ns5');
    
    analog_pin_idx = 97;
    sync_idx = 98;
    artifact_data = {};
    pwd = cd;
    window = [-2,8]; % ms
    pulse_width_1 = zeros(numel(file_list),1); % us
    pulse_width_2 = zeros(size(pulse_width_1)); % us
    interphase = 53; % us
    
    window_idx = window*30; % convert to data points

    
    for file_num = 1:numel(file_list)

        NS5 = openNSx([folderpath,file_list(file_num).name],'uV');
    
        artifact_data{file_num} = NS5.Data(analog_pin_idx,:);
        sync_line_data{file_num} = NS5.Data(sync_idx,:);
        
        % get pulse widths
        pw1_idx = strfind(file_list(file_num).name,'PW1');
        pw2_idx = strfind(file_list(file_num).name,'PW2');
        amp1_idx = strfind(file_list(file_num).name,'A1');
        amp2_idx = strfind(file_list(file_num).name,'A2');
        underscore_idx = strfind(file_list(file_num).name,'_');
        
        pulse_width_1(file_num) = str2num(file_list(file_num).name(pw1_idx+4:underscore_idx(find(underscore_idx > pw1_idx,1,'first'))-1));
        pulse_width_2(file_num) = str2num(file_list(file_num).name(pw2_idx+4:underscore_idx(find(underscore_idx > pw2_idx,1,'first'))-1));
        amp_1(file_num) = str2num(file_list(file_num).name(amp1_idx+3:underscore_idx(find(underscore_idx > amp1_idx,1,'first'))-1));
        amp_2(file_num) = str2num(file_list(file_num).name(amp2_idx+3:underscore_idx(find(underscore_idx > amp2_idx,1,'first'))-1));
    end
    cd(pwd);
    
%% pick a file (idx) and plot anodic and cathodic data
    
    peak_data = {};
    
    for file_num = 1:numel(file_list)
%         disp(file_list(file_num).name);
        
        stim_on=find(diff(sync_line_data{file_num}-mean(sync_line_data{file_num})>3)>.5);

        anodic_idx = 2:2:numel(stim_on);
        cathodic_idx = 1:2:numel(stim_on);

        plot_data = [];
        x_data = [window_idx(1):window_idx(2)]'/30 - (pulse_width_1(file_num) + pulse_width_2(file_num) + interphase)/1000;

        for st = 1:numel(stim_on)
            plot_data(:,st) = artifact_data{file_num}((stim_on(st)+window_idx(1)):(stim_on(st)+window_idx(2)));
        end
        
        %
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
        ylim([-1500,1500])
        xlabel('Time after stimulation offset (ms)');
        formatForLee(gcf)
        
        filtered_cathodic_data = mean(acausalFilter(plot_data(:,cathodic_idx))');
        filtered_anodic_data = mean(acausalFilter(plot_data(:,anodic_idx))');
        
        max_cathodic_peak = max(filtered_cathodic_data(find(x_data > 3,1,'first'):find(x_data > 6,1,'first')));
        min_cathodic_peak = min(filtered_cathodic_data(find(x_data > 3,1,'first'):find(x_data > 6,1,'first')));
        
        max_anodic_peak = max(filtered_anodic_data(find(x_data > 3,1,'first'):find(x_data > 6,1,'first')));
        min_anodic_peak = min(filtered_cathodic_data(find(x_data > 3,1,'first'):find(x_data > 6,1,'first')));
        
        [pks,locs,~,~] = findpeaks(abs(filtered_cathodic_data),'MinPeakHeight',max_cathodic_peak*0.2);
        locs(pks < max_cathodic_peak*0.2) = [];
        locs(x_data(locs) < 0.4 | x_data(locs) > 4) = [];
        
        peak_data{file_num}.cathodic.t_post_stim = x_data(locs);
        peak_data{file_num}.cathodic.peak = filtered_cathodic_data(locs);
        peak_data{file_num}.cathodic.max_peak = max_cathodic_peak;
        peak_data{file_num}.cathodic.min_peak = min_cathodic_peak;

        
        [pks,locs,~,~] = findpeaks(abs(filtered_anodic_data),'MinPeakHeight',max_anodic_peak*0.2);
        locs(pks < max_anodic_peak*0.2) = [];
        locs(x_data(locs) < 0.4 | x_data(locs) > 4) = [];
        
        peak_data{file_num}.anodic.t_post_stim = x_data(locs);
        peak_data{file_num}.anodic.peak = filtered_anodic_data(locs);
        peak_data{file_num}.anodic.max_peak = max_anodic_peak;
        peak_data{file_num}.anodic.min_peak = min_anodic_peak;
        
        peak_data{file_num}.amp1 = amp_1(file_num);
        peak_data{file_num}.amp2 = amp_2(file_num);
        peak_data{file_num}.pw1 = pulse_width_1(file_num);
        peak_data{file_num}.pw2 = pulse_width_2(file_num);
    end
    
    
%% plot amp of aux channel stim against time post stim

    colors = repmat(linspace(0,0.5,numel(peak_data))',1,3);
    
    
    for fieldname = {'cathodic','anodic'}
        figure(); hold on;
        for f = 1:numel(peak_data)
            ratio = peak_data{f}.(fieldname{1}).peak;
            ratio(ratio > 0) = ratio(ratio > 0)./peak_data{f}.(fieldname{1}).max_peak;
            ratio(ratio < 0) = ratio(ratio < 0)./peak_data{f}.(fieldname{1}).min_peak;
            
            
            
            plot(peak_data{f}.(fieldname{1}).t_post_stim,ratio,'color',colors(f,:))
        end
        
    end
    
%% make model to predict recovery as a function of time post stimulation, amplitude and polarity
    X = []; % n x m matrix
    Y = []; % n x 1 matrix

    for fieldname = {'cathodic'}
        if(strcmpi(fieldname,'cathodic'))
            is_cathodic = 1;
        else
            is_cathodic = 0;
        end
        for f = 1:numel(peak_data)
            ratio = peak_data{f}.(fieldname{1}).peak;
            ratio(ratio > 0) = ratio(ratio > 0)./peak_data{f}.(fieldname{1}).max_peak;
            ratio(ratio < 0) = ratio(ratio < 0)./peak_data{f}.(fieldname{1}).min_peak;
            t_post_stim = peak_data{f}.(fieldname{1}).t_post_stim;
            
            X = [X; t_post_stim, peak_data{f}.amp1 + zeros(size(t_post_stim)), peak_data{f}.pw1 + zeros(size(t_post_stim))];%, is_cathodic + zeros(size(t_post_stim))];
            Y = [Y;ratio'];
            
        end
        
    end
    

    % glm
    link_func = 'loglog';
    [b,dev,stats] = glmfit(X,Y,'poisson','link',link_func);
    y_fit = glmval(b,X,link_func);
    max_time_post_stim = max(t_post_stim);
    
    figure();
    plot(Y);
    hold on;
    plot(real(y_fit));
    
    