%% load in a ns5

    folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\StimArtifactData\20190923_gen2gen3_testing\';
%     folderpath = 'C:\Users\jts3256\Desktop\Duncan_stim_data\getIPI\';

        
    cd(folderpath);
    file_list = dir('*.ns5');
    
    analog_pin_idx = 1;
    sync_idx = 2;
    artifact_data = {};
    sync_line_data = {};
    pwd = cd;
    window = [-4,20]; % ms
    pulse_width_1 = zeros(numel(file_list),1); % us
    pulse_width_2 = zeros(size(pulse_width_1)); % us
    interphase = 53; % us
    
    window_idx = window*30; % convert to data points

    
    for file_num = 1%:numel(file_list)
        disp(file_list(file_num).name);
        NS5 = openNSx([folderpath,file_list(file_num).name],'uV');

        artifact_data{file_num} = NS5.Data(analog_pin_idx,:);
        sync_line_data{file_num} = NS5.Data(sync_idx,:);
%         
%         % get pulse widths
%         pw1_idx = strfind(file_list(file_num).name,'PW1');
%         pw2_idx = strfind(file_list(file_num).name,'PW2');
%         amp1_idx = strfind(file_list(file_num).name,'A1');
%         amp2_idx = strfind(file_list(file_num).name,'A2');
%         stim_idx = strfind(file_list(file_num).name,'stim');
%         chan_idx = strfind(file_list(file_num).name,'chan');
%         underscore_idx = strfind(file_list(file_num).name,'_');
%         
%         pulse_width_1(file_num) = str2num(file_list(file_num).name(pw1_idx+4:underscore_idx(find(underscore_idx > pw1_idx,1,'first'))-1));
%         pulse_width_2(file_num) = str2num(file_list(file_num).name(pw2_idx+4:underscore_idx(find(underscore_idx > pw2_idx,1,'first'))-1));
%         amp_1(file_num) = str2num(file_list(file_num).name(amp1_idx+3:underscore_idx(find(underscore_idx > amp1_idx,1,'first'))-1));
%         amp_2(file_num) = str2num(file_list(file_num).name(amp2_idx+3:underscore_idx(find(underscore_idx > amp2_idx,1,'first'))-1));
%         stim_chan(file_num) = str2num(file_list(file_num).name(chan_idx+4:stim_idx-1));
    end
    cd(pwd);
    
%% pick a file (idx) and plot anodic and cathodic data
    threshold = 0.1;
    peak_data = {};
    
    for file_num = 1% :numel(file_list)
        disp(file_list(file_num).name);
        
        stim_on=find(diff(sync_line_data{file_num}-mean(sync_line_data{file_num})>3)>.5);
        
    end
    
    %%
%         stim_on=stim_on(waveforms.waveSent == 1 & ...
%             cellfun(@isequal,waveforms.chanSent,mat2cell(analog_pin_idx+zeros(size(waveforms.chanSent)),ones(size(waveforms.chanSent)))));
%         
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
        
        
                %
%         f=figure();
%         f.Name = file_list(file_num).name(1:end-10);
%         
%         subplot(2,1,1)
%         plot(x_data,(plot_data(:,1:5)'),'linewidth',1.5);
%         xlim([-4,12])
%         ylim([-5000,5000])
%         ylabel('Voltage (\muV)');
%         formatForLee(gcf)
%         xlabel('Time after stimulation offset (ms)');
%         set(gca,'fontsize',14)
%         subplot(2,1,2)
%         plot(x_data,(acausalFilter(plot_data(:,1:5))'));
%         xlim([-4,12])
%         ylim([-1000,1000])
%         xlabel('Time after stimulation offset (ms)');
%         formatForLee(gcf)
        
        
    end
    
    
%% plot amp of aux channel stim against time post stim

    colors = repmat(linspace(0,0.5,numel(peak_data))',1,3);
    
    
    for fieldname = {'cathodic','anodic'}
        f=figure(); hold on;
        f.Name = ['AmplifierRecovery_',fieldname{1}];
        for f = 1:numel(peak_data)
            ratio = peak_data{f}.(fieldname{1}).peak;
            ratio(ratio > 0) = ratio(ratio > 0)./peak_data{f}.(fieldname{1}).max_peak;
            ratio(ratio < 0) = ratio(ratio < 0)./peak_data{f}.(fieldname{1}).min_peak;
            
            
            
            plot(peak_data{f}.(fieldname{1}).t_post_stim,ratio,'color',colors(f,:))
            xlabel('Time after stimulation offset (ms)');
            ylabel('Relative gain of amplifier');
            formatForLee(gcf)
            set(gca,'fontsize',14)
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
    
%% get average response and std response over all channels tested
    channels = unique(stim_chan);
    amps = unique(amp_1);
%     t_post_stim =
%     [0.68,0.75,1.01,1.11,1.35,1.41,1.68,1.75,2.01,2.08,2.31,2.41,2.65,2.75,2.98,3.09,3.31,3.41,3.65,3.71,3.98]; % duncan
    t_post_stim = [0.45,0.71,0.81,1.04,1.14,1.38,1.45,1.68,1.78,2.013,2.11,2.35,2.45,2.68,2.78,3.01,3.11,3.41,3.65,3.75,3.98]; % han
    
    gain_ratio_cathodic = nan(numel(amps),numel(t_post_stim),numel(channels));
    gain_ratio_anodic = nan(numel(amps),numel(t_post_stim),numel(channels));
    max_time_difference = 0.1;
    
    % for each peak data
    for p = 1:numel(peak_data)
        % find amp idx and chan idx
        amp_idx = find(amps == peak_data{p}.amp1);
        chan_idx = find(channels == peak_data{p}.chan);
        if(peak_data{p}.pw1 == 200 && peak_data{p}.pw2 == 200) % make sure this is a balanced phase case
            for t = 1:numel(peak_data{p}.cathodic.t_post_stim)
                % find time in t_post_stim closest to array
                time_differences_cath = abs(peak_data{p}.cathodic.t_post_stim(t) - t_post_stim);
                [~,min_idx_cath] = min(time_differences_cath);
                if(time_differences_cath(min_idx_cath) < max_time_difference)
                    if(peak_data{p}.cathodic.peak(t) < 0)
                        gain_ratio_cathodic(amp_idx, min_idx_cath,chan_idx) = peak_data{p}.cathodic.peak(t)/peak_data{p}.cathodic.min_peak;
                    else
                        gain_ratio_cathodic(amp_idx, min_idx_cath,chan_idx) = peak_data{p}.cathodic.peak(t)/peak_data{p}.cathodic.max_peak;
                    end
                end 
            end
            
            for t = 1:numel(peak_data{p}.anodic.t_post_stim)
                time_differences_anod = abs(peak_data{p}.anodic.t_post_stim(t) - t_post_stim);
                [~,min_idx_anod] = min(time_differences_anod);
                if(time_differences_anod(min_idx_anod) < max_time_difference)
                    if(peak_data{p}.anodic.peak(t) < 0)
                        gain_ratio_anodic(amp_idx,min_idx_anod,chan_idx) = peak_data{p}.anodic.peak(t)/peak_data{p}.anodic.min_peak;
                    else
                        gain_ratio_anodic(amp_idx,min_idx_anod,chan_idx) = peak_data{p}.anodic.peak(t)/peak_data{p}.anodic.max_peak;
                    end
                end
            end
        end


    end
        
        %%
    average_gain_cathodic = mean(gain_ratio_cathodic,3,'omitnan');
    average_gain_anodic = mean(gain_ratio_anodic,3,'omitnan');
    
    
    std_gain_cathodic = std(gain_ratio_cathodic,0,3);
    x_data = t_post_stim(1):0.001:t_post_stim(end);

    fit_gain_cathodic = zeros(size(gain_ratio_cathodic,1),numel(x_data));
    fit_gain_anodic = zeros(size(gain_ratio_anodic,1),numel(x_data));
    for a = 1:size(average_gain_cathodic,1)
        cath_fit = fit(t_post_stim',average_gain_cathodic(a,:)','smoothingspline');
        anod_fit = fit(t_post_stim',average_gain_cathodic(a,:)','smoothingspline');
        fit_gain_cathodic(a,:) = feval(cath_fit,x_data);
        fit_gain_anodic(a,:) = feval(anod_fit,x_data);
    end
    
    figure();
%     plot(t_post_stim,average_gain_cathodic(10,:))
%     hold on
%     plot(t_post_stim(1):0.01:t_post_stim(end),feval(cath_fit,t_post_stim(1):0.01:t_post_stim(end)))
    ax = subplot(2,1,1);
    imagesc(x_data,amps,fit_gain_cathodic, [0,1]);
    ax.YDir = 'normal';
    colormap(inferno);
    formatForLee(gcf);
    ylabel('Amplitude (\muA)');
    c= colorbar(); c.Label.String = 'Gain';
%     figure();
    ax = subplot(2,1,2);
    imagesc(x_data,amps,fit_gain_anodic, [0,1]);
    ax.YDir = 'normal';
    colormap(inferno);
    formatForLee(gcf);
    xlabel('Time after stim offset (ms)')
    ylabel('Amplitude (\muA)');
    c= colorbar(); c.Label.String = 'Gain';
%     plot(t_post_stim,average_gain_anodic)
