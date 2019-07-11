%% define a population of exponential decreasing curves 
% RT = a*exp(-b*I)+c
    a_bounds = [0.12,0.23];
    b_bounds = [0.001,0.05];
    c_bounds = [0.15,0.25];
    std_bounds = [0.01,0.05];
    
    min_I_bounds = [0,20];
    num_chans = 96;
    a_all = rand(num_chans,1)*(diff(a_bounds)) + a_bounds(1);
    b_all = rand(num_chans,1)*(diff(b_bounds)) + b_bounds(1);
    c_all = rand(num_chans,1)*(diff(c_bounds)) + c_bounds(1);
    min_I_all = rand(num_chans,1)*diff(min_I_bounds) + min_I_bounds(1);
    std_all = rand(num_chans,1)*diff(std_bounds) + std_bounds(1);
%% plot example RT vs amp curves (only for 2 electrodes)
    num_runs_per_chan = 100;

    figure();
    hold on;
    I_data = 1:1:100;
    RT_all_means = a_all.*exp(-b_all.*I_data)+c_all;  
%     RT_all_dots = repmat(RT_all_means,1,1,num_runs_per_chan) + normrnd(0,1,[size(RT_all_means),num_runs_per_chan]).*repmat(std_all,1,1,num_runs_per_chan);
    
    % std deviation bars
    color_offset=[2,2,1];
    for chan = 1:size(RT_all_means,1)
        x_data = [I_data(I_data>min_I_all(chan)),fliplr(I_data(I_data>min_I_all(chan)))];
        y_data = [RT_all_means(chan,I_data>min_I_all(chan))+std_all(chan,1),fliplr(RT_all_means(chan,I_data>min_I_all(chan))-std_all(chan,1))];
           
        fill(x_data,y_data,getColorFromList(1,chan+color_offset(chan)),'EdgeColor','none');
        alpha(0.5);
    end
    % fit lines
    for chan = 1:size(RT_all_means,1)
        
        plot(I_data(I_data>min_I_all(chan)),RT_all_means(chan,I_data>min_I_all(chan)),'linewidth',2,'linestyle','-','color',getColorFromList(1,chan+color_offset(chan)))

    end
    
    formatForLee(gcf)
    xlabel('Amplitude (\muA)');
    ylabel('RT (s)');
    set(gca,'fontsize',16);
    
% %     % all dots
% %     for chan = 1:size(RT_all_means,1)
% %         for i = 1:numel(I_data_plot)
% %             I_idx = find(I_data == I_data_plot(i));
% %             scatter(repmat(I_data_plot(i)+chan_offset(chan),num_runs_per_chan,1),squeeze(RT_all_dots(chan,I_idx,:)),10,'markeredgecolor','none','markerfacecolor',getColorFromList(1,chan));
% %             alpha(0.5);
% %         end
% %     end
% %     
% %     % means
% %     for chan = 1:size(RT_all_means,1)
% %         for i = 1:numel(I_data_plot)
% %             I_idx = find(I_data == I_data_plot(i));
% %             plot(I_data_plot(i)+chan_offset(chan),RT_all_means(chan,I_idx),'.','markersize',24,'color',getColorFromList(1,chan));
% %         end
% %     end
%% run my linear summation experiment
% sample N electrodes, compute RT based on the race model (sample each
% distribution, pick fastest). Store. Do this for different charges on each
% electrode and for different number of electrodes

    I_max = 100;
    num_elecs = [1:5:30];%3,6,12,24];
    total_charge = [360:120:1200];
    num_runs_per_condition = 1000;

    RT_out = zeros(numel(num_elecs),numel(total_charge),num_runs_per_condition);
    for e = 1:numel(num_elecs)
        for c = 1:numel(total_charge)
            % sample num_elecs(e)
            chan_idx = [];
            for n = 1:num_runs_per_condition
                chan_idx(n,:) = randperm(num_chans,num_elecs(e));
            end
            % sample RT from normal distribution based on channel parameters
            a_sample = a_all(chan_idx);
            b_sample = b_all(chan_idx);
            c_sample = c_all(chan_idx);
            I_sample = total_charge(c)/num_elecs(e);
            min_I_sample = min_I_all(chan_idx);
            std_sample = std_all(chan_idx);%min(1,(I_sample)/I_max)*(std_min_all(chan_idx)-std_max_all(chan_idx)) + std_max_all(chan_idx);

            RT_all_means = a_sample.*exp(-b_sample.*I_sample)+c_sample +randn(size(std_sample)).*std_sample;
            RT_all_means(min_I_sample > I_sample) = nan;
            % store min RT for each combo
            RT_out(e,c,:) = min(RT_all_means,[],2);
        end
    end

% plot RT_out data
    mean_RT = mean(RT_out,3,'omitnan');
    std_dev_RT = std(RT_out,[],3,'omitnan');
    colors = [228,26,28;55,126,184;77,175,74;152,78,163;255,127,0;255,255,51;166,86,40;247,129,191;153,153,153]/255;
    
    f=figure();
    hold on;
    for c=1:numel(total_charge)
        plot((num_elecs),squeeze(mean_RT(:,c)),'.','markersize',24,'color',colors(c,:))
    end
%     for c=1:numel(total_charge)
%         plot((repmat(num_elecs',1,2)'),(squeeze(mean_RT(:,c,:))+squeeze(std_dev_RT(:,c,:)).*[-1,1])','linewidth',1.5,'color',colors(c,:))
%     end
    

%% toy model to show normal distribution shift concept

num_samples = 20000000;
means = [0.2,0.21];
std_dev = [0.04,0.05];

data = normrnd(0,1,[num_samples,2]).*std_dev + means;
data(:,3) = min(data,[],2);
bE = [0:0.003:0.5];
bC = [];
%%

for i = 1:size(data,2)
    figure();
    hold on
    bC(i,:) = histcounts(data(:,i),bE)/num_samples;
    plot(bE(1:end-1)+mode(diff(bE))/2,bC(i,:),'linewidth',2,'color',getColorFromList(1,i))
    [~,max_idx] = max(bC(i,:));
    plot(bE(max_idx)+mode(diff(bE))/2+[0,0],[0,bC(i,max_idx)],'linewidth',2,'color',getColorFromList(1,i))
    
    
    formatForLee(gcf)
    xlabel('RT')
    ylabel('Density')
    ax = gca;
    ax.YTick = [];
    ax.YTickLabel = {};
    ax.YMinorTick = 'off';
    ax.XMinorTick = 'off';
    ax.XTick = [];
    ylim([0,0.035])
    xlim([0,0.4])
    ax.FontSize = 16;
    f = gcf;
    f.Units = 'inches';
    f.Position = [8.0104 5.6563 5.3646 3.5521];
    f.Name = ['RT_independenceModel_normalDistributions_chan',num2str(i)];

    fpath = 'C:\Users\jts3256\Desktop\Han_20180826_FCreactTime\';
    saveFiguresLIB(f,fpath,f.Name);
end

