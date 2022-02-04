% combine previously made duke and blackrock artifact figures
f_duke = openfig('D:\Lab\Data\stim_ephys_paper\Raw_Figures\FIG\noAuxChan61__chan4stim_A1-50_A2-50_PW1-200_PW2-200_interpulse53_po_raw.fig');
f_black = openfig('D:\Lab\Data\stim_ephys_paper\Raw_Figures\FIG\noAuxChan55_fs__chan4stim_A1-50_A2-50_PW1-200_PW2-200_interpulse53_po_raw.fig');


fig_keep = figure();
for i_sub = 1:2
    ax_keep=subplot(2,1,i_sub);
    switch i_sub
        case 1
            fig_move = f_black;
            line_style = '--';
        case 2 
            fig_move = f_duke;
            line_style = '-';
    end
    line_move = findobj(fig_move,'type','line');

    for l = 1:numel(line_move)
        new_handle = copyobj(line_move(l),ax_keep);
        new_handle.XData = line_move(l).XData;
        new_handle.LineStyle = line_style;
    end
    
    xlim([-4,12]);
    formatForLee(gcf);
    ylim([-8500,8500]);
    set(gca,'fontsize',14);
    ylabel('Voltage (\muV)');
    if(i_sub==2)
        xlabel('Time after stimulation offset (ms)');
    end
    
end


%%
fig_keep.Name = 'Han_duncan_artifactExamples_BlackrockDashed_DukeSolid_CathodicRed_AnodicBlue_50uA_chan4stimrec';

fpath = 'D:\Lab\Data\stim_ephys_paper\';
saveFiguresLIB(fig_keep,fpath,fig_keep.Name);

%%

f_duke = openfig('D:\Lab\Data\stim_ephys_paper\Raw_Figures\FIG\noAuxChan61__chan4stim_A1-50_A2-50_PW1-200_PW2-200_interpulse53_po_raw.fig');

ax = f_duke.Children;
line_list = ax.Children;

y_data_filt = zeros(numel(line_list),numel(line_list(1).XData));
for l = 1:numel(line_list)
    x_data = line_list(l).XData;
    y_data = line_list(l).YData;
    
    y_data_filt(l,:) = acausalFilter(y_data);
end

% anodic first is 1:3, cathodic first is 4:6

f=figure();
subplot(2,1,1) % cathodic first
plot(x_data,y_data_filt(4:6,:),'color',getColorFromList(1,0),'linewidth',1.25)

xlim([-4 12])
ylim([-500,500])
formatForLee(gcf);
set(gca,'fontsize',14)
ylabel('Voltage (\muV)')

subplot(2,1,2) % anodic first
plot(x_data,y_data_filt(1:3,:),'color',getColorFromList(1,1),'linewidth',1.25)

xlim([-4 12])
ylim([-500,500])
formatForLee(gcf);
set(gca,'fontsize',14)
xlabel('Time after stimulation offset (ms)');
ylabel('Voltage (\muV)')