% combine previously made duke and blackrock artifact figures
f_duke = openfig('D:\Lab\Data\stim_ephys_paper\Raw_Figures\FIG\noAuxChan61__chan4stim_A1-50_A2-50_PW1-200_PW2-200_interpulse53_po_raw.fig');
f_black = openfig('D:\Lab\Data\stim_ephys_paper\Raw_Figures\FIG\noAuxChan55_fs__chan4stim_A1-50_A2-50_PW1-200_PW2-200_interpulse53_po_raw.fig');


fig_keep = f_black;
ax_keep = f_black.Children;
fig_move = f_duke;


line_move = findobj(fig_move,'type','line');
line_keep = findobj(fig_keep,'type','line');

for l = 1:numel(line_keep)
    if(strcmpi(line_keep(l).LineStyle,'-'))
        line_keep(l).LineStyle = '--';
    end
end

for l = 1:numel(line_move)
    new_handle = copyobj(line_move(l),findobj(fig_keep,'type','axes'));
    new_handle.XData = line_move(l).XData;
end


%%
fig_keep.Name = 'Han_duncan_artifactExamples_BlackrockDashed_DukeSolid_CathodicRed_AnodicBlue_50uA_chan4stimrec';

fpath = 'D:\Lab\Data\stim_ephys_paper\';
saveFiguresLIB(fig_keep,fpath,fig_keep.Name);

