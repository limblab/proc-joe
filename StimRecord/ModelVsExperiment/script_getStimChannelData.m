% get model data -- stim channel response
    mdl_input_data = [];
    mdl_input_data.folderpath = 'C:\Users\Joseph\Box\Miller-Grill_S1-stim\ModelData\StimChannelResponse\with_Intrinsic_Activity\';
    mdl_input_data.diam_list = [1,2,3];
    mdl_input_data.amp_list = [5,10,15,20,25,30,40,50,100];
    mdl_input_data.get_axon_dendrite_locs = 0;
    mdl_input_data.cell_id_list = [6,11,16,21];
    % cell_id=6:10 %L23 PC, clones 1-5
    % cell_id=11:15 %L4 LBC, clones 1-5
    % cell_id=16:20 %L5 PC, clones 1-5
    % cell_id=21:25 %L6 PC, clones 1-5

    [mdl_data_all,diam_all,amp_all,cell_id_all] = getModelStimChannelData(mdl_input_data);
    
% get experimental data -- stim channel response
    [exp_data_all, amp_all] = getExperimentStimChannelData(exp_input_data);
    
    
    