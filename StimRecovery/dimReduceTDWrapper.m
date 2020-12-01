function [ td_move_pca, model_info, idx_plot ] = dimReduceTDWrapper( td, model_name, array_name )


    pca_params = []; pca_params.signals = {([array_name,'_spikes'])}; pca_params.do_plot = 0;
    pca_params.num_dims = 14;
    pca_params.marg_names = {'time','target','time/target interaction'}; 
    
    if(strcmpi(model_name,'dpca')==1)
        [td_move_pca,model_info] = runDPCA(td,'target_direction', pca_params);
        idx_plot = find(model_info.which_marg == 2,2,'first');
        [~,signal_means] = centerSignals(td,struct('signals',[array_name,'_spikes']));
        model_info.mu = signal_means{1};
    elseif(strcmpi(model_name,'pca')==1)
        pca_params.center_data = true;
        [td_move_pca,model_info] = dimReduce(td,pca_params);
        model_info.W = model_info.w(:,1:model_info.num_dims);
        model_info.num_comps = model_info.num_dims;
        idx_plot = [2,3];
    else
        error('model_name not implemented');
    end


end

