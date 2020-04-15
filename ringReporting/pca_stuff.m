%% subset td_all

    td_bump = trimTD(td_all,{'idx_bumpTime',0},{'idx_bumpTime',300});


% get PDs

    pd_params = [];
    pd_params.out_signals = 'LeftS1_spikes';
    pd_params.bootForTuning = 0;
    pd_params.num_boots = 1;
    pd_params.move_corr = 'vel';
    
    pd_all = getTDPDs(td_bump,pd_params);
    
    
% get PCA projection
    pca_params = [];
    pca_params.algorithm = 'pca';
    
    [td_bump,info] = dimReduce(td_bump,pca_params);

    
 %% plot pc space
    for tr = 1:numel(td_bump)
        
        proj = mean(td_bump(tr).LeftS1_pca(5:end,2:3));
        color_use = 0.75*(td_bump(tr).bumpDir + 180)/360 + [0,0,0];
        plot(proj(1),proj(2),'marker','.','color',color_use,'markersize',20)
        hold on
        
    end   
%% compare pd and pca_angle
    figure();
    pca_angle = atan2(info.w(:,3),info.w(:,2));

    plot(pd_all.velPD,pca_angle,'.','markersize',20)



    