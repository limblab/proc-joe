% runs stim column model. This tries to answer the question 'how much
% activation should we expect in a single column?' This simulates a single
% column of homogenous neurons (radius = r_col), then drops electrodes at
% different steps (r_step), and computes the area overlap
% between the stimulation sphere (based on Stoney's results, I = current, K
% = excitability constant) and the column. The function returns the percent overlaps (perc_overlap)
% at each step. 

%% run model for a single set of parameters
    input_data.r_col = 0.2; % mm
    input_data.I = 30; % uA
    input_data.K = 1292; % uA/mm^2
    input_data.test_vals = 0:0.01:1; % probabilities
    input_data.num_steps = 100;
    
    data = stimColumnModel(input_data);
    
    
%% run model for different values of I, store prob more overlap for all
% values of I

    input_data.r_col = 0.2; % mm
    input_data.K = 1292; % uA/mm^2
    input_data.test_vals = 0:0.001:1; % probabilities
    input_data.num_steps = 1000;
    
    I_vals = [10:10:100];
    prob_more_overlap = zeros(numel(input_data.test_vals),numel(I_vals));
    colors = inferno(numel(I_vals)+3);
    
    figure();
    for i = 1:numel(I_vals)
        input_data.I = I_vals(i);
        
        data = stimColumnModel(input_data);
        
        prob_more_overlap(:,i) = data.prob_move_overlap;
        plot(input_data.test_vals,prob_more_overlap(:,i),'color',colors(i,:),'linewidth',1.5); hold on
    end
    
    l=legend(num2str(I_vals')); set(l,'box','off','location','best')
    formatForLee(gcf); xlabel('Frac overlap'); ylabel('P(col \cap stim > Frac overlap)')
    set(gca,'fontsize',14)
    ylim([0,1])
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    