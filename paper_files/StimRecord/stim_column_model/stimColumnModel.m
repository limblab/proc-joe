function [output_data] = stimColumnModel(input_data)

    % This tries to answer the question 'how much
    % activation should we expect in a single column?' This simulates a single
    % column of homogenous neurons (radius = r_col), then drops electrodes at
    % different steps (r_step), and computes the volume overlap
    % between the stimulation sphere (based on Stoney's results, I = current, K
    % = excitability constant) and the cortical column (2D model)

    
    
    % get percent overlap for each step away from the center, for each
    % slice of a sphere
    r_stim_max = sqrt((input_data.I/input_data.K));
    vol_steps = linspace(-r_stim_max,r_stim_max,input_data.num_steps);
    
    vol_overlap_data.area_overlap = zeros(input_data.num_steps,input_data.num_steps);
    vol_overlap_data.prob_r = zeros(input_data.num_steps,1);
    vol_overlap_data.r_list = zeros(1,numel(vol_steps));
   
    for z = 1:numel(vol_steps)
        % get radius of circular slice
        input_data.r_stim = sqrt(r_stim_max^2 - vol_steps(z)^2 - 0^2); % set y = 0
        
        area_overlap_data = getStimColumnAreaOverlap(input_data);
        vol_overlap_data.area_overlap(:,z) = area_overlap_data.area_overlap;
    end
    vol_overlap_data.prob_r(:,1) = area_overlap_data.prob_r;

    % sum data along slices, not along different radius positions
    vol_overlap_data.percent_vol_overlap = (vol_overlap_data.area_overlap*(mode(diff(vol_steps))*ones(size(vol_overlap_data.area_overlap,1),1)))/(4*pi/3*r_stim_max^3);
    
    % compute P(A > x) for test_vals
    prob_more_overlap = zeros(numel(input_data.test_vals),1);
    
    for i = 1:numel(input_data.test_vals)
        percent_overlap_mask = vol_overlap_data.percent_vol_overlap >= input_data.test_vals(i);
        prob_more_overlap(i) = sum(vol_overlap_data.percent_vol_overlap.*vol_overlap_data.prob_r.*percent_overlap_mask);
    end
    
    % package output_data
    output_data.percent_overlap = vol_overlap_data.percent_vol_overlap;
    output_data.r_stim = r_stim_max;
    output_data.r_steps = area_overlap_data.r_steps;
    output_data.vol_steps = vol_steps;
    output_data.prob_r = vol_overlap_data.prob_r;
    
    output_data.prob_move_overlap = prob_more_overlap;
    output_data.prob_test_vals = input_data.test_vals;
    
end



function [output_data] = getStimColumnAreaOverlap(input_data)

    % The function returns amount of overlap between the stim and cortical
    % column based on input data (I, K, num_steps, r_col). Also returns the
    % probability of an electrode being placed at a radius away. 


    r_stim = input_data.r_stim;
    r_steps = linspace(0,input_data.r_col,input_data.num_steps);
    dr = mode(diff(r_steps));
    
    area_overlap = zeros(numel(r_steps),1);
    prob_r = zeros(size(area_overlap));
    
    % assume column is at (0,0), using polar coordinates
    
    for r = 1:numel(r_steps)
        
        % check if stim_col is fully within column
        if(r_stim + r_steps(r) <= input_data.r_col)
            area_overlap(r) = pi*r_stim^2;
        elseif(input_data.r_col+r_steps(r) <= r_stim) % column is fully within stim
            area_overlap(r) = pi*input_data.r_col^2; % ratio of areas.
        else
            area_overlap(r) = getAreaOverlapTwoCircles(input_data.r_col,r_stim,r_steps(r));
        end
            
        % get P(r) for each r used
        if(r_steps(r) == 0)
            prob_r(r) = 0;
        else
            prob_r(r) = (r_steps(r)^2 - (r_steps(r)-dr)^2)/(input_data.r_col^2);
        end
        
    end
    
    
    output_data.r_steps = r_steps;
    output_data.area_overlap = area_overlap;
    output_data.prob_r = prob_r;
end


function [area_overlap] = getAreaOverlapTwoCircles(r_col,r_stim,d)

    % r1 and r2 are the radii of the two circles, d is the distance between
    % the two

    d_col = r_col - (d^2 + r_col^2 - r_stim^2)/(2*d);
    d_stim = (r_stim+d) - (d^2 + r_col^2 - r_stim^2)/(2*d);

    A1 = r_col^2*acos((r_col - (d_col))/r_col) - (r_col - (d_col))*sqrt(2*r_col*(d_col) - (d_col)^2);
    A2 = r_stim^2*acos((r_stim - (d_stim))/r_stim) - (r_stim - (d_stim))*sqrt(2*r_stim*(d_stim) - (d_stim)^2);

    area_overlap = pi*r_stim^2 - (A2-A1); % A2 and A1 are the areas outside of the intersection, so do total area - what we found.
end