function [ out ] = pooled_std( std_list, num_trial_list )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    out = sqrt(sum((std_list.^2).*(num_trial_list-1))/...
            sum(num_trial_list-1));
end

