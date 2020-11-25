function [sig_out] = filtfilt_pad( b,a,sig,num_pad )
% only works on 1d arrays, must be horizontal



sig = [sig(1)*ones(1,num_pad),sig,sig(end)*ones(1,num_pad)];

sig_filt = filtfilt(b,a,sig);

sig_out = sig_filt(num_pad+1:numel(sig_filt)-num_pad);


end

