%% do gpfa on td 

% Fields for params input struct:
%   .arrays      : which arrays to use, e.g. 'M1' (put in cell for multiple)
%   .save_dir    : directory to save GPFA results in Byron's code. Pass in [] to skip (default)
%   .method      : gpfa, pca, etc. See Byron's code
%   .xdim        : assumed number of latent dimensions (default to 8, find optimal using CV)
%   .kernsd      : kernal width in s (default to 0.03, find optimal using CV)
%   .bin_w       : bin size desired for GPFA in s (default to 0.03)

params.arrays = 'LeftS1';
params.save_dir = [];
params.method = 'gpfa';
params.xdim = 8;
params.kernsed = 0.03;
params.bin_w = 0.03;

gpfa_mdl = runGPFA(td_all,params);