function x = D_dct(alpha,N)

opt.w = 64;      % size of patches
opt.q = opt.w/2; % controls overlap
opt.threshold_factor = 1; % reduce influency if <1
opt.dct_type = 'orthogonal4';
opt.dct_type = 'redundant';
opt.remove_lowfreq = 1; % force to 0 the lowfrequency
opt.n = N; % assuming rows=cols

x = callback_localdct(alpha,-1,opt);

