% ---------------------------------------------------------------------------- Model Setup ----------------------------------------------------------------------
% initialise random number generator
rng(0298,'twister')

%% General Parameter

% call parameter
p = parameter;

% call Stochiometry matrix for reaction
[p] = Stochiometry_matrix(p);

k_G = 1e-3;      % constant for prodcut crystal growth rate mm/s
kd = 1e-4;       % constant for donor crystal dissolution 1/s
k_3MAP = 1e-4;   % constant for miscibility rate 1/s
K_cat_f = 1.0e-2;% forward reaction rate constant 1/s
p.Volume = 0.1;  % reactor volume L

% parameter vector
para = [K_cat_f k_3MAP kd k_G];
p.n_p = length(para); % number of parameters

% Names (for plots)
Names_states = ["$n^{\prime\bullet}_\mathrm{IPA}$", '$n^{\prime\bullet}_\mathrm{3MPEA}$', '$n^{\prime\bullet}_\mathrm{Ac}$', '$n^{\prime\bullet}_\mathrm{3MAP}$', '$n^{\prime\bullet}_\mathrm{3DPPA}$'...
              '$n^{\prime\bullet}_\mathrm{E}$', '$n^{\prime\circ}_\mathrm{3MAP}$', '$\mu_{0}$', '$\mu_{1}$', '$\mu_{2}$', '$\mu_{3}$', '$n^{\prime\prime}_\mathrm{IPA-3DPPA}$','$n^{\prime\prime}_\mathrm{3MPEA-3DPPA}$'];
Names_para = ["$K_\mathrm{f}$", '$k_\mathrm{3MAP}$', '$k_\mathrm{d}$', '$k_\mathrm{G}$'];

Names_full = [Names_para, Names_states];