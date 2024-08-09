%% Bounds of parameters

%% for Profile Likelihood
% index of parameters in parameter vector
p.idx_para_est = [1:p.n_p];
p.n_p_est = length(p.idx_para_est);

% parameters, that are estimated
para_states_est = [para(p.idx_para_est)];

% number of estimated parameters
p.n_ps = p.n_p_est;

% always first parameters, then initial conditions
para_states_full = [para x0]';
idx_full = [p.idx_para_est];

% bounds of parameters
lb_para = para - 0.1*para;
ub_para = para + 0.1*para;

% bounds of parameters for profile likelihood
lb_para_full_pl = lb_para;
ub_para_full_pl = ub_para;


