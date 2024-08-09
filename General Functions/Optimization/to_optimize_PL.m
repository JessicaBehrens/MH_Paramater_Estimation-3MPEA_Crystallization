function [F] = to_optimize_PL(dec_var,full_para,x0,y_data,idx_free,idx_full,p)
% input
%      dec_var - decision variable 
%      full_para - full parameter vector
%      x0      - initial state
%      y_data  - measurement data
%      idx_free - index of decision variables
%      idx_full - index of all
%      p       - paramter structure

% output
%      F - objective function value

%%##################################################################################################################################################################################
%% Obj for profile likelihood parameter estimation
%%##################################################################################################################################################################################
 
    % Logarithmic input
    P_s = 10.^full_para(:);        % parameter
    Decision_var = 10.^dec_var(:); % decision variable

    % create whole parameter vector in right order
    Para_states = P_s(idx_full);
    Para_states(idx_free) = Decision_var;

    % parameters
    Para = P_s(1:p.n_p);
    Para(p.idx_para_est) = Para_states(1:p.n_p_est);

    % integrate ode system
    [~,x_int] = ode15s(@(t,c)ode_system_model(t,c,Para,p),p.tspan,x0,p.opt);

    % model output
    y_model = p.h(x_int);

    if size(x_int,1) < length(p.tspan) % check if ode solver failed
        F = 1e12;    % artificially set error to high number
    else
        % objective
        F = sum(sum((y_data - y_model).^2./(p.sigma_m).^2));
    end

 end