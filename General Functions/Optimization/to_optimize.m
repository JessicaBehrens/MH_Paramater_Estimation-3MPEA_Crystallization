function [F] = to_optimize(dec_var,x0,y_data,p)
% input
%      dec_var - decision variable 
%      x0      - initial state
%      y_data  - measurement data
%      p       - paramter structure

% output
%      F - objective function value

%%##################################################################################################################################################################################
%% Obj for full parameter estimation
%%##################################################################################################################################################################################
 
    % Logarithmic input
    P_s = 10.^dec_var(:);

    % parameters
    Para = P_s(1:p.n_p);

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