function [Obj,x_int] = obj_MHE_PL_input(decision_var,para_states_full_pl,x0_tH,Para_prev,y_tH,tH,P,idx_free,idx_full,x0_feed,idx_feed_tH,p)
%%##################################################################################################################################################################################
% objective of moving horizon estimator to calculate profile likelihoods in
% time horizon with re-adjustment of substances and filtering

% input
%     decision_var           - decision variable (scaled logarithmically)
%     para_states_full_pl    - full vector of parameter to get right values
%     x0_tH                  - states in time horizon
%     Para_prev              - previouse parameter estimate
%     y_tH                   - measurment data in time horizon
%     tH                     - time horizon
%     P                      - weighting matrix for arrival cost
%     idx_free               - index of variable free for optimization
%     idx_full               - index of all varibales
%     x0_feed                - value of states when states are re-adjusted during feeding
%     idx_feed_tH            - index, when feeding occures in time horizon
%     p                      - parameter structure

% output
%     Obj                    - value of objective
%     x_int                  - integrated states starting from x0_tH in time horizon

%%##################################################################################################################################################################################
%% Calculate MHE Objective
%%##################################################################################################################################################################################
 
    % Logarithmic input
    P_s = Para_prev;        % parameter
    Decision_var = 10.^decision_var(:); % decision variable

    % create whole parameter vector in right order
    Para_states = para_states_full_pl(idx_full);
    Para_states(idx_free) = Decision_var;

    % parameters
    Para = P_s(1:p.n_p);
    Para(p.idx_para_est) = Para_states(1:p.n_p_est);
    % ------------------------------------------------------------- arrival cost ------------------------------------------------------------------------------------------------------

    % arrival cost - error between previouse estimate and current estimate
    Arr_error = Para-Para_prev; 

    % calculate arrival cost weighted by P
    Obj = Arr_error' * inv(P) * Arr_error; 
    
    % ------------------------------------------------------------- current states and measurements ------------------------------------------------------------------------------------------------------
      [x_int] = int_ode_mit_input(Para,x0_tH,tH,idx_feed_tH,x0_feed,p);
    % ------------------------------------------------------------- current states and measurements ------------------------------------------------------------------------------------------------------
   
    % model output
    y_model = p.h(x_int);

   % ------------------------------------------------------------- output error ------------------------------------------------------------------------------------------------------
     % initialize error
    if size(x_int,1) < length(tH) % check if ode solver failed
        Obj = 1e12;    % artificially set error to high number
    else % no error

        % output error
        e_y = y_tH - y_model;

        for j = 1:length(tH)
           Obj = Obj + e_y(:,j)' * inv(p.V) * e_y(:,j); % calculate obj
        end
    end
end
