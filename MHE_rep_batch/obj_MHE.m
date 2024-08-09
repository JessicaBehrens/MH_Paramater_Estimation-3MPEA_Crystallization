function [Obj,x_int] = obj_MHE(decision_var,x0_tH,Para_prev,y_tH,tH,P,x0_feed,idx_feed_tH,p)
%%##################################################################################################################################################################################
% input
%     decision_var - decision variable (scaled logarithmically)
%     x0_tH        - states in time horizon
%     Para_prev    - previouse parameter estimate
%     y_tH         - measurment data in time horizon
%     tH           - time horizon
%     P            - weighting matrix for arrival cost
%     p            - parameter structure

% output
%     Obj          - value of objective
%     x_int        - integrated states starting from x0_tH in time horizon

%%##################################################################################################################################################################################
%% Calculate MHE Objective
%%##################################################################################################################################################################################
 
    % Logarithmic input - get real parameter values
    Para = 10.^decision_var;

    % ------------------------------------------------------------- arrival cost ------------------------------------------------------------------------------------------------------

    % arrival cost - error between previouse estimate and current estimate
    Arr_error = Para-Para_prev; 

    % calculate arrival cost weighted by P
    Obj = Arr_error' * inv(P) * Arr_error; 

    % ------------------------------------------------------------- current states and measurements ------------------------------------------------------------------------------------------------------
    [x_int] = int_ode_mit_input(Para,x0_tH,tH,idx_feed_tH,x0_feed,p);

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
