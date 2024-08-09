function [Obj,x_int] = obj_MHE(decision_var,x0_tH,Para_prev,y_tH,tH,P,p)
%%##################################################################################################################################################################################
% objective of moving horizon estimator

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

    % integrate ode system
   if length(tH) == 1
         x_int = x0_tH;
   else
    [t,x_int] = ode15s(@(t,c)ode_system_model(t,c,Para,p),tH,x0_tH,p.opt);

     if length(tH) == 2 % solver calculates to many points in case of tH = 2
         x_int(2:end-1,:) = []; % reduce to point in tH
     end
   end

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
