%% Moving Horizon Parameter Estimation - repetitive Batch
% Autor: Jessica Behrens
% Date: 09.08.2024
% Publication:  Moving Horizon Parameter Estimation for an Enzyme Catalyzed
% Transamination Reaction with Integrated Product Removal - Jessica Behrens, Sven Tiedemann, Tom Kunde, Prof. Dr. Jan von Langermann, Prof. Dr.-Ing. Achim Kienle 


%% clean up everything
clear;
clear path
close all;
clc;

% set format
format shortEng

% supress warning from fmincon - 'Matrix is close to singular'
warning('off','MATLAB:nearlySingularMatrix')

% add filepath of 'General Functions'
addpath(genpath('...\General Functions'));

% add filepath of 'main_PL_in_MHE_input'
file_path_PL_in_MH = '...Profile_Likelihood\PL_in_MHE_rep_batch';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Model Setup 
Model_Setup

%% #############################################################################################################
%% create artificial measurement for MHE
Generate_Data

%% weighting matrixes - measurement noise   
V_diag = ones(p.n_y,1).*p.sigma_m.^2;
p.V = diag(V_diag);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Indexing and Bounds
Bounds

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% horizon length
Hor_Length =  10;

% number of disturbance scenarios
Num_dist = 5;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% bounds
% bounds of parameters
lb_para = ones(length(para),1).*1e-12;
ub_para= 5.*para';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MHE Simulation

X_hat_save = {};

% initial parameter pertubation
Para_dist = linspace(-0.5,1.5,Num_dist);
Para_dist = Para_dist.*ones(p.n_p,length(Para_dist));
for l = 1:p.n_p
  rand_permute =  randperm(length(Para_dist));
Para_dist(l,:) = Para_dist(l,rand_permute);
end

idx_feed_tH = [];%
x0_feed = [];%
idx_feed = [];

for Sim_idx = 1:Num_dist

 if Sim_idx ~= 1
   % create noisy data
    rand_num = randn(p.n_y,length(p.tspan));
    y_data = p.h(x_data) + proz_std_mess.*rand_num.*mean_measurement;
    y_data(y_data<0) = 0; % non-negative measurements
 end
 
  % disturb initial state  
    para_dist = Para_dist(:,Sim_idx)';
    para_start_MHE = para+ para_dist.*para;

    % initial state covariance matrix
    Initial_state_cov = (para_states_est-para_start_MHE).^2;
    P = eye(p.n_p);

    % matrix to save noisy estimates
    X_hat = zeros(p.n_p,length(p.tspan)+1);
    
    % initialize the states from the measurement
    X0 = zeros(p.n_ps,Hor_Length);
    X0(:,1) = para_start_MHE;     
    
    % prediction for intial state
    Para_pred = para_start_MHE';
    x0_tH = x_data(1,:);
    
    % save initial value
    X_hat(:,1) = [Para_pred]';

   P_save = {}; 
   P_save{1} = P; 
   Para_prev_save = {};
   Para_prev_save{1} = Para_pred;
   x0_tH_save = {};
   x0_tH_save{1} = x0_tH';
   idx_feed_save = [];

    for k = 1:length(p.tspan)
       
        % initial index of time horizon
        if k<=Hor_Length
            init_idx_est = 1; 
        else
            init_idx_est = k-Hor_Length+1;
             if ~isempty(idx_feed_tH)
                idx_feed_tH = idx_feed_tH-1;
             end
        end

        % estimation index of time horizon
        idx_est = [init_idx_est:k]; 

        % time horizon
        tH = [p.tspan(idx_est)];

        % data in time horizon
        y_tH = y_data(:,idx_est);
        
        % Matrix for initial state covariance
        P_tH = P;
    
        % bounds
        lb_tH = log10(lb_para+1e-24); % add small amout in case value = 0
        ub_tH = log10(ub_para);    
    
        % initial parameter value
        Para0 = log10(Para_pred+1e-24); % add small amout in case value = 0
        
        % initialize flag of solver
        flag = 0;
    
        while flag <= 0 % repeat calc when solver fails
        % call solver for MHE prediction
          [Para_pred,Obj,flag] = fmincon(@(decision_var)obj_MHE(decision_var,x0_tH,Para_pred,y_tH,tH,P_tH,x0_feed,idx_feed_tH,p),Para0,[],[],[],[],lb_tH,ub_tH,[],p.opt_fmincon);
           Para0 = Para_pred+randn(size(Para_pred)).*Para_pred; % tweak initial condition in case solver did not converge
        end
    
        % recalc logarithmic para
        Para_pred = 10.^Para_pred;
    
        % calculate initial condition by propagating intial estimate through ode system
        x_int_tH = int_ode_rep_batch(Para_pred,x0_tH,tH,idx_feed_tH,x0_feed,p)';
        
  
        % calculate Model preditcion
        tH_pred = [tH tH(end)+diff([p.tspan(1:2)])];
        x_int_pred = int_ode_mit_input(Para_pred,x0_tH,tH_pred,idx_feed_tH,x0_feed,p)';
    
        % save estimated values
        X_hat(:,k+1) = [Para_pred'];  
    
        if k < Hor_Length % start always from initial conditions
    
           Para_pred = para_start_MHE';
           x0_tH = x_int_tH(:,1); 
    
        else % start with UKF update in time horizon
        
            %% UKF Update
            % -------------------------------------------UKF Parameter------------------------------------------------------------------------- 
            % number of states and parameters to estimate
            L = p.n_p;

            % % number of sigma points
            n_sigma_p = 2*L+1;

            % tuning parameter 
            kappa = 3-L;
            eta = sqrt(L+kappa);

            % a priori estimate
            X_pred_UKF = Para_pred;

            % allocate space
            rk = zeros(1,2*L);
            di = zeros(L,1);

            % matrix square root
            sk = sqrtm(P);
            for i = 1:L% check for bound violation
                    sk_i = sk(:,i);  
                    di1 = min(abs(ub_para-X_pred_UKF)./abs(sk_i));
                    di2 = min(abs(lb_para-X_pred_UKF)./abs(sk_i));
                    di(i) = min([sqrt(L+kappa),di1,di2]);
            end
                 
            % calculate how weights for sigma points need to be adapted
            Sr = 2*sum(di);
            Denom_1 = 2*(L+kappa);
            Denom_2 = (Sr-(2*L+1)*sqrt(L+kappa));
            a = (2*kappa-1)./(Denom_1*Denom_2);
            b = 1/Denom_1 - (2*kappa-1)./(2*sqrt(L+kappa).*Denom_2);
         
            
            % weighting matrices
            Wi = [b; a*[di; di] + b];
 
            % preallocate space for matrix with sigma points
            Xi_pred_UKF = zeros(p.n_ps,n_sigma_p);
            gamma_pred_UKF = zeros(p.n_y,n_sigma_p);
            P_x = zeros(p.n_ps,p.n_ps);
            P_y = zeros(p.n_y,p.n_y);
            P_xy = zeros(p.n_ps,p.n_y);
     
 
            % -------------------------------------------UKF Update-------------------------------------------------------------------------
           
            % calculate sigma points
            sqrt_term = eta.*sqrtm(P);
            Xi_est = [X_pred_UKF, X_pred_UKF+sqrt_term, X_pred_UKF-sqrt_term];

            Xi = [X_pred_UKF, X_pred_UKF+di.*sqrtm(P), X_pred_UKF-di.*sqrtm(P)];

            % map sigma points through ...
            for q = 1:n_sigma_p
               [~,x_int_UKF] =  ode15s(@(t,c)ode_system_model(t,c,Xi(:,q),p),[0 p.Ts],x0_tH,p.opt); % ...state equation
               Xi_pred_UKF(:,q) = Xi(:,q);
               %Xi_pred_UKF(find(Xi_pred_UKF<0)) = 0;

               gamma_pred_UKF(:,q) = p.h(x_int_UKF(end,:)); % ...measurement equation 
    
            end
    
            % predicted states and measurements - weighted
            x_UKF = Xi_pred_UKF*Wi;
            y_UKF = gamma_pred_UKF*Wi;
    
            % update covariance matrices
            for q = 1:n_sigma_p
                P_x = P_x + Wi(q).*(Xi_pred_UKF(:,q)- x_UKF)*(Xi_pred_UKF(:,q)- x_UKF).';
                P_y = P_y + Wi(q).*(gamma_pred_UKF(:,q)- y_UKF)*(gamma_pred_UKF(:,q)- y_UKF).';
                P_xy = P_xy + Wi(q).*(Xi_pred_UKF(:,q)- x_UKF)*(gamma_pred_UKF(:,q)- y_UKF).';
            end
            % add measurement noise
            P_y = P_y + p.V;
    
            % Kalman gain
            K = P_xy/P_y;
    
            % update covariance matrix for next iteration
            P = P_x - K*P_y*K.';
    
            % care for numerical issues
            if any(eig(P)<0)
                P = nearestSPD(P); % replace by nearest postivite semi-definite matrix
            end
               
            % Update predicted state - Filter scheme
            if isempty(idx_feed_tH) || idx_feed_tH > 2
               x0_tH = x_int_tH(:,2);
            else
               x0_tH = x0_feed;
               idx_feed_tH = [];
            end

        end

        if idx_feed == 0
            idx_feed = [];   
        end
       %% -------------------------------------------------------------------------------------------------------------------------

        % check if input is necessary
        if x_int_pred(p.idx_3MAP_u,end) < 1 % readjust reactor, when molar amount of 3MAp in non-polar phase is below 1mmol

           idx_feed = idx_est(end)+1; % in last interval discovered
           idx_feed_tH = length(tH)+1; % in last interval discovered
           x0_feed = x_int_tH(:,end);
           x0_feed(p.idx_3MAP_u,1) = 15;   % add 150 mM 3MAP
           x0_feed(p.idx_IPA_3DPPA,1) = 15;% add 100 mM IPA-3DPPA
           x0_feed(p.idx_mu_p) = x0(p.idx_mu_p);
           x_int_tH(:,idx_feed_tH) = x0_feed;%

           % time left
           t_left = p.tspan(idx_feed:end); 

           %% generate new data 
           % solve ode
           [~,x_int] =  ode15s(@(t,c)ode_system_model(t,c,para,p),t_left,x0_feed,p.opt);

           % measurements
           y_data(:,idx_feed:end) = p.h(x_int) + proz_std_mess.*rand_num(:,idx_feed:end).*mean_measurement;
           y_data(y_data<0) = 0; % non-negative measurements 

           x_data(idx_feed:end,:) = x_int;
           idx_feed_save = [idx_feed_save idx_feed];
        end

       %% -------------------------------------------------------------------------------------------------------------------------
        % save integrated states
        X_int_save(:,idx_est) = x_int_tH(:,1:length(idx_est));

       % Plot results
        Plot_Moving_Horizon

        % save results
        Para_prev_save{Sim_idx,k+1} = Para_pred;
        P_save{Sim_idx,k+1} = P;
        x0_tH_save{Sim_idx,k+1} = x0_tH;
        
    end

    Y_data_save{Sim_idx} = y_data;
    X_hat_save{Sim_idx} = X_hat;
end


% Plot results
t = p.tspan./3600;

for k = 1:Num_dist
    for idx_plot = 1:p.n_p
    figure(7)
    subplot(p.n_p/2,2,idx_plot)
    hold on
    X_plot = X_hat_save{k}(idx_plot,:);
    plot(t(1:end-1),X_plot)
    yline(para(idx_plot),'k-')
    end
end

%% ######################################################################################################
%% save results
%% ######################################################################################################

%% ----------------------------------------------------- Plot PL in MHE --------------------------------------------------------------------
  
MHE_Para_prev_save = ' Para_prev_save';
MHE_P_save = 'P_save';
MHE_x0_tH_save = 'PL_x0_tH_save';
MHE_xint_tH_save = 'PL_xint_tH_save';
MHE_Y_data = 'Y_data_save';
time = 'time';
MHE_idx_feed_save = 'MHE_idx_feed_save';
MHE_X_hat = 'X_hat_save';

MHE_Para_prev_save_i = strcat([file_path,'\',MHE_Para_prev_save,'.mat']);
MHE_P_save_i = strcat([file_path,'\',MHE_P_save,'.mat']);
MHE_x0_tH_save_i = strcat([file_path,'\',MHE_x0_tH_save,'.mat']);
MHE_Y_data_i = strcat([file_path,'\',MHE_Y_data,'.mat']);
time_i = strcat([file_path,'\',time,'.mat']);
MHE_idx_feed_save_i = strcat([file_path,'\',MHE_idx_feed_save,'.mat']);

save(MHE_Para_prev_save_i,'Para_prev_save')
save(MHE_P_save_i,'P_save')
save(MHE_x0_tH_save_i,'x0_tH_save');
save(MHE_Y_data_i,'Y_data_save');
save(MHE_idx_feed_save_i,'idx_feed_save');
t = p.tspan./3600;%
save(time_i,'t');
































