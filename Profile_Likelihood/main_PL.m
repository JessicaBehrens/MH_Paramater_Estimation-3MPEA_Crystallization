%% Profile Likelihood
% Autor: Jessica Behrens
% Date: 08.08.2024
% Publication:  Moving Horizon Parameter Estimation for an Enzyme Catalyzed
% Transamination Reaction with Integrated Product Removal - Jessica Behrens, Sven Tiedemann, Tom Kunde, Prof. Dr. Jan von Langermann, Prof. Dr.-Ing. Achim Kienle 

clear all; 
close all; 
clc;
clear path

% add filepath of 'General Functions'
addpath(genpath('...\General Functions'));

%% #############################################################################################################
%%  Model Setup 

Model_Setup

%% #############################################################################################################
%% get artificial data

Generate_Data

%% #############################################################################################################
%% get Bounds
Bounds


%% #############################################################################################################
% PLot measurements
Plot_Measurement

%% #############################################################################################################
%% caluclate the parameter values to be evalutated
% percentage of varying of fixed parameter
bound_of_parameters = [lb_para_full_pl; ub_para_full_pl]';
perc_bounds_phi = bound_of_parameters./para_states_est';

% number of steps in log likelihood parameter sweap
N_LL = 100;
dpara = abs(bound_of_parameters(:,1)-bound_of_parameters(:,2))/N_LL; % step size per parameter
var_para = zeros(p.n_ps,N_LL+1);

for l = 1:N_LL
    for m = 1:p.n_ps
    var_para(m,l) = bound_of_parameters(m,1) + (l-1)*dpara(m); % Parameters to be evaluated
    end
end

% find optimal solution
N_multistart = 16;
Obj_opt_save = zeros(N_multistart,1);
Para_opt_save = zeros(N_multistart,p.n_p);

parfor i = 1:N_multistart
    decision_var = log10(para);
    decision_var = decision_var +rand(size(decision_var)).*decision_var;

    [para_opt,Obj_opt_mst,flag]  = fmincon(@(decision_var)to_optimize(decision_var,x0,y_data,p),decision_var,[],[],[],[],log10(lb_para),log10(ub_para),[],p.opt_fmincon);

    Obj_opt_save(i,1) = Obj_opt_mst;
    Para_opt_save(i,:) = 10.^para_opt;

end
idx_opt_mutistart = find(min(Obj_opt_save));

Obj_opt = Obj_opt_save(idx_opt_mutistart);
para_opt = Para_opt_save(idx_opt_mutistart,:)';

% save optimal value to varying parameter vector
var_para(:,N_LL+1) = para_opt;

% sort parameters in ascending order
var_para = sort(var_para,2);

%% #############################################################################################################
%% Calculate profile Likelihoods
% preallocate matrices
Obj = zeros(p.n_ps,size(var_para,2));

Para_opt_PL = zeros(p.n_ps-1,N_LL+1,p.n_ps);
Para_opt = zeros(p.n_ps,N_LL+1,p.n_ps);

% marke when while iteration is left
while_broke = zeros(p.n_ps,size(var_para,2));

% index of profile likelihood parameter
for idx_PL_para = 1:p.n_ps  % (parameter to be fixed)

    % idx of free parameter - decision variables
    idx_free = [1:p.n_ps];
    idx_free(find(idx_free ==idx_PL_para)) = [];

    % Parameter bounds
    lbw = lb_para(idx_free); % lower bounds
    ubw = ub_para(idx_free); % upper bounds

    %% Solve Profile Likelihood
    parfor k = 1:size(var_para,2) % vary one parameter at a time

        % fixed parameter for profile likelihood
        para_states_full_pl = para;
        para_states_full_pl(idx_full(idx_PL_para)) = var_para(idx_PL_para,k);

        % set solver sattus to zero for while loop
        solver_success = false;

        % varying parameter
        para_var = para(idx_free);    

        % initialize indexing
        idx_while = 0;

        % set 
        P_var_sol = [];

               while solver_success ~= 1  

                    % varying phi
                     decision_var = log10(para_var);
                     decision_var(find(decision_var==0)) = 1e-16;   

                    % call fmincon
                   [P_var_sol,obj,flag]  = fmincon(@(decision_var)to_optimize_PL(decision_var,log10(para_states_full_pl),x0,y_data,idx_free,idx_full,p),decision_var,[],[],[],[],log10(lbw),log10(ubw),[],p.opt_fmincon);

                   % in case of solver error - tweak initial condion a bit
                   para_var = decision_var + (2*(rand(size(decision_var)))-1).*decision_var;
                   para_var(para_var < 0) = 0;  
                   para_var(isnan(para_var)) = 0;

                    if flag > 0
                       solver_success = 1; 
                    end

                    % count number of while iterations
                    idx_while = idx_while + 1;

                    if idx_while >10 % stop if more than 10 iterations and save idx
                        while_broke(idx_PL_para,k) = 1;
                        break
                    end
               end

        % save value of objective function
        Obj(idx_PL_para,k) = obj;

        % save value of optimized parameters
        Para_opt_PL(:,k,idx_PL_para) = 10.^P_var_sol; 

    end
idx_PL_para
end

% save all optimized values in right order
for idx_PL_para = 1:p.n_ps  % (parameter to be fixed)

    % idx of free parameter - decision variables
    idx_free = [1:p.n_ps];
    idx_free(find(idx_free ==idx_PL_para)) = [];
    Para_opt(idx_free,:,idx_PL_para) = Para_opt_PL(:,:,idx_PL_para); 
    Para_opt(idx_PL_para,:,idx_PL_para) = var_para(idx_PL_para,:);

    % find optimal values
    idx_opt(idx_PL_para,1) = find(Obj(idx_PL_para,:) == min(Obj(idx_PL_para,:)));
    para_opt(idx_PL_para,1) = Para_opt(idx_PL_para,idx_opt(idx_PL_para,1),idx_PL_para);
    Obj_opt(idx_PL_para,1) = Obj(idx_PL_para,idx_opt(idx_PL_para,1));
end


%% #############################################################################################################
%% Confidenz interval

% profile
Log_Likelihood = Obj./2;

% Value at maximum Likelihoof estimate
Log_Likelihood_opt = Obj_opt./2;
p.Log_Likelihood_opt_m = Log_Likelihood_opt;

% L_hat function
l_hat = Log_Likelihood - Log_Likelihood_opt;

% confidence interval - 95%
Alpha_quantil = chi2inv(0.95,1)/2;

% preallocate space
Obj_Ci = {};
Conf_Int_save = {};
Para_Conf_Int = zeros(p.n_p,2);

for m = 1:p.n_ps
    idx_phi = find(var_para(m,:)==para_states_est(m));
    idx_while_broke = find(while_broke(m,:)==1);

    % index of confidenz interval
    idx_CI = find(l_hat(m,:)<= Alpha_quantil);% in 95% confidenz interval?

    lb_dalpha_C_int = min(idx_CI); 
    ub_dalpha_C_int = max(idx_CI);

    % Preliminary confidence interval based on PL points
    Para_Conf_Int(m,1) = var_para(m,lb_dalpha_C_int); % lower bound
    Para_Conf_Int(m,2) = var_para(m,ub_dalpha_C_int); % upper bound

  
    % Interval for linear interpolation
    idx_lb_dalpha_C_int = lb_dalpha_C_int-1;
    idx_ub_dalpha_C_int = ub_dalpha_C_int+1; 

    % correct bounds if bounds of parameter are met    
    if lb_dalpha_C_int-1 < 1
        idx_lb_dalpha_C_int = 1;
    end
    if ub_dalpha_C_int+1 > size(var_para,2)
        idx_ub_dalpha_C_int = size(var_para,2);
    end

    Bounds_CI_lb(m,:) = [var_para(m,idx_lb_dalpha_C_int) var_para(m,lb_dalpha_C_int),];
    Bounds_CI_ub(m,:) = [var_para(m,ub_dalpha_C_int) var_para(m,idx_ub_dalpha_C_int)];

    Obj_CI_lb(m,:) = [Obj(m,idx_lb_dalpha_C_int) Obj(m,lb_dalpha_C_int),];
    Obj_CI_ub(m,:) = [Obj(m,ub_dalpha_C_int) Obj(m,idx_ub_dalpha_C_int)];

end

%% lower bound CI
% value of Obj at CI
Obj_CI = Log_Likelihood_opt.*2 + Alpha_quantil*2;

% Lineare Geradengl.
m_lb = diff(Obj_CI_lb,[],2)./diff(Bounds_CI_lb,[],2);
n_lb = Obj_CI_lb(:,1)- m_lb.*Bounds_CI_lb(:,1);

% Para value at lower bound of CI
Conf_Int_lb = (Obj_CI-n_lb)./m_lb ;

%% upper bound CI
% Lineare Geradengl.
m_ub = diff(Obj_CI_ub,[],2)./diff(Bounds_CI_ub,[],2);
n_ub = Obj_CI_ub(:,1)- m_ub.*Bounds_CI_ub(:,1) ;

% Para value at lower bound of CI
Conf_Int_ub = (Obj_CI-n_ub)./m_ub ;

% CI is on bounds
for m = 1:p.n_ps  
    if Bounds_CI_lb(m,1) == Bounds_CI_lb(m,2)
        Conf_Int_lb(m) = Bounds_CI_lb(m,1);
    end 
    if Bounds_CI_ub(m,1) == Bounds_CI_ub(m,2)
        Conf_Int_ub(m) = Bounds_CI_ub(m,1);
    end     
end

% full CI
Conf_Int = [Conf_Int_lb  Conf_Int_ub];

%% ######################################################################################################
%% Plot of Confidence Intervals
%% ######################################################################################################

Plot_CI

%% ######################################################################################################
%% Test confidenz interval
%% ######################################################################################################
% calculate the test cases yourself or load("5000_para_est.mat")

N_Test_samples = 5000; % number of test samples

% %% ---------------------------------------------------------------------------------------------------------
% Save_Para_opt = zeros(p.n_ps,N_Test_samples);
% Obj_test = zeros(N_Test_samples,1);
% 
% % bounds for optimization
% lb_para = ones(length(para),1).*1e-12;
% ub_para= 5.*para';
% 
% % generate data sets to test CI
% parfor i = 1:N_Test_samples
%     % new data
%     y_data = p.h(x_data) + proz_std_mess.*randn(p.n_y,length(p.tspan)).*mean_measurement;
% 
%     P_var_sol = [];
% 
%     % varying parameter
%     para_var = para;    
% 
%     idx_while = 0;
%     obj = 1e12;
%     solver_success = false;
% 
%      for w = 1:1
%               while solver_success ~= 1  
% 
%                   %  varying phi
%                     decision_var2 = log10(para_var);
%                     decision_var2(find(decision_var2==0)) = 1e-16;   
% 
%                     % call fmincon
%                    [P_var_sol,obj,flag]  = fmincon(@(decision_var2)to_optimize(decision_var2,x0,y_data,p),decision_var2,[],[],[],[],log10(lb_para),log10(ub_para),[],p.opt_fmincon);
% 
%                    % in case of solver error
%                    para_var = decision_var2 + (2*(rand(size(decision_var2)))-1).*decision_var2;
%                    para_var(find(para_var < 0)) = 0;
% 
%                     if flag > 0
%                        solver_success = 1; 
%                     end
%                     idx_while = idx_while + 1;
% 
%                     if idx_while >1
%                         while_broke(i) = 1;
%                         break
%                     end
%               end
%               para_var = 10.^(P_var_sol + 0.1.*P_var_sol.*randn(size(P_var_sol)));
%      end
% 
%         % save value of objective function
%         para_var = 10.^P_var_sol + (2.*rand(size(P_var_sol))-1).*10.^P_var_sol;
% 
%         % save value of objective function
%         Obj_test(i) = obj;
% 
%         % save value of optimized parameters
%         Save_Para_opt(:,i) = 10.^P_var_sol; 
% 
% end

%% ---------------------------------------------------------------------------------------------------------
load("5000_para_est.mat")
%% ---------------------------------------------------------------------------------------------------------

% count number of parameter estimates in CI
count_in_conf = [0;0;0];
    for j = 1:p.n_ps
        % index of optimal parameter
        count_in_conf(j,1) = sum(and(Save_Para_opt(j,:)<= Conf_Int(j,2),Save_Para_opt(j,:)>= Conf_Int(j,1))); 
    end

% percentage in CI
in_CI = count_in_conf./N_Test_samples.*100
