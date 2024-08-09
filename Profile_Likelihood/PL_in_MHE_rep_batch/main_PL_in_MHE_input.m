%% Profile Likelihood in moving time horizon with re-adjustment of substances and filtering
% Autor: Jessica Behrens
% Date: 09.08.2024
% Publication:  Moving Horizon Parameter Estimation for an Enzyme Catalyzed
% Transamination Reaction with Integrated Product Removal - Jessica Behrens, Sven Tiedemann, Tom Kunde, Prof. Dr. Jan von Langermann, Prof. Dr.-Ing. Achim Kienle 


clear all; 
close all; 
clc;
clear path

warning('off','MATLAB:singularMatrix')

% add filepath of 'General Functions'
addpath(genpath('...\General Functions'));

%% #############################################################################################################
%%  Model Setup 

Model_Setup

%% #############################################################################################################
%% get artificial data

Generate_Data

% standart deviation
std_mess = p.sigma_m;

%% weighting matrixes - measurement noise   
V_diag = ones(p.n_y,1).*std_mess.^2;
p.V = diag(V_diag);

%% #############################################################################################################
%% get Bounds

Bounds
% bounds of parameters
lb_para = [0.005 0 0 0];
ub_para = [0.025 2.5e-4 5e-4 2.5e-3];

% bounds of parameters for profile likelihood
lb_para_full_pl = lb_para;
ub_para_full_pl = ub_para;


% horizon length
Hor_Length =  10;

load('P_save.mat')
load('Para_prev_save.mat')
load("PL_x0_tH_save.mat")
load("Y_data_save.mat")
load("MHE_idx_feed_save.mat")
x0_tH_save = cell2mat(x0_tH_save);

%% #############################################################################################################
%% caluclate the parameter values to be evalutated
% percentage of varying of fixed parameter
bound_of_parameters = [lb_para_full_pl; ub_para_full_pl]';
perc_bounds_phi = bound_of_parameters./para_states_est';

% number of steps in log likelihood parameter sweap
N_LL = 30;

% find optimal solution
N_multistart = 16;

Para_opt_save = {};
Conf_Int_save = {};
Obj_save = {};
Obj_opt_save = {};
Obj_CI_save= {};
Var_para_save = {};
idx_feed_tH = [];
x0_feed = [];

for q = 1:length(p.tspan)

 % initial index of time horizon
        if q<=Hor_Length
            init_idx_est = 1; 
        else
            init_idx_est = q-Hor_Length+1;
             if ~isempty(idx_feed_tH)
                idx_feed_tH = idx_feed_tH-1;
             end     
        end

        if idx_feed_tH == 0
            idx_feed_tH = [];   
        end  

        % estimation index of time horizon
        idx_est = [init_idx_est:q];
    
        % time horizon
        tH = [p.tspan(idx_est)];
       
        % data in time horizon
        y_tH = Y_data_save{1}(:,idx_est);
        x_tH_0 = x0_tH_save(:,idx_est(1));

        % check if input is necessary
        if idx_est(end) == idx_feed_save
           idx_feed_tH = length(tH); % in last interval discovered
           x0_feed = x0_tH_save(:,idx_feed_save);
           x0_feed(p.idx_3MAP_u,1) = 15;   % add 150 mM 3MAP
           x0_feed(p.idx_IPA_3DPPA,1) = 15;% add 100 mM IPA-3DPPA
           x0_feed(p.idx_mu_p) = x0(p.idx_mu_p);
        end

%% #############################################################################################################
        
        Obj_opt = zeros(N_multistart,1);
        Para_opt = zeros(N_multistart,p.n_p);
        
        % previouse parameter estimate
        Para_prev = Para_prev_save{idx_est(1)};

        % weighting matrix arrival cost
        P = P_save{idx_est(1)};

        parfor i = 1:N_multistart
            decision_var = log10(para_states_full(1:4));
            decision_var = decision_var +rand(size(decision_var)).*decision_var;
        
            [para_opt,Obj_opt_mst,flag]  = fmincon(@(decision_var)obj_MHE_input(decision_var,x_tH_0,Para_prev,y_tH,tH,P,x0_feed,idx_feed_tH,p),decision_var,[],[],[],[],log10(lb_para),log10(ub_para),[],p.opt_fmincon);
        
            Obj_opt(i,1) = Obj_opt_mst;
            Para_opt(i,:) = 10.^para_opt;
        
        end
        idx_opt_mutistart = find(min(Obj_opt));
        
        Obj_opt = Obj_opt(idx_opt_mutistart);
        para_opt = Para_opt(idx_opt_mutistart,:)';
      
        % number of steps in log likelihood parameter sweap
        var_para = zeros(p.n_ps,N_LL);
        
        for m = 1:p.n_ps
            var_para(m,:) = linspace(bound_of_parameters(m,1),bound_of_parameters(m,2),N_LL); % Parameters to be evaluated
        end
       
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
                para_states_full_pl = Para_prev;%para;
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
                           [P_var_sol,obj,flag]  = fmincon(@(decision_var)obj_MHE_PL_input(decision_var,para_states_full_pl,x_tH_0,Para_prev,y_tH,tH,P,idx_free,idx_full,x0_feed,idx_feed_tH,p),decision_var,[],[],[],[],log10(lbw),log10(ubw),[],p.opt_fmincon);
        
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

        idx_opt = {};
        % save all optimized values in right order
        for idx_PL_para = 1:p.n_ps  % (parameter to be fixed)
        
            % idx of free parameter - decision variables
            idx_free = [1:p.n_ps];
            idx_free(find(idx_free ==idx_PL_para)) = [];
            Para_opt(idx_free,:,idx_PL_para) = Para_opt_PL(:,:,idx_PL_para); 
            Para_opt(idx_PL_para,:,idx_PL_para) = var_para(idx_PL_para,:);
        
            % find optimal values
            idx_opt{idx_PL_para} = find(Obj(idx_PL_para,:) == min(Obj(idx_PL_para,:)));
            para_opt(idx_PL_para,1) = Para_opt(idx_PL_para,idx_opt{idx_PL_para}(1),idx_PL_para);
            Obj_opt(idx_PL_para,1) = Obj(idx_PL_para,idx_opt{idx_PL_para}(1));
        end
        
        
        %% #############################################################################################################
        %% Confidenz interval
        
        % profile
        Log_Likelihood = Obj/2;
        
        % Value at maximum Likelihoof estimate
        Log_Likelihood_opt = Obj_opt./2;
        
        % L_hat function
        l_hat = Log_Likelihood - Log_Likelihood_opt;
        
        % confidence interval - 95%
        Alpha_quantil = chi2inv(0.95,p.n_p-1)/2;
        
        % preallocate space
        Obj_Ci = {};
        Para_Conf_Int = zeros(p.n_p,2);
        
        Obj_CI_ub = zeros(p.n_ps,2);
        Obj_CI_lb = zeros(p.n_ps,2);
        Bounds_CI_lb = zeros(p.n_ps,2);
        Bounds_CI_ub = zeros(p.n_ps,2);
 
        parfor m = 1:p.n_ps
            idx_phi = find(var_para(m,:)==para_states_est(m));
            idx_while_broke = find(while_broke(m,:)==1);
        
            % index of confidenz interval
            idx_CI = find(l_hat(m,:)<= Alpha_quantil);% in 95% confidenz interval?
        
            lb_dalpha_C_int = min(idx_CI); 
            ub_dalpha_C_int = max(idx_CI);
        
            % % Preliminary confidence interval based on PL points
            % Para_Conf_Int(m,1) = var_para(m,lb_dalpha_C_int); % lower bound
            % Para_Conf_Int(m,2) = var_para(m,ub_dalpha_C_int); % upper bound
            % 
            % 
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
        Obj_CI = Log_Likelihood_opt.*2 + Alpha_quantil.*2;
        
        % linear equation
        m_lb = diff(Obj_CI_lb,[],2)./diff(Bounds_CI_lb,[],2);
        n_lb = Obj_CI_lb(:,1)- m_lb.*Bounds_CI_lb(:,1);
        
        % Para value at lower bound of CI
        Conf_Int_lb = (Obj_CI-n_lb)./m_lb ;
        
        %% upper bound CI
        % linear equation
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
% save everything
Para_opt_save{q} = para_opt;
Conf_Int_save{q} = Conf_Int;
Obj_save{q} = Obj;
Obj_CI_save{q} = Obj_CI;
Obj_opt_save{q} = Obj_opt;
Var_para_save{q} = var_para;
end

%% ######################################################################################################
%% Plot of Pl in MHE Window
%% ######################################################################################################

Obj = {};
Obj_opt = {};
Var_para = {};
t = p.tspan./3600;


for q = 1:length(t)
    for m = 1:p.n_p     
        Obj{m}(q,:) = log10(Obj_save{q}(m,:));
        Obj_opt{m}(1,q) = log10(Obj_opt_save{q}(m));
        Var_para{m}(q,:) = Var_para_save{q}(m,:);
    end
end

fig1 = figure(1);
tiledlayout(ceil(p.n_p/2),2);

for m = 1:p.n_p     

   nexttile
   [xq,yq] = meshgrid(linspace(0,t(end),100), [linspace( Var_para{m}(2,1),Var_para{m}(2,end),100)]);

   
    Obj_m = Obj{m};
    Obj_opt_m = Obj_opt{m};

    l_hat = Obj_m - Obj_opt_m';
    [t_mesh,para_mesh] = meshgrid(t(2:end),Var_para{m});
    hold on
    surf(t.*ones(size(l_hat')),Var_para{m}',l_hat','EdgeColor','none','FaceColor','interp')%
    
     xlabel('$t/h$','interpreter','latex')
     ylabel(Names_para(m),'interpreter','latex')
     zlabel('log$_{10}$($\chi^2_{LSE}(\theta)$)','interpreter','latex')
     ylim([lb_para(m) ub_para(m)])

end

figure(1)
cb = colorbar;
cb.Layout.Tile = 'east';
ylabel(cb,'$\mathrm{log_{10}}(\mathrm{PL}_\Omega(\mathbf{\theta},\mathbf{t_\mathrm{H}}))$','interpreter','latex');

