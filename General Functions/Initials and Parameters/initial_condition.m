function [initial_cond] = initial_condition(p)
% input
%      p - paramter structure

% output
%      initial_cond - initial conditions


%%##################################################################################################################################################################################
%% Initial conditions
%%##################################################################################################################################################################################

%% molar amounts and concentrations - solid and liquid phase
%_____________________________________________________________ IPA ________________________________________________________________________________________________________________________


% IPA - isopropylamine
c0_IPA =  100; % mM

%_____________________________________________________________ 3MPEA ________________________________________________________________________________________________________________________

% 3MPEA - (S)-1-(3-methoxyphenyl)ethylamine 
c0_3MPEA = 0; % mM

%_____________________________________________________________ Ac ________________________________________________________________________________________________________________________

% Ac - Acetone
c0_Ac = 0; % mM

%_____________________________________________________________ 3DPPA ________________________________________________________________________________________________________________________

% 3DPPA - 3,3-diphenylpropionic acid 
c0_3DPPA = 55; % mM 

%_____________________________________________________________ Enzyme ________________________________________________________________________________________________________________________

% calculate Enzym mass 
U_ml = 80; % activity in solution [U/ml]
A_E = 170;  % acivity of enzym [U/mg]
c0_E = (U_ml/A_E)/p.M_E;  % mM

%_____________________________________________________________ 3MAP ________________________________________________________________________________________________________________________

% 3MAP - 3-methoxyacetophenone - total amount
c0_3MAP_total =  150; % mM 

% 3MAP - unpolar, calculate part of aquaouse phase
c0_3MAP_aq = 0; % aquaouse phase
c0_3MAP_u = c0_3MAP_total; % unpolar phase 

% molar mass
n0_3MAP_u = c0_3MAP_u.*p.Volume; 

%____________________________________________________________ polar phase ________________________________________________________________________________________________________________________

% polar phase - molar amounts
initial_cond_p = [c0_IPA,c0_3MPEA,c0_Ac,c0_3MAP_aq,c0_3DPPA,c0_E].*p.Volume;  % mmol

%_______________________________________________________ Product salt 3MPEA-3DPPA________________________________________________________________________________________________________________________

%% initial conditions
% 3MPEA-3DPPA -(S)-1-(3-Methoxyphenyl)ethylammonium 3,3-diphenylpropionate (product salt)

% number of initial particles
N0 = 1e2;

% mean length of initial particles
L_0 = 0.01;

% zeroth moment- number of particles
mu_0 = N0; 

% first moment = total length [mm]
mu_1 = L_0*N0;

% second moment = total surface [mm^2]
mu_2 = L_0^2*N0/p.shape_fac(3);

% third moment = total volume [mm^3]
mu_3 = L_0^3*N0/p.shape_fac(4); 

moments_p = [mu_0 mu_1 mu_2 mu_3];

%_______________________________________________________ Donor salt IPA-3DPPA________________________________________________________________________________________________________________________

% IPA-3DPPA - isopropylammonium 3,3-diphenylpropionate (donor salt)
n_IPA_3DPPA_ext = 150*p.Volume; % mmol

%____________________________________________________ concatenate initial conditions ________________________________________________________________________________________________________________________

% initial conditions
initial_cond = [initial_cond_p,n0_3MAP_u,moments_p,n_IPA_3DPPA_ext]; 

end