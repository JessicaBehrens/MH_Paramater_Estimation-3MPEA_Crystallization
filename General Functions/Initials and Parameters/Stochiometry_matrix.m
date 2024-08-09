function [p] = Stochiometry_matrix(p)
% input
%      p - paramter structure

% output
%      p.F - Stochiometric matrix


% liquid - polar phase
                % 1- IPA - Isopropylamine
                % 2- 3MPEA - (S)-1-(3-methoxyphenyl)ethylamine 
                % 3- Ac - Acetone
                % 4- 3MAP (aq) - 3-methoxyacetophenone
                % 5- 3DPPA - 3,3-diphenylpropionic acid 
                % 6- Biocatalyst(E)
% unpolar phase
                % 7- 3MAP (nonpolare phase) - 3-methoxyacetophenone   

% solid phase - product salt   
                % 8- mu0
                % 9- mu1
                % 10- mu2
                % 11- mu3

% solid phase - donor salt  
                % 12- IPA-3DPPA - Isopropylamine 3,3-diphenylpropionic acid

%%##################################################################################################################################################################################
%% Matrix for flow directions
%%##################################################################################################################################################################################
 
% chemical reaction - in polar phase
nue_p = [-1;1;1;-1;0;0];
nue_n = [0];
nue = [nue_p; nue_n]; 

% diffusion of 3MAP - from polar to unpolar phase
nue_3MAP_p = [0;0;0;1;0;0];
nue_3MAP_n = [-1];
nue_3MAP = [nue_3MAP_p; nue_3MAP_n]; 

% IPA Feed
nue_IPA_3DPPA_p = [1;0;0;0;1;0];
nue_IPA_3DPPA_n = [0];
nue_IPA_3DPPA = [nue_IPA_3DPPA_p; nue_IPA_3DPPA_n]; 

% evaporation
nue_evap_p = [0;0;-1;0;0;0];
nue_evap_n = [0];
nue_evap = [nue_evap_p; nue_evap_n ];  

% Product salt formation
nue_Cry_p_p = [0;-1;0;0;-1;0];
nue_Cry_p_n = [0];
nue_Cry = [nue_Cry_p_p; nue_Cry_p_n]; 

% Overall matrix
p.F = [nue,nue_Cry,nue_3MAP,nue_IPA_3DPPA,nue_evap]';

end
