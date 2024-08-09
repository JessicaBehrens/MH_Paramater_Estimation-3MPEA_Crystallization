function [dX] = ode_system_model(t,X,para,p)
        % dynamic model of a reactive crystallization process
        % input
        %      t - time
        %      X - states 
        %      para - model parameters 
        %      p - paramter structure
        
        % output
        %     dX - differential change of model states
        
        
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

        %% ########################################################### Parameter for optimization ###################################################################################################
        % parameters
        K_cat_f = para(1); % forward reaction rate constant [1/s]
        K_cat_r = K_cat_f/sqrt(p.K_eq/((p.K_m_MPEA*p.K_m_Ac)/(p.K_m_IPA*p.K_m_MAP))); % reverse reaction rate constant [1/s]
        k_3MAP = para(2);  % constant for miscibility rate [1/s]
        kd = para(3);      % constant for donor crystal dissolution [1/s]
        k_G = para(4);     % constant for prodcut crystal growth rate [mm/s]

        %% ####################################################################################################################################################################################

        % liquid phase concentrations
        n = X(1:p.n_comp_l,1); % [mM]
 
        % polar phase molar masses
        n1i = n(p.idx_IPA:p.idx_E);

        % moments
        mu_p = X(p.idx_mu_p);

        % liquid phase concentrations
        n_d = X(p.idx_mu_p(end)+1);

        %% ____________________________________________________ molar concentrations ________________________________________________________________________________________________________________________
        % polar phase
        Volume = p.Volume;

        % molar concentration  [mM] = [mmol/L]
        IPA = n1i(p.idx_IPA)/Volume;
        MPEA = n1i(p.idx_3MPEA )/Volume;
        Ac = n1i(p.idx_Ac)/Volume;
        MAP_aq = n1i(p.idx_3MAP_aq)/Volume; 
        DPPA = n1i(p.idx_3DPPA)/Volume; 
        E = n1i(p.idx_E)/Volume;

        % molar amounts  [mmol]
        n_MAP_aq = n1i(p.idx_3MAP_aq); 
        n_MAP_u = n(p.idx_3MAP_u); 

        %% ####################################################################################################################################################################################
        %                                                                      Reaction kinetics
        % ####################################################################################################################################################################################
        % forward and reverse velocity
        v_f = K_cat_f*E;
        v_r = K_cat_r*E;

        % nominator
        N = v_f*v_r*(IPA*MAP_aq-Ac*MPEA/p.K_eq);

        % denominator
        D = v_r*p.K_m_IPA*MAP_aq + v_r*p.K_m_MAP*IPA + v_f*p.K_m_MPEA/p.K_eq*Ac + v_f*p.K_m_Ac/p.K_eq*MPEA ....
               + v_r*IPA*MAP_aq + v_f*p.K_m_MPEA/(p.K_i_IPA*p.K_eq)*IPA*Ac + v_r*p.K_m_IPA/p.K_i_MPEA*MPEA*MAP_aq + v_f/p.K_eq*MPEA*Ac;

        % reaction rate
        gamma = N/D*Volume;

        %% ####################################################################################################################################################################################
        %                                                        3MAP Diffusion from 2nd liquid Phase
        % ####################################################################################################################################################################################

        S_3MAP = MAP_aq/p.Cs_3MAP;
        dec_3MAP = (S_3MAP-1);
      
        d_3MAP =  k_3MAP*(1-S_3MAP)*((p.sig(dec_3MAP))*n_MAP_aq + (1-p.sig(dec_3MAP))*n_MAP_u);
    
        %% ####################################################################################################################################################################################
        %                                                             Product Salt Crystallization
        % ####################################################################################################################################################################################

        [dmoment_p,r_dot_p] = Crystal_dynamics_moment(MPEA,DPPA,mu_p,k_G,p);

        %% ####################################################################################################################################################################################
        %                                                                    IPA & 3DPPA - control
        % ####################################################################################################################################################################################
        %----------------------------------------------------- IPA , 3DPPA ------------------------------------------------------------------------------------------------------
        Cs_IPA_3DPPA = p.Cs_IPA_3DPPA;

        % ################################################################################################################################
        % donor salt dissolution
        
        dec_min_IPA_3DPPA = (IPA-DPPA)/Cs_IPA_3DPPA;
        min_IPA_3DPPA = (1-p.sig(dec_min_IPA_3DPPA))*IPA +(p.sig(dec_min_IPA_3DPPA))*DPPA;
        S_IPA_3DPPA = min_IPA_3DPPA/Cs_IPA_3DPPA;
        dec_IPA_3DPPA= S_IPA_3DPPA-1;
        
        r_IPA_3DPPA  =  kd*(1-S_IPA_3DPPA)*(p.sig(dec_IPA_3DPPA )*min_IPA_3DPPA*Volume + (1-p.sig(dec_IPA_3DPPA))*n_d);
        dn_IPA_3DPPA_extdt = -r_IPA_3DPPA;

        %% ####################################################################################################################################################################################
        %                                                             gasous Phase
        % ####################################################################################################################################################################################       

        % evaporation
        N_dot_i = gamma; 

        %% ####################################################################################################################################################################################
        %                                                            liquid phase mol balance
        % ####################################################################################################################################################################################
        
        dM_scalar = [gamma,r_dot_p,d_3MAP, r_IPA_3DPPA,N_dot_i];
        dM = (dM_scalar*p.F)';
 
        %% ####################################################################################################################################################################################
        %                                                             Concatenate all dynamic variables
        % ####################################################################################################################################################################################
      
        dX = [dM;dmoment_p';dn_IPA_3DPPA_extdt];  
         
  end



  function [dmoment,r_dot] = Crystal_dynamics_moment(c1,c2,mu,k_G,p)

        % calc condition for growth
        dec_cry = (c2-c1)/p.c_s;
        min_cry = (1-p.sig(dec_cry))*c2+p.sig(dec_cry)*c1;
        S = min_cry/p.c_s;
        Act_growth = S-1;
   
        % growth rate - oversaturation (S>1) triggers growth
        G =  p.sig(Act_growth)*k_G*(S-1); % [g/(s*mm^2) / g/mm^3] = [mm/s]
 
        % calc moment dynamics
        dmoment(1) = 0;        % zeroth moment
        dmoment(2) = 1*G*mu(1); % first moment
        dmoment(3) = 2*G*mu(2); % second moment
        dmoment(4) = 3*G*mu(3); % third moment

        % flux from polar liquid to solid phase
        r_dot = p.density/p.molar_mass*dmoment(4)*p.shape_fac(4); % [mmol/s]
  end




