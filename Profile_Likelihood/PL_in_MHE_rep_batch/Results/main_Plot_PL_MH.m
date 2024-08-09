%% Plot Results


clear all; 
close all; 
clc;

addpath(genpath('\\sthome.ovgu.de\jesbehre\Desktop\Promotion\FOR5538\Veröffentlichungen\CIT_2024\Matlab\Result_Plots'));
Names_settings

time = 'time';
PL_Theta = 'PL_Theta';
PL_Para_opt_save = 'PL_Para_opt_save';
PL_Conf_Int_save = 'PL_Conf_Int_save';
PL_Obj_save = 'PL_Obj_save';
PL_Obj_opt_save = 'PL_Obj_opt_save';
PL_Obj_CI_save = 'PL_Obj_CI_save';

load(time)
load(PL_Theta)
load(PL_Para_opt_save)
load(PL_Conf_Int_save)
load(PL_Obj_save)
load(PL_Obj_opt_save)
load(PL_Obj_CI_save)

% bounds of parameters
lb_para = [0.005 0 0 0];
ub_para = [0.025 2.5e-4 5e-4 2.5e-3];

p.n_p = 4;

% Names
Names_para = ["$K_\mathrm{f}$", '$k_\mathrm{3MAP}$', '$k_\mathrm{d}$', '$k_\mathrm{G}$'];

Obj = {};
Obj_opt = {};
Var_para = {};

for q = 1:length(t)-1
    for m = 1:p.n_p     
        Obj{m}(q,:) = log10(Obj_save{q}(m,:));
        Obj_opt{m}(1,q) = log10(Obj_opt_save{q}(m));
        Var_para{m}(q,:) = Var_para_save{q}(m,:);
    end
end


t = t./3600;
fig1 = figure(1);
tiledlayout(ceil(p.n_p/2),2);

for m = 1:p.n_p     

   nexttile
   [xq,yq] = meshgrid(linspace(0,t(end),100), [linspace( Var_para{m}(2,1),Var_para{m}(2,end),100)]);

    %subplot(ceil(p.n_p/2),2,m)
    Obj_m = Obj{m};
    Obj_opt_m = Obj_opt{m};

    l_hat = Obj_m - Obj_opt_m';
   
    [t_mesh,para_mesh] = meshgrid(t(2:end),Var_para{m});
   % vq = griddata(t_mesh,para_mesh,l_hat',xq,yq);
    hold on
    %surf(xq,yq,vq)
    surf(t(2:end).*ones(size(l_hat')),Var_para{m}',l_hat','EdgeColor','none','FaceColor','interp')%
    %surf(t_mesh,para_mesh,Obj);
    
     xlabel('$t/h$')
     ylabel(Names_para(m))
     zlabel('log$_{10}$($\chi^2_{LSE}(\theta)$)')
     ylim([lb_para(m) ub_para(m)])

end
figure(1)
    cb = colorbar;
    cb.Layout.Tile = 'east';
       ylabel(cb,'$\mathrm{log_{10}}(\mathrm{PL}_\Omega(\mathbf{\theta},\mathbf{t_\mathrm{H}}))$','Interpreter','latex');

saveas(fig1,strcat([file_path_save,'\','PL_MH_mit_P_fed_batch.png']))