set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');


% delete plot entries to generate moving horizon window
 if k>1
             delete(plt_1); delete(plt_2); delete(plt_3);delete(plt_4); delete(plt_5); delete(plt_6); delete(plt_7); delete(plt_8);
             delete(plt_9); delete(plt_10); delete(plt_11);delete(plt_12); delete(plt_13); delete(plt_14); delete(plt_15); delete(plt_16);
             delete(plt_17); delete(plt_18); delete(plt_19);delete(plt_20); delete(plt_21); delete(plt_22); delete(plt_23); delete(plt_24);
             delete(plt_25); delete(plt_26); delete(plt_27);delete(plt_28); delete(plt_29); delete(plt_30); delete(plt_31); delete(plt_32);
             delete(plt_33); delete(plt_34); delete(plt_35);delete(plt_36);delete(plt_37); delete(plt_40);delete(plt_41);delete(plt_42);delete(plt_43);
             delete(lh) ;delete(lh2);   

 end     
 
% measurement in time horizon- for Plot
t_max_in_h = p.tspan./3600;
Para_tH = X_hat(:,k);
t_H_plot = tH./3600;
x_int_tH_plt = x_int_tH(:,1:length(idx_est));
x_mess_double = x_data';

% Plot of states
   fig2 = figure(1);
       subplot(6,2,1)
       hold on 
           idx_state_plt = 1;
           plot(t_max_in_h(1:k),x_mess_double(idx_state_plt,1:k),'r-');
           plt_1 = plot(t_H_plot,x_int_tH_plt(idx_state_plt,:),'kx');hold on
           plt_2 = plot(t_max_in_h(1:k),y_data(idx_state_plt,1:k),'bo');hold on
           plt_3 = plot(t_max_in_h(1:k),X_int_save(idx_state_plt,1:k),'k-','Color','#642EFE'); 
           ylabel(Names_states(idx_state_plt))
           xlim([0 t_H_plot(end)+0.1])
       subplot(6,2,2)
       hold on 
           idx_state_plt = 2;
           plot(t_max_in_h(1:k),x_mess_double(idx_state_plt,1:k),'r-');
           plt_4 = plot(t_H_plot,x_int_tH_plt(idx_state_plt,:),'kx');hold on
           plt_5 = plot(t_max_in_h(1:k),y_data(idx_state_plt,1:k),'bo');hold on
           plt_6 = plot(t_max_in_h(1:k),X_int_save(idx_state_plt,1:k),'k-','Color','#642EFE'); 
           ylabel(Names_states(idx_state_plt))
           xlim([0 t_H_plot(end)+0.1])
       subplot(6,2,3)
       hold on 
           idx_state_plt = 3;
           plot(t_max_in_h(1:k),x_mess_double(idx_state_plt,1:k),'r-');
           plt_7 = plot(t_H_plot,x_int_tH_plt(idx_state_plt,:),'kx');hold on
           plt_8 = plot(t_max_in_h(1:k),y_data(idx_state_plt,1:k),'bo');hold on
           plt_9 = plot(t_max_in_h(1:k),X_int_save(idx_state_plt,1:k),'k-','Color','#642EFE'); 
           ylabel(Names_states(idx_state_plt))
           xlim([0 t_H_plot(end)+0.1])
       subplot(6,2,4)
       hold on 
           idx_state_plt = 4;
           plot(t_max_in_h(1:k),x_mess_double(idx_state_plt,1:k),'r-');
           plt_10 = plot(t_H_plot,x_int_tH_plt(idx_state_plt,:),'kx');hold on
           plt_11 = plot(t_max_in_h(1:k),y_data(idx_state_plt,1:k),'bo');hold on
           plt_12 = plot(t_max_in_h(1:k),X_int_save(idx_state_plt,1:k),'k-','Color','#642EFE'); 
           ylabel(Names_states(idx_state_plt))
           xlim([0 t_H_plot(end)+0.1])         
       subplot(6,2,5)
       hold on 
           idx_state_plt = 5;
           plot(t_max_in_h(1:k),x_mess_double(idx_state_plt,1:k),'r-');
           plt_13 = plot(t_H_plot,x_int_tH_plt(idx_state_plt,:),'kx');hold on   
           plt_14 = plot(t_max_in_h(1:k),y_data(idx_state_plt,1:k),'bo');hold on
           plt_15 = plot(t_max_in_h(1:k),X_int_save(idx_state_plt,1:k),'k-','Color','#642EFE'); 
           ylabel(Names_states(idx_state_plt))
           xlim([0 t_H_plot(end)+0.1])
        subplot(6,2,6) 
        hold on 
           idx_state_plt = 6;
           plot(t_max_in_h(1:k),x_mess_double(idx_state_plt,1:k),'r-');
           plt_16 = plot(t_H_plot,x_int_tH_plt(idx_state_plt,:),'kx');hold on 
           plt_17 = plot(t_max_in_h(1:k),y_data(idx_state_plt,1:k),'bo');hold on
           plt_18 = plot(t_max_in_h(1:k),X_int_save(idx_state_plt,1:k),'k-','Color','#642EFE'); 
           ylabel(Names_states(idx_state_plt))
           xlim([0 t_H_plot(end)+0.1])
        subplot(6,2,7) 
        hold on 
           idx_state_plt = 7;
           plot(t_max_in_h(1:k),x_mess_double(idx_state_plt,1:k),'r-');
           plt_19 = plot(t_H_plot,x_int_tH_plt(idx_state_plt,:),'kx');hold on  
           plt_20 = plot(t_max_in_h(1:k),y_data(idx_state_plt,1:k),'bo');hold on
           plt_21 = plot(t_max_in_h(1:k),X_int_save(idx_state_plt,1:k),'k-','Color','#642EFE'); 
           ylabel(Names_states(idx_state_plt))
           xlim([0 t_H_plot(end)+0.1])
        subplot(6,2,8) 
        hold on 
           idx_state_plt = 8;
           plot(t_max_in_h(1:k),x_mess_double(idx_state_plt,1:k),'r-');
           plt_22 = plot(t_H_plot,x_int_tH_plt(idx_state_plt,:),'kx');hold on 
           plt_23 = plot(t_max_in_h(1:k),y_data(idx_state_plt,1:k),'bo');hold on
           plt_24 = plot(t_max_in_h(1:k),X_int_save(idx_state_plt,1:k),'k-','Color','#642EFE'); 
           ylabel(Names_states(idx_state_plt))
           xlim([0 t_H_plot(end)+0.1])
       subplot(6,2,9) 
       hold on 
           idx_state_plt = 9;
           plot(t_max_in_h(1:k),x_mess_double(idx_state_plt,1:k),'r-');
           plt_25 = plot(t_H_plot,x_int_tH_plt(idx_state_plt,:),'kx');hold on 
           plt_26 = plot(t_max_in_h(1:k),y_data(idx_state_plt,1:k),'bo');hold on
           plt_27 = plot(t_max_in_h(1:k),X_int_save(idx_state_plt,1:k),'k-','Color','#642EFE'); 
           ylabel(Names_states(idx_state_plt))
           xlim([0 t_H_plot(end)+0.1])
       subplot(6,2,10)
       hold on 
           idx_state_plt = 10;
           plot(t_max_in_h(1:k),x_mess_double(idx_state_plt,1:k),'r-');
           plt_28 = plot(t_H_plot,x_int_tH_plt(idx_state_plt,:),'kx');hold on 
           plt_29 = plot(t_max_in_h(1:k),y_data(idx_state_plt,1:k),'bo');hold on
           plt_30 = plot(t_max_in_h(1:k),X_int_save(idx_state_plt,1:k),'k-','Color','#642EFE'); 
           ylabel(Names_states(idx_state_plt))
           xlim([0 t_H_plot(end)+0.1])
       subplot(6,2,11) 
       hold on 
           idx_state_plt = 11;
           plot(t_max_in_h(1:k),x_mess_double(idx_state_plt,1:k),'r-');
           plt_31 = plot(t_H_plot,x_int_tH_plt(idx_state_plt,:),'kx');hold on 
           plt_32 = plot(t_max_in_h(1:k),y_data(idx_state_plt,1:k),'bo');hold on
           plt_33 = plot(t_max_in_h(1:k),X_int_save(idx_state_plt,1:k),'k-','Color','#642EFE'); 
           ylabel(Names_states(idx_state_plt))
           xlim([0 t_H_plot(end)+0.1])
           xlabel('$t/h$')
       subplot(6,2,12) 
       hold on 
           idx_state_plt = 12;
           plt_37 = plot(t_max_in_h(1:k),x_mess_double(idx_state_plt,1:k),'r-');
           plt_34 = plot(t_H_plot,x_int_tH_plt(idx_state_plt,:),'kx');hold on  
           plt_35 = plot(t_max_in_h(1:k),y_data(idx_state_plt,1:k),'bo');hold on
           plt_36 = plot(t_max_in_h(1:k),X_int_save(idx_state_plt,1:k),'k-','Color','#642EFE','HandleVisibility','off'); 
           ylabel(Names_states(idx_state_plt))
           xlim([0 t_H_plot(end)+0.1])
           xlabel('$t/h$')
           lh = legend('$xk$','$\hat{x}k$','$yk$','interpreter','latex','Location','east'); 
 
         
% Plot of parameter estimates 
       figure(2) 
        hold on
        subplot(p.n_p/2,2,1)
        hold on
        idx_plot_para = 1;
        plt_40 = plot(t_max_in_h(1:k),X_hat(idx_plot_para,2:k+1),'k-','Color','#642EFE','Marker','.'); 
        yline(para(idx_plot_para),'k-');
        ylabel(Names_para(idx_plot_para))
        xlim([0 t_H_plot(end)+0.1])
        ylim([0.005  0.025])

        subplot(p.n_p/2,2,2)
        hold on
        idx_plot_para = 2;
        plt_41 = plot(t_max_in_h(1:k),X_hat(idx_plot_para,2:k+1),'k-','Color','#642EFE','Marker','.');
        yline(para(idx_plot_para),'k-');
        ylabel(Names_para(idx_plot_para))
        xlim([0 t_H_plot(end)+0.1])
        ylim([0  2.5e-4])


        subplot(p.n_p/2,2,3)
        hold on
        idx_plot_para = 3;
        plt_42 = plot(t_max_in_h(1:k),X_hat(idx_plot_para,2:k+1),'k-','Color','#642EFE','Marker','.');
        xlim([0 t_H_plot(end)+0.1])
        yline(para(idx_plot_para),'k-');
        ylabel(Names_para(idx_plot_para))
        ylim([0  4.5e-4])
        xlabel('$t/h$')

        subplot(p.n_p/2,2,4)
        hold on
        idx_plot_para = 4;
        plt_43 = plot(t_max_in_h(1:k),X_hat(idx_plot_para,2:k+1),'k-','Color','#642EFE','Marker','.');
        xlim([0 t_H_plot(end)+0.1])
        yline(para(idx_plot_para),'k-');
        ylabel(Names_para(idx_plot_para))
        xlabel('$t/h$')
        ylim([0 2.5e-3])
        lh2 = legend('$\hat{\theta}k$','$\theta$','interpreter','latex','Location','east'); 


pause(0.1)
