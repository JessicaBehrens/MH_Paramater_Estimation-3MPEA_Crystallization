
% convert casadi to double
y_mess_double = y_data;
x_mess_double= x_data'; 


% Plot measured concentration and actual concentration
t_plot = p.tspan./3600; % time in h

fig1 = figure(1) ;
 for i = 1:p.n_s
       subplot(p.n_s/2,2,i) 
       hold on
       plot(t_plot,x_mess_double(i,:),'k-');
       plot(t_plot,y_mess_double(i,:),'bo');
       ylabel(Names_states(i))
       xlabel('t/h')     
 end
lh = legend('$xk$','$yk$','interpreter','latex','Location','east');

% plot molar amount of product salt + 3MPEA in liquid phase
MPEA_sim= x_data(:,2) + p.conv_mue3_2_n(x_data(:,11));
MPEA_measure = y_data(2,:) + p.conv_mue3_2_n(y_data(11,:));

% Plot Measurments
figure(2)
plot(t_plot,MPEA_sim,'k-.')
hold on
plot(t_plot,MPEA_measure,'bo');
ylabel(Names_states(2))
xlabel('t/h')
lh = legend('$xk$','$yk$','interpreter','latex','Location','east');

% Plot Parameter value
figure(3) 
for i = 1:p.n_p
  subplot(p.n_p/2,2,i) 
   hold on
   yline(para(i),'b');
   ylabel(Names_para(i))
   xlabel('Zeit [s]')
   ylim([0, para(i)*2])
end
