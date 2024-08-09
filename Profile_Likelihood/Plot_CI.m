%% Plot of CI


for m = 1:p.n_p

     parameter_val = var_para(m,:);

     figure(4);
     subplot(ceil(p.n_p/2),2,m)
     plot(parameter_val,log10(Obj(m,:)),'b-')
     hold on 
     plot(para_opt(m),log10(Obj_opt(m)),'r*')
     %plot(Conf_Int(:,1),Conf_Int(:,2),'k-')
     y2 = xline(Conf_Int(m,1),'k-',num2str(Conf_Int(m,1)),'LabelHorizontalAlignment','left'); 
     y3 = xline(Conf_Int(m,2),'k-',num2str(Conf_Int(m,2)),'LabelHorizontalAlignment','right'); 
end


for warum_auch_immer = 1:2 % titel is printed in second call.... for whatever reason
    for m = 1:p.n_p
     %% Plot
     hold on
     subplot(ceil(p.n_p/2),2,m)
     xlabel(Names_para(m))
     ylabel('log$_{10}$(PL$(\mathbf{\theta})$)')

    end
end