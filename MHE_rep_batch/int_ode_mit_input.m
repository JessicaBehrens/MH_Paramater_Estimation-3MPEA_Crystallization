function [x_int] = int_ode_mit_input(Para,x0_tH,tH,idx_feed_tH,x0_feed,p)


  if ~isempty(idx_feed_tH)
        % time left
        t_new_feed = tH(idx_feed_tH:end); 
           if  length(t_new_feed) ==1
               x_int_feed = x0_feed';
           else
            %% generate new data 
            % solve ode
            [~,x_int_feed] =  ode15s(@(t,c)ode_system_model(t,c,Para,p),t_new_feed,x0_feed,p.opt);
              if length(t_new_feed) == 2 % solver calculates to many points in case of tH = 2
                  x_int_feed(2:end-1,:) = []; % reduce to point in tH
              end
           end
          t_no_feed = tH(1:idx_feed_tH-1);
          if length(t_no_feed) == 1
               x_int_no_feed = x0_tH';
          elseif isempty(t_no_feed)
              x_int_no_feed = [];
          else
               % integrate ode system
               [t,x_int_no_feed] = ode15s(@(t,c)ode_system_model(t,c,Para,p),t_no_feed,x0_tH,p.opt);
              if length(t_no_feed) == 2 % solver calculates to many points in case of tH = 2
                  x_int_no_feed(2:end-1,:) = []; % reduce to point in tH
              end
          end
          x_int = [x_int_no_feed;x_int_feed];
     else
       if length(tH) == 1
         x_int = x0_tH;
       else
        % integrate ode system
        [t,x_int] = ode15s(@(t,c)ode_system_model(t,c,Para,p),tH,x0_tH,p.opt);
       end
     end
     if length(tH) == 2 % solver calculates to many points in case of tH = 2
         x_int(2:end-1,:) = []; % reduce to point in tH
     end


end