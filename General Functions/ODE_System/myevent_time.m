function [values,isterminal,direction] = myevent_time(t,x,tstart)
% stop integration if ode iteration exceeds ...
 values = toc(tstart) < 10; % ... seconds
 
 isterminal = true(size(values));
 direction = zeros(size(values));
end