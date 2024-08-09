%% Artificial data generation
% ---------------------------------------------------------------------- Stats ----------------------------------------------------------------------

% end time
t_end = 30;

% sampling rate
t_sample = 30; % (every t_sample min)
p.Ts = t_sample/60; % h

% number of samples
p.N_sample = t_end/(p.Ts);

% timespan
p.tspan = linspace(0,t_end*3600,p.N_sample+1);

% ------------------------------------------------------------- integrate system ------------------------------------------------------------------------------------------------------

% initial condition
x0 = initial_condition(p);
p.n_s = length(x0); % numer of states

% solve ode
[~,x_int] =  ode15s(@(t,c)ode_system_model(t,c,para,p),p.tspan,x0,p.opt);

% -------------------------------------------------------- create measurement eq function handle ------------------------------------------------------------------------------------------------------

%% measurement equation
p.h = @(x)Measurement_Eq(x,p); % create fucntion handle
p.n_y = length(p.h(x0)); % number of output equations

% ------------------------------------------------------------- artificial noisy measurement ------------------------------------------------------------------------------------------------------

% standard deviation
proz_std_mess = 0.05; % in %

% states
x_data = x_int;
p.real_x_data = x_data;

% mean of measurements
mean_measurement = mean(Measurement_Eq(x_int,p),2);

% create noisy data
rand_num = randn(p.n_y,length(p.tspan));
y_data = p.h(x_data) + proz_std_mess.*rand_num.*mean_measurement;
y_data(y_data<0) = 0; % non-negative measurements

% save standard deviation
p.sigma_m = proz_std_mess*mean_measurement;
p.sigma_m(p.sigma_m==0) = 1e-24;