function [t,x] = concentrations (k, K, S, mu, T, x_0)
%-------------------------------------------------     
%Division of the time interval
t_span = 0:T/mu:T;
%-------------------------------------------------
x_0(x_0==0) = realmin;
%-------------------------------------------------
f = @(tau,xi)mathematical_model(tau, k, K, xi);
%-------------------------------------------------
%Species concentrations
[t,x] = ode23tb(f, t_span, x_0);
%-------------------------------------------------
x = x';
end
