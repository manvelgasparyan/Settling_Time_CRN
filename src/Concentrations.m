function [t,x] = Concentrations (k, K, S, mu, T, x_0)
%-------------------------------------------------    
    function [dx_dt] = Model (~, x)
       %---------------------------
       %Substrate composition matrix
       Delta  = - min(S,0);
       %---------------------------
       %Diagonal of rate constants
       Gamma_k = diag(k);
       %---------------------------
       %Substrate expression function
       phi = exp(Delta'*log(x));
       %---------------------------
       %Rational terms in reaction rates
       g = [1/(1 + x(1)/K(1) + x(2)/K(2) + x(3)/K(3));
            1/(1 + x(1)/K(1) + x(2)/K(2) + x(3)/K(3));
            1/(1 + (x(3)*x(4))/(K(4)*K(5)) + x(3)/K(4) + x(4)/K(5));
            1/(1 + (x(3)*x(4))/(K(4)*K(5)) + x(3)/K(4) + x(4)/K(5))];
       %---------------------------
       %Diagonal of rational terms
       Gamma_g = diag(g);
       %---------------------------
       %Reaction rates
       v = Gamma_k*Gamma_g*phi;
       %---------------------------
       %The balance laws
       dx_dt = S*v;
    end
%-------------------------------------------------
%Division of the time interval
t_span = 0:T/mu:T;
%-------------------------------------------------
x_0(x_0==0) = realmin;
%-------------------------------------------------
f = @(tau,xi)Model(tau, xi);
%-------------------------------------------------
%Species concentrations
[t,x] = ode23tb(f, t_span, x_0);
%-------------------------------------------------
x = x';
end
