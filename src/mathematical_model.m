function [dx_dt] = mathematical_model (~, k, K, x)
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
       g = Denominator(k,K,x);
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
