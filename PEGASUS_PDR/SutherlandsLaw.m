function [mu, k] = SutherlandsLaw(Tatm)
    T = Tatm; 
    mu_o = 1.716e-5; %Air viscosity 
    To = 273;
    S_mu = 111;
    mu = mu_o * (T/To) ^(3/2) *(To+S_mu) / (T+S_mu) ; 

    ko = 0.0241; 
    S_k = 194; 
    k = ko * (T/To) ^(3/2) *(To+S_k) / (T+S_k) ; 
end
