function [f] = biophys_eqs_v3(p)
% This function creates a structure file f containing biophysical
%   functions (currents and activation variables)
% A synaptic train is built here (from array of specific or random start
%   times)
% All currents return as: uS * (mV-mV) = nA.
    
% I_Leak
    f.Ileak = @(v,g,E) g.*(v-E);                                            

% I_Kd:
    f.n_inf_Kd = @(v,Vh,k) 1./(1 + exp((v - Vh) ./k));
    f.tau_Kd = @(v,tau_lo,tau_hi,Vh,k) tau_lo + tau_hi .*(1 - f.n_inf_Kd(v,Vh,k));
    f.dn_dt = @(v,n,tau_lo_Kd,tau_hi_Kd,Vh,k) (f.n_inf_Kd(v,Vh,k) - n) ./f.tau_Kd(v,tau_lo_Kd,tau_hi_Kd,Vh,k);
    f.I_Kd = @(v,g,n) g .*(n .^4) .*(v - p.E_K);

% I_Na
    f.m_inf_Na = @(v,Vh,k) 1 ./ (1 + exp((v - Vh) ./k));
    f.h_inf_Na = @(v,Vh,k) 1 ./ (1 + exp((v - Vh) ./k));
    f.I_Na = @(v,g,m_Na,h_Na) g .*(m_Na .^3) .*h_Na .*(v-p.E_Na);
    f.dm_dt = @(v,m,Vh,k,tau_m) (f.m_inf_Na(v,Vh,k) - m) ./tau_m;
    f.dh_dt = @(v,h,Vh,k,tau_h) (f.h_inf_Na(v,Vh,k) - h) ./tau_h;
    f.I_Na = @(v,g,m,h) g .*(m .^3) .*h .*(v - p.E_Na);

% I_A
    f.m_inf_A = @(v,Vh,k) 1 ./ (1 + exp((v - Vh) ./ k));
    f.h_inf_A = @(v,Vh,k) 1 ./ (1 + exp((v - Vh) ./ k));
    f.I_A = @(v,g,m_A,h_A) g .*(m_A .^2) .*h_A .*(v - p.E_K);

% I_MI
    f.m_inf_MI = @(v,Vh,k) 1 ./ (1 + exp((v - Vh) ./ k));
    f.I_MI = @(v,g,m_MI) g .*(m_MI .^2) .*(v - p.E_MI);

% I_MIT
    f.m_inf_MIT = @(v,Vh,k) 1 ./ (1 + exp((v - Vh) ./ k));
    f.h_inf_MIT = @(v,Vh,k) 1 ./ (1 + exp((v - Vh) ./ k));
    f.I_MIT = @(v,g,m_MIT, h_MIT) g .*(m_MIT .^3) .*h_MIT .* (v - p.E_MIT);
    f.taum_MIT = @(v) 2  + 150 ./ cosh(0.25 .* (v + 20));
    f.tauh_MIT = @(v) 100  + 600 ./ cosh(0.1 .* (v + 30)) + 200 ./ (1 + exp(-0.2 .* (v + 30)));

% I_syn
    f.I_syn = @(v,t,br,j) p.g_syn(br,j) .*p.ns(br,t) .*(v - p.E_syn);              % uS * (mV-mV) = nA    

end    % ------------------ end of function --------------------