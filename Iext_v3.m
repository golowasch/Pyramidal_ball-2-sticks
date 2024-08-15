function [H,p] = Iext_v3(p,idx)  
% External current injection function. t is time, p is a structure
% parameter file, and idx is an index for pulse number
% Iext is made of 5 epochs: 10% of tspan is baseline at the begining and at the end
% A, B and C are epochs, each of which has amplitude and can be pulse or ramp
% If ramp, each epoch starts at the end of the previous.
% If p.subpulses = 1, override pulse structure and do a train that last t
% minus 10% at the begining and at the end

% Future: inject a current in each segment. For this define all I_ext
% params for each segmentand then run a for-loop through all segments and
% create an I_ext output for each. Redefine I_ext within the integrator for
% this.
% The variable H could be set to vary between [0 1] and be used to ramp
% conductances or currents or anything.

nsp = zeros(1, length(p.t));                                             % Array for subpulse structure
pls =  zeros(1, length(p.t));
H = zeros(1, length(p.t));
t = p.t;

pulse = @(t,ti,tf) heaviside(t-ti).*heaviside(tf-t);     

if p.do_subpulses==1
    % build a train of subpulses (e.g. to test Rin)
    sptrain = (p.tf-p.ti)/10: p.subpulse_dur*2 : p.tf*0.9-p.subpulse_dur;                                     % A train is built

    for lspt = 1:length(sptrain)
        sp = @(t) heaviside(t-sptrain(lspt)) .*heaviside(sptrain(lspt)+p.subpulse_dur-t);              % the train of pulses
        sp = sp(t)/max(sp(t));
        nsp = nsp+sp;
    end
    nsp = nsp*p.iapp_subpulse;   

%% Apply pulses or ramps (5 epochs per episode)
elseif p.do_pulses==1 && p.num_pulses>0

    tAi = p.tspan / 20;                                                     % First epoch starts at 5% from the begining
    if tAi<0.05*p.tspan; tAi = 0.05*p.tspan; end                            % Warning. Not earlier than 10% of p.tspan
    tAf = p.tspan * 0.1;                                                    % and ends here
    if tAf<=tAi; tAf = tAi+0.1*p.tspan; end                                 % Warning. Pulse must end after it begins
    tBf = p.tspan * 0.767;                                                    % End of second epoch
    if tBf<=tAf; tBf = tAf+0.1*p.tspan; end                                 % Warning. Pulse must end after it begins
    tCf = p.tspan * 0.82;                                                    % End of third epoch
    if tCf<=tBf; tCf = tBf+0.1*p.tspan; end                                 % Warning. Pulse must end after it begins
    tDf = p.tspan * 0.9;                                                    % Last epoch starts at 10% from the end
    if tDf<=tCf
        tDf=tCf+0.1*p.tspan;                                                % Warning. Pulse must end after it begins
    elseif tDf > 0.95*p.tspan
        tAi=0.05*p.tspan; tAf=0.3*p.tspan; tBf=0.5*p.tspan; tCf=0.7*p.tspan; tDf=0.95*p.tspan;
    end
       
    % Amplitudes (all in nA/cm2)
    A_lo = 0;                                           % Start value
    B_lo = 00;
    C_lo = 0;
    D_lo = 0;
     
    p.A_stp = 00;       % 0
    p.B_stp = 00;      %200                                   
    p.C_stp = 00;       %200
    p.D_stp = 00;       %200
    % Current step size; if there is a step >< 0, then the high level is x_stp if p.num_pulses = 1, else it increments

    p_A = 1;                                            % Activate pulses (0 or 1)
    p_B = 1;
    p_C = 1;
    p_D = 1;
    
    r_A = 0;                                            % Activate ramps (0 or 1)
    r_B = 0;
    r_C = 0;
    r_D = 0;
    
    if p_A==1 && r_A==1
        r_A = 0;
    elseif p_B==1 && r_B==1
        r_B = 0;
    elseif p_C==1 && r_C==1
        r_C = 0;
    elseif p_D==1 && r_D==1
        r_D = 0;
    end
    
    epi_base = A_lo.*pulse(t,t(1),tAi);                                      % The baseline
    epi_A = p_A.*(A_lo + p.A_stp*(idx)).*pulse(t,tAi,tAf)...
        + r_A.*(A_lo + p.A_stp*(idx).*pulse(t,tAi,tAf).*(t-tAi)./(tAf-tAi));
    
    epi_B = p_B.*(B_lo + p.B_stp*(idx)).*pulse(t,tAf,tBf)...
        + r_B.*(A_lo+p.A_stp.*(idx) + p.B_stp*(idx).*(t-tAf)./(tBf-tAf)).*pulse(t,tAf,tBf);
    
    epi_C = p_C.*(C_lo + p.C_stp*(idx)).*pulse(t,tBf,tCf)...
        + r_C.*(B_lo+p.B_stp.*(idx) + p.C_stp*(idx).*(t-tBf)./(tCf-tBf)).*pulse(t,tBf,tCf);
            
    epi_D = p_D.*(D_lo + p.D_stp*(idx)).*heaviside(t-tCf).*heaviside(tDf-t)...
        + r_D.*(C_lo+p.C_stp.*(idx) + p.D_stp*(idx).*(t-tCf)./(tDf-tCf)).*pulse(t,tCf,tDf);
            
    pls = epi_base + epi_A + epi_B + epi_C + epi_D;
    
 
end
H = (pls + nsp);

% figure(4); plot(p.t, H)
      
end   % -------------------- End function ----------------------