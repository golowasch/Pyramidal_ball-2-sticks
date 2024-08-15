function [p] = cell_parameters_v4(Cm, dt_syn)
% This function creates parameters for a ball-and-multiple-stick model and passes parameters in the 
% structure file 'p'.
% Branch 1, and segment 1 of branch 1 (1,1) is the soma, First branch is numbered 2, with segments 1,2, ...
% Define independent synaptic inputs (random or regular trains, single) for each branch; one synapse/branch only.

%% Time stuff
    p.ti = 0;
    p.dt = 0.0008;
    p.tf = 105;
    p.t = p.ti:p.dt:p.tf;
    p.tspan = p.tf-p.ti;

%% Current Parameters 
    p.do_subpulses = 0;                                                     % 0 or 1
    p.subpulse_dur = 30;                                                    % sub-pulse duration
    p.num_subpulses = 1;                                                    % number of pulses for a train to test Rin
    p.iapp_subpulse = 100;                                                  % nA/cm2
    p.do_pulses = 1;                                                        % 0 or 1
    p.num_pulses = 2;
    if p.do_pulses == 0; p.num_pulses=1; end
 
% Check if subpulses fit in time window (p.tf-p.ti)
    if p.subpulse_dur * 2 > p.tspan*0.8                                     % *2 for a full subpulse cycle; *0.8 to leave 10% edges
       p.subpulse_dur = 0.5*p.tspan*0.8;
       p.num_subpulses = 1;
    elseif p.subpulse_dur*2 < p.dt * 100
       p.subpulse_dur = p.dt*100;
       p.num_subpulses = round(p.tspan*0.8/(p.subpulse_dur*2));
    end
    
    if p.num_subpulses > (0.5*p.tspan*0.8)/p.subpulse_dur                   % Warning
       p.num_subpulses = floor(0.5*p.tspan*0.8/p.subpulse_dur);
    elseif p.num_subpulses <= 1
        p.num_subpulses = 1;
    end
    
%% Specific electrical properties
    p.Ra = 80 * 1e-6;                           % Ra in MOhm*cm; makes axial conductances into uS
    p.Rm = 6000 * 1e-6; % 12500                       % Rm in MOhm.cm^2 / default = 40000; this produces g_leak in uS/cm^2
    p.Cm = Cm;                                  %750;                          % uF/cm2

%% Geometry 
    % Remember: 1um = 1e-4 cm
    p.branch_num = 2;                                                       % Nr of branches: 0 - inf
    p.withcable= 1;                                                         % 0 for no dendrite attached
    
    p.A = zeros(p.branch_num+1, 1);                                         % Areas of each branch
    p.cm = zeros(p.branch_num+1, 1);                                        % Branch capacitance
    p.rad = zeros(p.branch_num+1, 1);                                       % Radius of each branch
    p.dend_l = zeros(p.branch_num+1, 1);                                    % Dendrite length
    p.lambda = zeros(p.branch_num+1, 1);                                    % Branch lambda
    p.lenseg = zeros(p.branch_num+1, 1);                                    % Segment length
    p.gleak = zeros(p.branch_num+1, 1);                                     % Branch gleak 
    p.ra = zeros(p.branch_num+1, 1);                                        % Branch axial resistance 
    p.nums = zeros(p.branch_num+1, 1);                                      % Number of segments per branch

% The soma
    p.rad(1) = 15 * 1e-4;                                                   % soma radius in microm * cm/microm = cm
    p.A(1) = 4 * pi * p.rad(1)^2;                                           % sphere's surface
    p.nums(1) = 1;                                                          % Soma

% Here we deal with branch and leave soma out. Soma is i=1.
    p.rad(2:end) = 1 * 1e-4;                                                % dendrite radius in um * cm/um = cm
    p.dend_l(2) = 180 * 1e-4;                                          % dendrite length in um * cm/um = cm
    p.dend_l(3) = 150 * 1e-4;
    p.A(2:end) = p.dend_l(2:end).*2.*pi.*p.rad(2:end);                      % Areas of the whole dendrites

% Passive electrical Properties, and dendrite compartments
    p.lambda(2:end) = sqrt(p.Rm.*p.rad(2:end)*2./p.Ra/4);                   % length constant (cm)
    p.nums(2:end) = ceil(p.dend_l(2:end) ./ (p.lambda(2:end)/40));          % Round # of segments up @ 40/Lambda
    p.lenseg(2:end) = p.dend_l(2:end) ./ p.nums(2:end);
    p.ra(2:end) = p.Ra .* p.lenseg(2:end)./(pi.*p.rad(2:end).^2);           % ra of each segment (MOhm)
    
    p.A(:) = p.A(:) ./ p.nums(:);                                           % Adjust all areas for segment length; for soma p.nums(1)=1
    p.cm(:) = p.Cm .* p.A(:);                                               % dendritic segment cm (nF)

% I_Kd
    p.g_Kd = zeros(p.branch_num+1, 1); 
    p.Vh_n_Kd = zeros(p.branch_num+1, 1);                              
    p.k_n_Kd = zeros(p.branch_num+1, 1);                               
    p.tau_lo_Kd = zeros(p.branch_num+1, 1);                            
    p.tau_hi_Kd = zeros(p.branch_num+1, 1);                            

% I_Na
    p.g_Na = zeros(p.branch_num+1, 1);                                 
    p.Vh_m_Na = -zeros(p.branch_num+1, 1);                            
    p.k_m_Na = zeros(p.branch_num+1, 1);                               
    p.Vh_h_Na = zeros(p.branch_num+1, 1);                             
    p.k_h_Na = zeros(p.branch_num+1, 1);                               
    p.tau_m_Na = zeros(p.branch_num+1, 1);                             
    p.tau_h_Na = zeros(p.branch_num+1, 1);

% I_A
    p.g_A = zeros(p.branch_num+1, 1);                             
    p.Vh_m_A = zeros(p.branch_num+1, 1);                         
    p.k_m_A = zeros(p.branch_num+1, 1);                           
    p.Vh_h_A = zeros(p.branch_num+1, 1);                         
    p.k_h_A = zeros(p.branch_num+1, 1);                          
    p.tau_m_A = zeros(p.branch_num+1, 1);                         
    p.tau_h_A = zeros(p.branch_num+1, 1);

%% I_MI:
    p.g_MI = zeros(p.branch_num+1, 1);                                
    p.Vh_m_MI = zeros(p.branch_num+1, 1);                                        
    p.k_m_MI = zeros(p.branch_num+1, 1);                              
    p.tau_m_MI = zeros(p.branch_num+1, 1);                          

%% I_MIT:
    p.g_MIT = zeros(p.branch_num+1, 1);
    p.Vh_m_MIT = zeros(p.branch_num+1, 1);
    p.k_m_MIT = zeros(p.branch_num+1, 1);
    p.Vh_h_MIT = zeros(p.branch_num+1, 1);
    p.k_h_MIT = zeros(p.branch_num+1, 1);
    p.tau_m_MIT = zeros(p.branch_num+1, 1);
    p.tau_h_MIT = zeros(p.branch_num+1, 1);

%% Set default values
 % I_leak
    p.E_leak = -69;
    p.gleak = (1./p.Rm).*p.A(:);                                             % leak conductance (uS)
    
 %% I_Kd:
    p.E_K = -90;
    p.g_Kd(:) = 17500.*p.A(:);%.*(1.02-0.04*rand(1,1));                      % 3500 for Cm1200
    p.Vh_n_Kd(:) = -16.*(1.05-0.1*rand(1,1));                              % in mV / -12 for Vthresh=-35, -32 for Vthresh=-56
    % p.Vh_n_Kd(:) = -16;                              % in mV / -12 for Vthresh=-35, -32 for Vthresh=-56
    p.k_n_Kd(:) = -7;                                                      % in mV / -10 for Vthresh=-35
    p.tau_lo_Kd(:) = 1;                                                     % in msec / 1 msec for a.p. duration 2msec
    p.tau_hi_Kd(:) = 12;                                                    % in msec / 12 msec for a.p. duration 2msec
    
 %% I_Na:
    p.E_Na = 60;
    p.g_Na(:) = 35000.*p.A(:);%.*(1.02-0.04*rand(1,1));                     % 20000 .* p.A(:); for Cm1200
    p.Vh_m_Na(:) = -33.*(1.05-0.1*rand(1,1));                              % in mV / Default = -25 for Vthresh =-35mV
    % p.Vh_m_Na(:) = -33;
    p.k_m_Na(:) = -8.0;                                                     % in mV / Default = -5 for Vthresh =-35mV
    p.Vh_h_Na(:) = -39;                                                     % in mV / Default = -28 for Vthresh =-35mV, -50 for Vthresh=-56
    p.k_h_Na(:) = 3.8;                                                      % in mV / Default = 3.7 for Vthresh =-35mV
    p.tau_m_Na(:) = 0.01;                                                   % 0.01 for a.p. of 2 msec
    p.tau_h_Na(:) = 0.6;                                                    % 0.6 for a.p. of 2 msec
    
 %% I_A:
    p.g_A(:) = 00 .* p.A(:);                                                % Dendrites. Default ~300
    p.Vh_m_A(:) = -23;                                                      % for m / Default = -18 (or -25)
    p.k_m_A(:) = -15.5;                                                     % Default = -12.5
    p.Vh_h_A(:) = -43;                                                      % for h / Default = -28
    p.k_h_A(:) = 12;                                                        % Default = 7.7
    p.tau_m_A(:) = 2;                                                       % time constants in msec / Default = 2
    p.tau_h_A(:) = 40;                                                      % Default = 5
    
 %% I_MI:
    p.g_MI(:) = 0 .* p.A(:);                                                % Default ~300
    p.E_MI = 7;
    p.Vh_m_MI(:) = -30;                                                     % for m / Default = -18 (or -25)
    p.k_m_MI(:) = -55;                                                      % Default = -12.5
    p.tau_m_MI(:) = 2;                                                      % time constants in msec / Default = 2
    
 %% I_MIT:
    p.g_MIT(:) = 00 .* p.A(:);                                              % Default ~300
    p.E_MIT = 7;
    p.Vh_m_MIT(:) = -13;                                                    % for m / Default = -18 (or -25)
    p.k_m_MIT(:) = -4;                                                      % Default = -12.5
    p.Vh_h_MIT(:) = -57;                                                    % for h / Default = -28
    p.k_h_MIT(:) = 11;                                                      % Default = 7.7
    p.tau_m_MIT(:) = 10;                                                    % time constants in msec / Default = 2
    p.tau_h_MIT(:) = 40;                                                    % Default = 5
    
 %% Synapses (for now these are external inputs for chemical and cell-to-cell if electrical)
    p.syn_branch = zeros(p.branch_num+1, 1);                                % Build an array for branches with synapses
    p.pos = zeros(p.branch_num+1, 1);
    syn_seg = zeros(p.branch_num+1, max(p.nums)+1);                         % Build an array for synaptic positions
    p.g_syn = zeros(p.branch_num+1, max(p.nums)+1);                         % Synaptic conductances   
       
    p.syn_branch = [0; 1; 1];                                               % Branch with syn = 1; branch 1 is soma
    p.pos = [1; 50; 50];                                                  % Position of syn along ea branch: 1 is base, 100 is tip,   
    p.pos_percent = p.pos;
   
    % g_syn = [0 2820 2810 ];      % [0 2820 2810] is Vthr for tau_syn=5 msec, apical and basal syns at pos=50%; Vthresh ~-41 mV, Cm=600nF
    % g_syn = [0 2150 2160];      % [0 2150 2160] is ~80% of the maximal slope that triggers an a.p. independently for each synapse and w/Cm=600nF, pos=66%

    % g_syn = [0 3730 3720];      % [0 3730 3720] is Vthr for tau_syn=5 msec, apical and basal syns at pos=50%; Vthresh ~-41 mV, Cm=1200nF
    % g_syn = [0  2830 2830];      % [0 2830 2830] is ~80% of the maximal slope that triggers an a.p. independently for each synapse and w/Cm=1200nF, pos=66%
   
    g_syn = [0 6900 6750];      % [0 6900 6750] is Vthr for tau_syn=5 msec, apical and basal syns at pos=50%, respectively; Vthresh =-41 mV, Cm=4000nF
    g_syn = [0 5550 5400];      % [0  5550 5400] is ~80% of the maximal slope that triggers an a.p. independently for each synapse and w/Cm=4000nF, pos=50%

    g_syn = [0 4180 4120];      % [0 4180 4120] is Vthr for tau_syn=5 msec, apical and basal syns at pos=50%, respectively; Vthresh =-41 mV, Cm=2000nF
    g_syn = [0 3240 3250];      % [0 3240 3250] is ~80% of the maximal slope that triggers an a.p. independently for each synapse and w/Cm=2000nF, pos=50%

     % g_syn = [0 0 0];

    p.E_syn = 0;                % in mV
    
    for j=1:length(p.syn_branch)
        if p.syn_branch(j) > 1
            p.syn_branch(j) = 1;
        elseif p.syn_branch(j) < 0
            p.syn_branch(j) = 0;
        end
    end
    p.pos = p.syn_branch .* p.pos;                                          % This ensureS that if a branch has no syn, p.pos=0
                                                             
    if length(p.syn_branch) > p.branch_num + 1
        p.syn_branch = p.branch_num + 1;                                    
    elseif length(p.syn_branch) < 1
        p.syn_branch = 1;
    end
    
    p.pos = ceil(p.pos.*p.nums./101)+1;                                     % Put syn in here. Position '1' is for soma
   
    for j = 1:length(p.syn_branch)
        syn_seg(j,p.pos(j)) = p.syn_branch(j);                              % 1s where there is a syn, 0s elsewhere
    end
    
    for j = 1:length(p.syn_branch)
        p.g_syn(j,:) = g_syn(j).*syn_seg(j,:);
        p.g_syn(j,:) = p.g_syn(j,:) .*p.A(j); 
    end
    
    p.syn_train = 1;                                                        % if 1 use array of start times p.t_syn; 2 use random train in msec
    p.tau_syn = 5;                                                           
    train_len = 2;                                    
    p.ini = 5;
    p.t_syn = zeros(p.branch_num+1,train_len);
    
    if p.syn_train ==1                                                      % This is for train of non-random start times
        for j = 1:train_len
            if p.syn_branch(2)==1 && p.syn_branch(3)==0
                p.t_syn(2,j) = p.ini + (j-1)*dt_syn;
            elseif p.syn_branch(3) == 1 && p.syn_branch(2)==0
                p.t_syn(3,j) = p.ini + (j-1)*dt_syn;
            else
                p.t_syn(2,j) = p.ini ;                                      % Apical dendrite first
                p.t_syn(3,j) = p.ini + dt_syn;
            end
        end              
        p.dt_syn = ['_even_' char(string(dt_syn(1)))];
    elseif p.syn_train == 2 && Cm > p.Cm(1)                                 % For train with random start times
        n = train_len;
        a = rand(n); a = sort(a(1,:));
        p.t_syn(2,:) = p.tspan*0.1+p.tspan*0.8 * a;                         % Random train
        a = rand(n); a = sort(a(1,:));
        p.t_syn(3,:) = p.tspan*0.1+p.tspan*0.8 * a;                         % Random train
        p.dt_syn = ['_rdm_' char(string((p.t_syn(n)-p.t_syn(1))/(n-1)))];   % return average dt_syn
    else
        p.t_syn = NaN;
    end
        
end