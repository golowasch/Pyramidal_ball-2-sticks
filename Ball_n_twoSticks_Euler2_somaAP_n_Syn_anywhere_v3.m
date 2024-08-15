% 2023-12-04 JG
    % Basic ball-and-stick cell and simple Euler order 2 integrator 
    %   which returns (only) the voltages of all segments
    % All parameters in a structure p file create by cell_parameters_v2.m,
    % including:
    %   - building cable, 
    %   - setting current injection setting for subpulses,
    %   - synaptic settings (conductances, position in the cable, whether 
    %     trains or random, etc),
    %   - all biophysical parameters.
    % Active currents only in soma to produce action potentials, but can
    % expand to any segment:
    %   - All currents and variables are built and calculated in
    %     biophys_eqs_v1.m,
    %   - Synapses are alpha functions.
    % All conductances and Iapp are densities.
    % Currents (train, pulses or ramps) injected into soma by call to Iext.m function
        % Independent synapses (trains or single) into multiple branches, but one segment/branch at a time for now.
        
% initialize all parameters
%clear loc p1 p2 p3 p1r p.gsyn1 p.gsyn2 gs1 gs2 Vhm_Na Vhn_K gNa gKd
clear all 
spp.leg0 = [];
dt_syn = 1;
Cm = 1200;                                                           % nF/cm^2
g_syn = [0 2765 1345];
p = cell_parameters_v3(Cm(1), dt_syn(1));

dt_syn = 60:5:80
 dt_syn= [20]

Cm = 2000; %:300:1200;
p = cell_parameters_v3(Cm(1), dt_syn(1));

pf = NaN(length(Cm),length(dt_syn));

jj = NaN(1,length(Cm)*length(dt_syn)); p.dtsyn = NaN(1,length(Cm)*length(dt_syn));
c = NaN(1,length(Cm)*length(dt_syn)); d = NaN(1,length(Cm)*length(dt_syn)); 
pp1 = NaN(1,length(Cm)*length(dt_syn)); 
pp2 = NaN(1,length(Cm)*length(dt_syn)); 
pp3 = NaN(1,length(Cm)*length(dt_syn)); 
p1 = zeros(1,length(dt_syn));
p2 = zeros(1,length(dt_syn));
p3 = zeros(1,length(dt_syn));

VG_plots = figure(1); clf  

for k=1:length(Cm)
    p1 = NaN(1,length(dt_syn));
    p2 = NaN(1,length(dt_syn));
    p3 = NaN(1,length(dt_syn));
  
    for j=1:length(dt_syn)                                              % Increments in dt_syn
        for pn = 1:p.num_pulses
            tic
            p = cell_parameters_v3(Cm(k), dt_syn(j));                                  
            p.ns = build_syn_train(p);
            f = biophys_eqs_v3(p);  
            i_ext = Iext_v3(p, pn);                                     % Calculate the current (in nA/cm2)
            [V,I] = Euler2_v3(p, f, i_ext);                             % Integrate the voltages
                 
            tau = p.Rm*p.Cm; 
            tr = '[Cm: '+string(p.Cm)+', Tau: '+ string(tau)+']';
            spp.leg0 = [spp.leg0, tr];
                
            subplot(2,1,1)                        
            plot(p.t, V(1,:))
            hold on
% %%
       %     b=diff(V(1,:));
       % %   plot(p.t(1:25/p.dt), 100*b(1:25/p.dt)/(max(b)-min(b)))     % this is the derivative of V
       % 
       %     [a,aa]=max(b(1:20/p.dt))                      % Find the max of the derivative (a) and its position (b)
       %    aveb= mean(b(1,aa-5:aa+5))                    % Average of derivative around its max
       % 
       %    bb=b(ceil((p.ini+3)/p.dt))*1000                     % Derivative at 3msec after pulse onset
       %    bbb=bb*.8                                     % 80% of that
           j,pn
%%
            %legend(spp.leg0(1), 'dVm/dt', 'Syn')
            ylabel('Vm (mV)')
            ylim auto;            
            leg1 = 'Soma Vm + Derivative of Vm. Synapses at '+ string(p.pos_percent(2)) + ' & '+ string(p.pos_percent(3)) + '%' + spp.leg0(1);
            title(leg1);
            lim = axis; ylim_range = lim(4)-lim(3);
              
            subplot(2,1,2)
            hold on
            if sum(sum(p.g_syn))==0
                plot(p.t, i_ext*p.A(1))
                xlabel('Time [ms]')
                ylabel('Current [nA]')      
                ylim auto;
                leg1 = 'I_ext applied'; 
                title (leg1);
            else
                plot(p.t, p.ns(2,:)*p.g_syn(2,p.pos(2))*1000, p.t, p.ns(3,:)*p.g_syn(3,p.pos(3))*1000)
                xlabel('Time [ms]')
                ylabel('Conductance [nS]')
                ylim auto;
                leg1 = 'Synaptic conductance'; 
                title (leg1);     
            end

            % Find first threshold crossing
            x = find(V(1,:) > 0, 1); 
            if length(x)<1
                x=NaN;
            end            
            p1(j) = x*p.dt;
            p1r(j,pn) = x*p.dt;
            gs1(p.num_pulses*(j-1)+pn) = p.g_syn(2,p.pos(2));             % Create var to write g_syns to file
            gs2(p.num_pulses*(j-1)+pn) = p.g_syn(3,p.pos(3));                 % "

            Vhm_Na(p.num_pulses*(j-1)+pn) = p.Vh_m_Na(1);
            Vhn_Kd(p.num_pulses*(j-1)+pn) = p.Vh_n_Kd(1);
            gNa(p.num_pulses*(j-1)+pn) = p.g_Na(1);
            gKd(p.num_pulses*(j-1)+pn) = p.g_Kd(1);

            % %[pks,locs]=findpeaks(V(1,:),'MinPeakHeight', -10);
            %   [pks,locs]=findpeaks(V(1,:));
            % loc{j}=locs;
            % 
            % % This is to prevent errors in writing to table below:
            % if length(loc{j})<3
            %     loc{j}(length(loc{j})+1:3) = NaN;
            % end
            % 
            % p1(j)=loc{j}(1).*p.dt;       
            % p2(j)=loc{j}(2).*p.dt;         
            % p3(j)=loc{j}(3).*p.dt;      
            if pn==1 && j==1
                toc
                ['Estimated time to end ='  string(toc*length(dt_syn)*p.num_pulses/60) ' min']
            end
        end   
        

    end
    subplot(2,1,1)
    hold on
    scatter(p.ini+dt_syn, lim(3)+ylim_range*0.02, 14, 'filled', '^')        % Marke the time of the second syn
    
    syn_delay = dt_syn;    
    
    for m = 1:length(dt_syn)
       jj = length(dt_syn)*(k-1) + m;
       c(jj) = Cm(k);                    % These dummy variables are necessary because struct arrays can't easily be added on
       d(jj) = dt_syn(m); 
       if sum(p.syn_branch) == 1
          pp1(jj) = p1(m)-(syn_delay(m)+p.ini);
          %pf(k,:) = p1-syn_delay-p.ini;            % Delay from first syn input
          pf(k,:) = p1 - p.ini;                     % Delay from second syn input
       elseif sum(p.syn_branch) == 2
          pp1(jj) = p1(m)-(syn_delay(m)+p.ini);
          pf(k,:) = p1-syn_delay-p.ini;
       end
       if pp1(jj)<0; pp1(jj)=NaN; end
       
       pp2(jj) = p2(m)-syn_delay(m);
       if pp2(jj)<0; pp2(jj)=NaN; end
       
       pp3(jj) = p3(m)-syn_delay(m);
       if pp3(jj)<0; pp3(jj)=NaN; end   
    end
 
end

% % % dt_delays = figure(2);
% % % plot(syn_delay, pf(1:length(Cm),:))
% % % xlabel('dtsyn (msec)')
% % % ylabel('Time-to-threshold (msec)')
% % % xlim([0 dt_syn(end)*1.1])
% % % %ylim([0 30])
% % % legend(string((Cm)));

% for j=1:length(dt_syn) 
    % p.dV = (V(1,p.tspan * 0.9/p.dt)-V(1,50));                           % dV from Vrest to Vthresh
    % p.percent_Vthresh(1:length(Cm)) = p.dV/32.4*100;                    % 32.4 is only for this model 
    % p.Rin(1:length(Cm)) = p.dV/i_ext(p.tspan * 0.9/p.dt);               % Rin in MOhm
    % p.Rm(1:(length(dt_syn))) = p.Rm;
    % p.gsyn1(1:length(Cm)*(length(dt_syn))) = p.g_syn(2,p.pos(2));
    % p.gsyn2(1:length(Cm)*(length(dt_syn))) = p.g_syn(3,p.pos(3));
    % p.tau_syn(1:length(Cm)*(length(dt_syn))) = p.tau_syn;
% end

for j=1:length(dt_syn)
    p.Cm_w(p.num_pulses*(j-1)+1:p.num_pulses*j) = Cm;                          % Repeat for as many pulses
    p.dtsyn(p.num_pulses*(j-1)+1:p.num_pulses*j) = d(j);
    p.p1_w(p.num_pulses*(j-1)+1:p.num_pulses*j) = p1r(j,:);
    p.Vm_Na = Vhm_Na;
    p.Vm_Kd = Vhn_Kd;
    p.gNa=gNa;
    p.gKd=gKd;
    p.gsyn1 = gs1;
    p.gsyn2 = gs2;

    incr=10;
    for i=1:incr:p.num_pulses-incr                              % Count not NaN entries in groups of 'incr'
        a(floor(i/incr)+1,j)=sum(~isnan(p1r(j,i:i+incr)));
    end

% p.p2_w = pp2;
% p.p3_w = pp3;

end

load train        
sound(y,Fs)
 
%% Save data in files
selpath = datafolder
filename = makefilename(p);                     % indicate variable being changed
filename = [filename '_Active_PairedPulse-withvarISI'];
writedata(p, a, [filename '.xlsx'])

exportgraphics(VG_plots,[selpath '\Vm&gsyn.emf'],'ContentType','vector')

%exportgraphics(dt_delays,[filename '_Cm-vs_ISI.emf'],'ContentType','vector')
% % Prepare to write struct file

% p.t=[]; p.g_syn=[]; p.t_syn=[]; p.ns=[];        % Empty these to avoid crash of writestruct call (p.t is just too long)
% fn = [filename '_CellParams.xml'];              % they are written above anyway
% writestruct(p,fn)

function selpath=datafolder
    %% Select where data will go
    selpath = uigetdir(path, 'Select a directory to save your data');
    addpath(selpath);
    cd(selpath);
end

function filename = makefilename(p)
    % and get time and date to build a file name
    d=datetime('now');    
    d = char(d);
    d = replace(d,'2024','24');                                             % remove colons
    d = d(1:9);                     
    
    g1 = truncate(p.g_syn(2,p.pos(2)),3);
    g2 = truncate(p.g_syn(3,p.pos(3)),3);
   
    filename = ['Cm-effect_' d '_gsyn1_' g1 '_gsyn2_' g2 '_tau_syn_' char(string(p.tau_syn(1))) 'Cm_' char(string(p.Cm))];
end

function writedata(p, a, filename, selpath)   
     names = ["Cm"; "1st_ap del"; "g_Syn1"; "g_Syn2"; 'dt_syn'; 'gNa'; 'gKd'; 'V_m_Na'; 'V_m_Kd'];
     T = table(p.Cm_w', p.p1_w', p.gsyn1', p.gsyn2', p.dtsyn', p.gNa', p.gKd', p.Vm_Na', p.Vm_Kd', 'VariableNames', names);
     writetable(T, filename, 'WriteVariableNames', true)

     names = ["dt_syn"];
     T = table(a, 'VariableNames', names);
     writetable(T, ['Counts_per_10_Cm_' char(string(p.Cm)) '.xlsx'], 'WriteVariableNames', true)
     
    % names = ["Cm"; "1st_ap"; "2nd_ap"; "3rd ap"; "g_Syn1"; "g_Syn2"; "dt_syn"; "tau_syn"];
    % T = table(p.Cm_w', p.p1_w', p.p2_w', p.p3_w', p.gsyn1', p.gsyn2', p.dtsyn', p.tau_syn', 'VariableNames', names);
    % writetable(T, filename, 'WriteVariableNames', true)
    % 
    % Anum = length(p.Cm_w) + 4;    
    % pos = ['A' char(string(Anum))];
    % rows = ["Branch 1" ; "Branch 2"; "Branch 3"];
    % T = table(p.g_syn, 'RowNames', rows);
    % writetable(T, filename, 'WriteRowNames', true ,'Range', pos)
    % 
    % pos = ['B' char(string(Anum))];
    % a=1:1:max(size(p.g_syn));   
    % T = table(a);                             
    % writetable(T, filename, 'WriteVariableNames', false, 'Range', pos)
    % 
    % Anum = Anum - 1;    
    % names = "gSyn Dendrite segments";
    % pos = ['A' char(string(Anum))];
    % T = table("-", 'VariableNames', names);                                 % Just a title
    % writetable(T, filename, 'Range', pos)
    % 
    % s = size(p.g_syn);
    % Anum = Anum + s(1) + 5;                                                 % Go past the p.g_syn array
    % pos = ['A' char(string(Anum))];
    % T = table(p.t_syn, 'RowNames', rows);
    % writetable(T, filename, 'WriteRowNames', true, 'WriteVariableNames', false, 'Range', pos)
    % 
    % a = 1:1:length(p.t_syn);    
    % Anum = Anum - 1;                                                 % Go past the p.g_syn array
    % pos = ['B' char(string(Anum))];
    % T = table(a);                             
    % writetable(T, filename, 'WriteVariableNames', false, 'Range', pos)
    % 
    % Anum = Anum - 1;
    % names = "Synaptic intervals";
    % pos = ['A' char(string(Anum))];
    % T = table("-", 'VariableNames', names);
    % writetable(T, filename, 'Range', pos)    
   
end

function plotx_infs(p,f)
    figure
    v=-100:50;
    m_Na_inf = f.m_inf_Na(v,p.Vh_m_Na,p.k_m_Na);
    h_Na_inf = f.h_inf_Na(v,p.Vh_h_Na,p.k_h_Na);
    n_Kd_inf = f.n_inf_Kd(v,p.Vh_n_Kd,p.k_n_Kd);
    plot(v,m_Na_inf, v,h_Na_inf, v,n_Kd_inf, LineWidth=2)
    legend%('m_Na_inf', 'h_Na_inf', 'n_Kd_inf')
end

function r4=truncate(r,n)      % Truncate the value of r to n decimal points, returns a characters 
    r1=round(r);
    r2=string(r1);
    r3=string(r);
    l=strlength(r2);
    l=l+n+1;                   % n decimal positions and the .
    r3=char(r3);
    if r2==r3
        l=length(r3);
    end
    r4=r3(1:l);
end

function subplotting(p,spp,I,V,i_ext) 
  if spp.num_subplots==2
    % Always plot soma voltage 
    spp.leg1 = '                                Soma                                   Cm / Tau';
    subplot(2,1,1)
    plot(p.t,V(1,:))                                                        % Soma V
    ylabel('Vm (mV)')
    ylim auto;
    title (spp.leg1);
    legend(spp.leg0);
 
    % and always plot soma ionic currents
    subplot(2,1,2)
    plot(p.t,I(1:3,:))                                                        % Soma V
    xlabel('Time [ms]')
    ylabel('Current (pA)')
    ylim auto;
    title ('Ionic currents');
    hold on

  elseif spp.num_subplots==3 && spp.gsyn && spp.iext && spp.ionic==0
        spp.leg1 = '                                Soma                                   Cm / Tau';
        subplot(3,1,1)
        plot(p.t,V(1,:))                                                        % Soma V
        ylabel('Vm (mV)')
        ylim auto;
        title (spp.leg1);
        legend(spp.leg0);
       
    hold on

        subplot(3,1,2)
        plot(p.t, I(4,:))
        ylabel('Synaptic Current [pA]')
        ylim auto;
        leg1 = 'Isyn'; %, Dendrite length: '+ string(p.dend_l(p.syn_branch(end))*1e4) + 'um, Synapse position: ' + string(p.dend_l(p.syn_branch(end))/(p.pos/max(p.pos))*1e4) + 'um';
        leg2 = [];%'  gSyn= ' + string(1000*p.g_syn(p.syn_branch(end), p.pos(end))) + ' nS';
        title (leg1);
        hold on

        subplot(3,1,3)
        plot(p.t, i_ext*1e3)
        xlabel('Time [ms]')
        ylabel('Current [pA]')
        ylim auto;
        leg1 = 'Soma injection'; max(abs(i_ext)*1e3);
        title (leg1);
        hold on

  elseif spp.num_subplots==3 && spp.gsyn==1 && spp.iext==0
        spp.leg1 = '                              Soma                                 Cm / Tau';
        subplot(3,1,1)
        plot(p.t,V(1,:))                                                        % Soma V
        ylabel('Vm (mV)')
        ylim auto;
        title (spp.leg1);
        legend(spp.leg0);
        hold on

        subplot(3,1,2)
        plot(p.t,I(1:3,:))                                                        % Soma V
        ylabel('Current (pA)')
        ylim auto;
        title ('Ionic currents');
        hold on
    
        subplot(3,1,3)
        plot(p.t(1:length(p.t)), I(4,:))
        xlabel('Time [ms]')
        ylabel('Current [pA]')
        ylim auto;
        leg1 = []; %'Isyn, Dendrite length: '+ string(p.dend_l(p.syn_branch(end))*1e4) + 'um, Synapse position: ' + string(p.dend_l(p.syn_branch(end))/(p.pos/max(p.pos))*1e4) + 'um';
        leg2 = [];%'  gSyn= ' + string(1000*p.g_syn(p.syn_branch(end), p.pos(end))) + ' nS';
        title (leg1+leg2);
        hold on

    elseif spp.num_subplots==3 && spp.gsyn==0 && spp.iext==1
        spp.leg1 = '                                Soma                                 Cm / Tau';
        subplot(3,1,1)
        plot(p.t,V(1,:))                                                        % Soma V
        ylabel('Vm (mV)')
        ylim auto;
        title (spp.leg1);
        legend(spp.leg0);
        hold on
    
        subplot(3,1,2)
        plot(p.t,I(1:3,:))                                                        % Soma V
        ylabel('Current (pA)')
        ylim auto;
        title ('Ionic currents');
        hold on

        subplot(3,1,3)
        plot(p.t,i_ext*1e3)
        xlabel('Time [ms]')
        ylabel('Current [pA]')
        ylim auto;
        leg1 = 'Soma injection'; max(abs(i_ext)*1e3);
        title (leg1);
        hold on

    elseif spp.num_subplots==4
        spp.leg1 = '                               Soma                                 Cm / Tau';
        subplot(4,1,1)
        plot(p.t,V(1,:))                                                        % Soma V
        ylabel('Vm (mV)')
        ylim auto;
        title (spp.leg1);
        legend(spp.leg0);
        hold on
    
        subplot(4,1,2)
        plot(p.t,I(1:3,:))                                                        % Soma V
        ylabel('Current (pA)')
        ylim auto;
        title ('Ionic currents');
        hold on

        subplot(4,1,3)
        plot(p.t(1:length(p.t)), I(4,:))   
        ylabel('Current [pA]')
        ylim auto;
        leg1 = 'Isyn, Dendrite length: '+ string(p.dend_l(p.syn_branch(end))*1e4) + 'um, Synapse position: ' + string(p.dend_l(p.syn_branch(end))/(p.pos/max(p.pos))*1e4) + 'um';
        hold on       

        subplot(4,1,4)
        plot(p.t,i_ext*1e3)
        xlabel('Time [ms]')
        ylabel('Current [pA]')
        ylim auto;
        leg1 = 'Soma injection'; max(abs(i_ext)*1e3);
        title (leg1);
        hold on
    end 
end

%% Build a synaptic train
function pns=build_syn_train(p)
    p.ns = zeros(p.branch_num+1,length(p.t));
    for br = 1:p.branch_num + 1
        s = zeros(1,max(length(p.t)));
        ns = zeros(1,length(p.t));
        for j=1:length(p.t_syn(br,:))
            s = heaviside(p.t-p.t_syn(br,j)) .*(p.t-p.t_syn(br,j))./p.tau_syn .*exp(-p.t/p.tau_syn);  % alpha function as synapse
            if max(s)>0; s = s./max(s); end
            ns = ns+s;
        end
        pns(br,:) = ns;    
        pns(br,:) = pns(br,:)/length(p.t_syn(br,:));
    end
 end