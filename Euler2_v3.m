function [V,I] = Euler2_v3(p,f,i_ext)
% Simple Euler order 2 integrator for soma and multiple branches, but branches emerge only from soma for now.
% Returns V and currents (ionic and synaptic) (I). It needs a short dt when a cable is attached (e.g.<=0.005ms).
% Consider changing this to an internal ODE solver
       
%% Initializing variables for the ODE solver (Euler order 2)
   % make as many V arrays as soma (V1) + number of dendritic segments
        
    V = zeros(sum(p.nums),length(p.t));            
    n = zeros(1,length(p.t));
    m = zeros(1,length(p.t));
    h = zeros(1,length(p.t));
    I = zeros(3+length(p.branch_num),length(p.t));
    k1v = zeros(1,sum(p.nums));                                             % initialize integration variables
    k2v = zeros(1,sum(p.nums));
    av = zeros(1,sum(p.nums));
    icable = 0;

   % Initialize variables    
    V(:,:) = p.E_leak;                                                      % initialize first time point at E_leak
    n(:) = f.n_inf_Kd(p.E_leak,p.Vh_n_Kd(1),p.k_n_Kd(1));
    m(:) = f.m_inf_Na(p.E_leak,p.Vh_m_Na(1),p.k_m_Na(1));
    h(:) = f.h_inf_Na(p.E_leak,p.Vh_h_Na(1),p.k_h_Na(1));
      
    for j=1:length(p.t)-1
        k1v(1) = (i_ext(j)*p.A(1) - f.Ileak(V(1,j),p.gleak(1),p.E_leak) - f.I_syn(V(1,j),j,1,1) ...
            - f.I_Na(V(1,j),p.g_Na(1),m(1,j),h(1,j)) - f.I_Kd(V(1,j),p.g_Kd(1),n(1,j)))/ p.cm(1);
        if p.withcable
            for i = 2:p.branch_num+1                 % Calculate currents going into all the cables attached to the soma
                b = 2 + (i-2)*p.nums(i-1);
                icable = icable + (V(1,j)-V(b,j))/p.ra(i);
            end
        end
        k1v(1) = k1v(1) - icable/p.cm(1);
        k1m(1) = f.dm_dt(V(1,j),m(1,j),p.Vh_m_Na(1),p.k_m_Na(1),p.tau_m_Na(1));
        k1h(1) = f.dh_dt(V(1,j),h(1,j),p.Vh_h_Na(1),p.k_h_Na(1),p.tau_h_Na(1));
        k1n(1) = f.dn_dt(V(1,j),n(1,j),p.tau_lo_Kd(1),p.tau_hi_Kd(1),p.Vh_n_Kd(1),p.k_n_Kd(1));

        seg1 = 3;                                         % First segment of first branch (p.num_branch+1)
        if p.withcable                                    % if withcable=0 there is only a soma
            for br = 2:p.branch_num+1
                k = seg1-1;                               % 1st pos of each cable is a special case when p.num_branch>1
                k1v(k) = (- f.Ileak(V(k,j),p.gleak(br),p.E_leak) - f.I_syn(V(k,j),j,br,k-(seg1-3)) ... 
                        + (V(1,j)-V(k,j))/p.ra(br) - (V(k,j)-V(k+1,j))/p.ra(br)) / p.cm(br);
                for k = seg1 : seg1 + p.nums(br) - 3
                    k1v(k) = (- f.Ileak(V(k,j),p.gleak(br),p.E_leak) - f.I_syn(V(k,j),j,br,k-(seg1-3)) ... 
                        + (V(k-1,j)-V(k,j))/p.ra(br) - (V(k,j)-V(k+1,j))/p.ra(br)) / p.cm(br);
                end
            
                % Sealed end BC (seg1+p.nums(b)-2 is the index of the last segment in branch b)
                k = seg1 + p.nums(br) - 2;
                k1v(k) = (- f.Ileak(V(k,j),p.gleak(br),p.E_leak) - f.I_syn(V(k,j),j,br,k-(seg1-3)) ...
                    + (V(k-1,j)-V(k,j))/p.ra(br)) / p.cm(br);
    
                for k = seg1-1 : seg1 + p.nums(br) - 2
                    av(k) = V(k,j)+k1v(k)*p.dt;
                end
                seg1 = seg1 + p.nums(br);                                      % This is the first position of branch b > p.branch_num+1
            end
        end

        av(1) = V(1,j)+k1v(1)*p.dt;
        am = m(1,j)+k1m(1)*p.dt;
        an = n(1,j)+k1n(1)*p.dt;
        ah = h(1,j)+k1h(1)*p.dt;
        
        icable = 0;
        k2v(1) = (i_ext(j)*p.A(1) - f.Ileak(av(1),p.gleak(1), p.E_leak) - f.I_syn(av(1),j,1,1) ...
            - f.I_Na(av(1),p.g_Na(1),am,ah) - f.I_Kd(av(1),p.g_Kd(1),an)) / p.cm(1);
        %   - f.I_Kir(av(1),p.g_Kir(1),an)/p.cm(1);
        if p.withcable
            for i = 2:p.branch_num+1                                       % All the cables attached to the soma-1
                br = 2 + (i-2)*p.nums(i-1);
                icable = icable + (av(1)-av(br))/p.ra(i);
            end
        end
        k2v(1) = k2v(1) - icable/p.cm(1);
        k2m(1) = f.dm_dt(av(1),am,p.Vh_m_Na(1),p.k_m_Na(1),p.tau_m_Na(1));
        k2h(1) = f.dh_dt(av(1),ah,p.Vh_h_Na(1),p.k_h_Na(1),p.tau_h_Na(1));
        k2n(1) = f.dn_dt(av(1),an,p.tau_lo_Kd(1),p.tau_hi_Kd(1),p.Vh_n_Kd(1),p.k_n_Kd(1));

        seg1 = 3; 
        if p.withcable
           for br = 2:p.branch_num + 1
               k = seg1 - 1;                                            % 1st pos of each cable is a special case when p.num_branch>1
               k1v(k) = (- f.Ileak(V(k,j),p.gleak(br),p.E_leak) - f.I_syn(V(k,j),j,br,k-(seg1-3)) ... 
                        + (av(1)-av(k))/p.ra(br) - (av(k)-av(k+1))/p.ra(br)) / p.cm(br);
              
                for k = seg1 : seg1 + p.nums(br) - 3
                    k2v(k) = (- f.Ileak(av(k),p.gleak(br),p.E_leak) - f.I_syn(av(k),j,br,k-(seg1-3))...
                        + (av(k-1)-av(k))/p.ra(br) - (av(k)-av(k+1))/p.ra(br)) / p.cm(br);
                end

                %Sealed end BC
                k = seg1 + p.nums(br) - 2;  
                k2v(k) = (- f.Ileak(av(k),p.gleak(br),p.E_leak) - f.I_syn(av(k),j,br,k-(seg1-3)) ...
                    + (av(k-1)-av(k))/p.ra(br)) / p.cm(br);     
    
                for k = seg1-1:seg1 + p.nums(br) - 2 
                    V(k,j+1) = V(k,j) + (k1v(k) + k2v(k))/2*p.dt;
                end
                seg1 = seg1 + p.nums(br);                                      % This is the first position of branch b > p.branch_num+1
            end
        end

        V(1,j+1) = V(1,j) + (k1v(1) + k2v(1))/2*p.dt;
        m(1,j+1) = m(1,j) + (k1m + k2m)/2*p.dt;
        h(1,j+1) = h(1,j) + (k1h + k2h)/2*p.dt;
        n(1,j+1) = n(1,j) + (k1n + k2n)/2*p.dt;

        I(1,j) = f.Ileak(V(1,j),p.gleak(1), p.E_leak)*1000;                 % converted to pA
        I(2,j) = f.I_Kd(V(1,j), p.g_Kd(1), an)*1000;
        I(3,j) = f.I_Na(V(1,j), p.g_Na(1), am, ah)*1000;        
        % I(4,j) = f.I_syn(V(1,j), j, 2, p.nums(2)+1)*1000;
        % I(5,j) = f.I_syn(V(1,j), j, 3, p.nums(3)+1)*1000;
        I(6,j) = i_ext(j)*1000;

    end     % End time loop
end         % ----------------- end of function ------------------