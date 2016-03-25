%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to calculate hopping probabilities for TAT conduction
%
% COMPUTATION:
% Calculates Vo-Vo hopping probabilities from Mott hopping model
% Calculates Electrode-Vo and Vo-Electrode tunneling rates via WKB
% Tunneling from electrode to Vo accounts for integrating energies from Vo
% trap level to the conduction band too.
% Returns all tunneling rates in GHz (1/ns)
%
% INPUT VARIABLES:
% num_v:  the number of vacancy sites
% N,M:    the number of mesh points along y and x directions.
% vindex: N*M array. vindex(ii,jj) corresponds to a position at
%         (ii*dx, jj*dy). pos(ii,jj)==-1 if there is no vacancy at this site.
%         0<vindex(ii,jj)<num_v is the index of the vacancy.
%          Sites indexing scheme:
%          (N,1) (N,2) (N,3) ... (N,M)   y
%          ......                        ^
%          (2,1) (2,2) (2,3) ... (2,M)   |
%          (1,1) (1,2) (1,3) ... (1,M)   |---------->x
% vlist:  num_v*2 array. vlist(nn,1:2) is the (ii,jj) position index of the 
%         n'th vacancy.
% dx, dy: in units of nm. Mesh size along the x and y directions.
% a0:     in units of nm. Decaying length.
% P0:     in units of GHz (1/ns). Vibration frequency. 
% tunnel_rate0: in units of GHz (1/ns). Base tunneling rate.
% meff:   in units of m0. Effective mass.
% NtL/R:  Density of states in the L/R electrodes. in units of /eV.  
% EFL0 =  Fermi level of zero-biased left electrode
% EFR0 =  Fermi level of zero-biased right electrode
% V:      in units of V. N*M array. Electric potential inside the cell.
%         V(:,1) is the potential of the left electrode, while V(:,M) is
%         the potential of the right electrode.
% T:      in units of K. Device temperature. Scalar.
%
% Returned variables:
% P:      in units of GHz (1/ns). num_v*num_v array. P(nn,mm) is the
%         hopping rate from the nn'th site to the mm'th site. The diagonal
%         of P is always zero. 
% PiL:    in units of GHz (1/ns). num_v*1 array. PiL(nn) is the hopping-in
%         rate from the left electrode to the nn'th site.
% PiR:    in units of GHz (1/ns). num_v*1 array. PiR(nn) is the hopping-in
%         rate from the right electrode to the nn'th site.
% PoL:    in units of GHz (1/ns). num_v*1 array. PoL(nn) is the hopping-out
%         rate from the nn'th site to the left electrode.
% PoR:    in units of GHz (1/ns). num_v*1 array. PoR(nn) is the hopping-out
%         rate from the nn'th site to the right electrode.
% Note that the calculation of P is expensive in both time and memory. So
% this routine should be called only once before the Newton-Raphson
% interation starts.
% Ref:    ./doc/Current_continuity.docx.
% Author: Ximeng Guan
% Last modified: 1/15/2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [P, PiL, PoL, PiR, PoR] = hoppingProbability_uptocb(num_v, N, M, vindex, vlist, ...
    dx, dy, a0, P0, tunnel_rate0, meff, NtL, NtR, EFL0, EFR0, V, T)

% Initialize parameters
Ptol = 1e-4; % Ptol*electron = 1.6e-3 pA; Probability error tolerance
kB = 1.3806503e-23; % m^2kgs^{-2}K^{-1}
electron = 1.60217646e-19; % C
eV = 1.60217646e-19; % J
hbar = 1.05457148e-34; %m^2kg/s
m0 = 9.10938188e-31; % kg;
KBT = kB*T/electron; % eV

% Occupied and unoccupied state energy (base energy), energy zero is the Ec
% of a unbiased HfO2.
Eempty = -1.85; %-1.83; % -1.77eV, trap energy below conduction band, Vo charge 1+ state
Efilled = -1.95; %-1.97; % -2.01eV, trap energy below conduction band, Vo charge neutral state
%EFL0 = -1.9; %-2.3; % -1.9eV, Fermi level of the zero-biased left electrode
%EFR0 = -1.9; %-2.3; % eV, Fermi level of the zero-biased right electrode

% Position of the left and right electrode
% xL = dx; xR = M*dx; 

% Initialize return variables to 0's
P = zeros(num_v,num_v);
PiL = zeros(num_v,1); PoL = zeros(num_v,1);
PiR = zeros(num_v,1); PoR = zeros(num_v,1);

% Iterate over all vacancies
for nn = 1:num_v
    
    % Get coordinates of current vacancy
    ii_n = vlist(nn,1); jj_n = vlist(nn,2);
    
    for mm = 1:num_v
        if mm == nn
            P(mm,nn) = 0;
            continue;
        else
            % extract positions of both Vo's
            ii_m = vlist(mm,1); jj_m = vlist(mm,2);
            
            % Get Vo-Vo distance
            rn = [jj_n*dx; ii_n*dy];
            rm = [jj_m*dx; ii_m*dy];
            rnm = norm(rm - rn); % nm
            
            % Get potential difference between the Vo's
            Vmn = V(ii_m,jj_m) - V(ii_n,jj_n);
            
            % Get tunneling rate between Vo's
            % Tunneling rate is set to saturate at Vmn = 0.5 V
            if Vmn>=0.5
                P(nn,mm) = P0*exp(-rnm/a0+0.5/KBT);
            else
                P(nn,mm) = P0*exp(-rnm/a0+(Vmn-0.0)/KBT);
            end
        end
    end
%     xn = jj_n*dx;
%     VnL = V(ii_n, jj_n)-VL;
%     VnR = V(ii_n, jj_n)-VR;
%     if VnL>=0

    %%%%%%%% Calculate tunneling rate from left electrode into Vo %%%%%%%%%
    % modified by Ximeng, jan 15, 2012
    EFL = EFL0-V(ii_n,1);  % Fermi level of the left electrode
    NElist  = 32;          % Number of discrete energy levels above trap energy to list
    Elist = linspace(0, max(0, -Eempty), NElist); % Discrete energy levels from trap energy up to conduction band minimum
    fL = zeros(NElist, 1); % Values of something for each discrete energy level
    PiLlist = zeros(NElist, 1); % Hopping rates for each discrete energy level?
    
    % Iterate over each discrete energy level from Et to the CBM
    for jj=1:NElist
        % Absolute energy of current energy level above Et, accounting for
        % Et and potential differences, as a function of distance from
        % left electrode
                   % Scalar Scalar    1*jj_n vector  Scalar      
        Ebarrier = (-Eempty-Elist(jj)-V(ii_n,1:jj_n)+V(ii_n,jj_n)); % maximum barrier over jj will be Ec.
        
        % Nrefine = 20;
        % xx = linspace(0,jj_n*dx-dx,jj_n); xxx = linspace(0,jj_n*dx-dx,Nrefine);
        % Vfinegrid = spline(xx, V(ii_n,1:jj_n), xxx);
        % Ebarrier = (-Eempty-Elist(jj)-Vfinegrid+V(ii_n,jj_n));  % maximum barrier over jj will be -EFL+30KBT below Ec==0.
        
        % Reduce Ebarrier to vector of all entries >= 0
        Ebarrier = Ebarrier(Ebarrier>=0); 
        
        % Calculate k = sqrt(2mE)/hbar
        action_iL = 1/hbar*sqrt(2*meff*m0*Ebarrier*eV); % m
        
        % Sum up the number of positive barriers
        nbarrier = length(action_iL);
        
        % Create trapezoidal interpolation of all barriers
        TiL = exp(-2.0*sum((action_iL(2:nbarrier)+action_iL(1:nbarrier-1))*0.5)*dx*1e-9);
        
        % TiL = exp(-2.0*sum((action_iL(2:nbarrier)+action_iL(1:nbarrier-1))*0.5)*(jj_n-1)*dx/(Nrefine-1)*1e-9);
        
        % Get absolute energy for a discrete level above the current vacancy
        Enempty = Eempty+Elist(jj)-V(ii_n,jj_n); % eV
        
        % Get Fermi level of the left electrode
        EFL = EFL0-V(ii_n,1); 
        
        % Fermi distribution for tunneling from electrode fermi level this energy level
        fL(jj) = 1/(1+exp((Enempty-EFL)/KBT));
        
        % Tunneling rate from left electrode for this energy level
        % rate = base_rate*DOS_electrode*fermi_distribution*tranmission_probability
        PiLlist(jj) = NtL*tunnel_rate0/NElist*fL(jj)*TiL;
        
    end
    
    % PiL contains total tunneling rate = sum over all energy levels above
    % the trap level
    PiL(nn) = sum(PiLlist(1:NElist));
    
%     if nn==1
%         length(Ebarrier)
%         ss = sprintf('VL=%4.2fV, TiL(NElist)=%5.3e, fL(NElist)=%5.3e, PiL(1)=%5.3e',V(ii_n,1), TiL,fL(NElist),PiL(1));
%         disp(ss);
%         ss = sprintf('ToL(1)=%e, FL_out(1)=%e, PoL(1)=%e',ToL,FL_out,PoL(1));
%         disp(ss);
%     end
    % endofmodification Jan 15, 2012
    
    % Heights for all barriers between filled trap and left electrode
    Ebarrier = (-Efilled-V(ii_n,1:jj_n)+V(ii_n,jj_n));
    
    % Eliminate all barriers < 0
    Ebarrier = Ebarrier(Ebarrier>=0);
    
    % Get k = sqrt(2mE)/hbar
    action_oL = 1/hbar*sqrt(2*meff*m0*Ebarrier*eV); % m
    
    % Number of positive barriers
    nbarrier = length(action_oL);
    
    % Transmission probability from Vo to left electrode
    ToL = exp(-2.0*sum((action_oL(2:nbarrier)+action_oL(1:nbarrier-1))*0.5)*dx*1e-9);
    
    % Get absoluted energy for current vacancy
    Enfilled = Efilled-V(ii_n,jj_n); % eV
%     fL = 1/(1+exp((Enfilled-EFL)/KBT));

    % Fermi distribution for tunneling from Vo to electrode
    FL_out = KBT*log( 1 + exp((-EFL+Enfilled)/KBT) ); % eV, modified 5/11/2011.
%     PoL(nn) = NtL*tunnel_rate0*(1-fL)*ToL;

    % Total tunneling rate from filled Vo to left electrode
    PoL(nn) = NtL*tunnel_rate0*FL_out*ToL; % eV, modified 5/11/2011.
%     if nn==1
%         ss = sprintf('ToL(1)=%e, FL_out(1)=%e, PoL(1)=%e',ToL,FL_out,PoL(1));
%         disp(ss);
%     end
        
     %%%%%%% Calculate tunneling rate from right electrode into Vo %%%%%%%%
    % modified by Ximeng, jan 15, 2012
    EFR = EFR0-V(ii_n,M);       % Fermi level of the right electrode
    NElist  = 32;               % Sum over 32 discrete energy levels above the Vo trap level
    Elist = linspace(0, max(0, -Eempty), NElist); % net energy above the targeted vacancy level up to the conduction band min
    fR = zeros(NElist, 1);      % Fermi distribution at each energy level
    PiRlist = zeros(NElist, 1); % Tunneling rates at each energy level
    for jj=1:NElist
        %Ebarrier = (-Eempty-V(ii_n,jj_n:M)+V(ii_n,jj_n));
        
        % Absolute energy of current energy level above Et, accounting for
        % Et and potential differences, as a function of distance from
        % right electrode
                   % Scalar Scalar    1*jj_n vector  Scalar 
        Ebarrier = (-Eempty-Elist(jj)-V(ii_n,jj_n:M)+V(ii_n,jj_n)); % maximum barrier over jj will be Ec.
        Ebarrier = Ebarrier(Ebarrier>=0); % Delete all barriers<0
        
        % Get k = sqrt(2mE)/hbar
        action_iR = 1/hbar*sqrt(2*meff*m0*Ebarrier*eV); % m
        
        % Number of positive barriers
        nbarrier = length(action_iR);
        
        % Transmission probability from right electrode to Vo
        TiR = exp(-2.0*sum((action_iR(2:nbarrier)+action_iR(1:nbarrier-1))*0.5)*dx*1e-9);
        
        % Get absoluted energy for current energy level above current Vo
        Enempty = Eempty+Elist(jj)-V(ii_n,jj_n); % eV
        
        % Fermi level of the right electrode
        EFR = EFR0-V(ii_n,M);  
        
        % Fermi distribution for tunneling from right electrode to Vo
        fR(jj) = 1/(1+exp((Enempty-EFR)/KBT));
%         FR_in = EFR-Enempty + KBT*log( 1 + exp((Enempty-EFR)/KBT) ); % eV, modified 5/11/2011.

        % Tunneling rate from right electrode for this energy level
        % rate = base_rate*DOS_electrode*fermi_distribution*tranmission_probability
        PiRlist(jj) = NtR*tunnel_rate0/NElist*fR(jj)*TiR;
    end
    %PiR(nn) = NtR*tunnel_rate0*fR*TiR;
    %PiR(nn) = NtR*tunnel_rate0*FR_in*TiR;  % eV, modified 5/11/2011.
    
    % PiR contains total tunneling rate = sum over all energy levels above
    % the trap level
    PiR(nn) = sum(PiRlist(1:NElist));
    % endofmodification Jan 15, 2012
    
    %%%%%%% Calculate tunneling rate into right electrode %%%%%%%%
    % Heights for all barriers between filled trap and right electrode
    Ebarrier = (-Efilled-V(ii_n,jj_n:M)+V(ii_n,jj_n));
    Ebarrier = Ebarrier(Ebarrier>=0); % Only keep entries > 0
    
    % Get k = sqrt(2mE)/hbar
    action_oR = 1/hbar*sqrt(2*meff*m0*Ebarrier*eV); % m
    
    % Number of positive barriers
    nbarrier = length(action_oR);
    
    % Transmission probability from Vo to right electrode
    ToR = exp(-2.0*sum((action_oR(2:nbarrier)+action_oR(1:nbarrier-1))*0.5)*dx*1e-9);
    
    % Get absoluted energy for current Vo
    Enfilled = Efilled-V(ii_n,jj_n); % eV
    
%     fR = 1/(1+exp((Enfilled-EFR)/KBT));

    % Fermi distribution for tunneling from Vo to right electrode 
    FR_out = KBT*log( 1 + exp((-EFR+Enfilled)/KBT) ); % eV, modified 5/11/2011.
    
%     PoR(nn) = NtR*tunnel_rate0*(1-fR)*ToR;

    % Total tunneling rate from filled Vo to right electrode
    PoR(nn) = NtR*tunnel_rate0*FR_out*ToR;  % eV, modified 5/11/2011.
end
end

% filter out the small parts to improve condition
% idex = find(P<Ptol); P(idex) = 0.0;
% idex = find(PiL<Ptol); PiL(idex) = 0.0;
% idex = find(PoL<Ptol); PoL(idex) = 0.0;
% idex = find(PiR<Ptol); PiR(idex) = 0.0;
% idex = find(PoR<Ptol); PoR(idex) = 0.0;