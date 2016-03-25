%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to update ion positions
%
% COMPUTATION:
%
% INPUT VARIABLES:
% vindex   = NxM array; == -1 if there is no Vo at (i,j)
%                       == n if there is a Vo, where n is the index of the Vo
% vlist   = num_Vx2 array where num_V is the total number of Vo's
%           vlist(n,:) = [i j] coordinates of vacancy n
% vmatrix = NxM array; vm(i,j) == 0 if no Vo; == 1 if there is a Vo
% num_v   = scalar; total number of Vo's
% N = grid height = device width
% M = grid width = device height
% tstep = s; length of time step; Probabilities are proportional to
%            tstep/t0, where t0 is 1/vibration_frequency
% cycle = 
% ion_list   = nionx2 array where nion = total number of ions
%              ion_list(n,:) == (i,j) coordinates of ion n
% ion_matrix =  NxM array 
%            == 0 if no Oi; == 1 if Oi
% num_ion    = Scalaer; total number of oxygen ions
% V = NxM array of potential at each grid point
% T = NxM array of temperature at each lattice point
% operation = 0 for form; 1 for reset; 2 for set
%
% RETURN VARIABLES:
% vindex = 
% vlist = 
% vmatrix = 
% num_v = 
% ion_list = nionx2 array where nion = total number of ions
%            ion_list(n,:) == (i,j) coordinates of ion n
% ion_matrix =  NxM array 
%            == 0 if no Oi; == 1 if Oi
% num_ion = 
% 
% TO DO:
% The "Create New Ion" section, around line 194, is biased toward upward
% motion. It should take a max over all the probabilities and use that one
% instead.
% Vectorize Vo deletion around line 276 and ion deletion at line 289
% Ion migration code around line 380 has bugs similar to deletion above
% Vo self-migration code is similarly biased
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vindex,vlist,vmatrix,num_v,ion_list,ion_matrix,num_ion]=DC_MC(vindex,vlist,vmatrix,num_v,N,M,tstep,cycle,ion_list,ion_matrix,num_ion,V,T,operation)

% num_v_before=num_v;
% num_ion_initial=num_ion;

lattice = 2.5; % A,  lattice constant in Angstroms
gamma=4; % field enhance factor, eA. should be 20

if operation==1 %% reset
    Ea=2.5;   % Vo generation; 4.4 for reset, 1.3 for set
    Em0=1.5;  % ion migration out of the GB; 2.2
    Em1=1.3;  % ion migration into the GB; 1.3
    Em2=1.5;  % Vo migration; 2.2 for reset, 4 for set
    Em3=1.5;  % interface_right
    Em4=6;    % interface_left (barrier into interface I think)
    Er=0.9;   % recombination, 0.8
    t0=1E-12; % s
elseif operation==0 %% forming
    Ea=2;     % Vo generation, 2.5 for forming, 1 for set
    Em0=0.25; % ion migration out of the GB, 1.5
    Em1=0.25; % ion migration in the GB, 1
    Em2=4;    % Vo migration, 3
    Em3=1.5;  % interface_right
    Em4=3;    % interface_left
    Er=6;     % recombination, 0.5 for reset, 3 for set
    t0=1E-12; % s
elseif operation==2 %% set
    Ea=1.1;   % Vo generation, 2.5 for forming, 1 for set
    Em0=0.8;  % ion migration out of the GB, 1.5
    Em1=0.3;  % ion migration in the GB, 1
    Em2=4;    % Vo migration, 3
    Em3=1.5;  % interface_right
    Em4=6;    % interface_left
    Er=6;     % recombination, 0.5 for reset, 3 for set
    t0=1E-12; % s
end

kT=T/297*0.0259;  

% cycle = number of times to repeat DC_MC()
for kk=1:cycle
        
    %%%%%%%%%% Vo generation
    test_prob=rand(N,M); % NxM array of initial probabilities for Vo generation, 1 for each mesh point
                         % Lower test_prob(ii,jj) --> higher chance of generating Vo
    for ii=1:N 
        for jj=round((M-1)/10)+1:M-1  % Iterate from 5 to 39 for case of M=40 (seems to skip an interfacial layer)
            
            %%%%%%%%%% Calculate field for one mesh point %%%%%%%%%%%%%%%%%
            % Approximates field as the difference in V between adjacent points
            % Takes the effective field for a point to be the maximum field
            % due to the potential difference with any of the 4 neighboring
            % points (up, down, left, right)
            if ii==1
                % Calculate field specifically for left edge of device
                F_local=max([abs(V(ii,jj+1)-V(ii,jj)),abs(V(ii,jj-1)-V(ii,jj)),abs(V(ii+1,jj)-V(ii,jj))])/lattice;
            elseif ii==N
                % Calculate field specifically for right edge of device
                F_local=max([abs(V(ii,jj+1)-V(ii,jj)),abs(V(ii,jj-1)-V(ii,jj)),abs(V(ii-1,jj)-V(ii,jj))])/lattice;
            else
                % Calculate field anywhere inside device between the out edges
                F_local=max([abs(V(ii,jj+1)-V(ii,jj)),abs(V(ii,jj-1)-V(ii,jj)),abs(V(ii+1,jj)-V(ii,jj)),abs(V(ii-1,jj)-V(ii,jj))])/lattice;
            end
            
            % Calculate Vo generation probability at mesh point
            % Prob_gen depends exponentially on local field F and 
            % activation energy Ea
            prob_gen(ii,jj)=(1-vmatrix(ii,jj))*exp(-(Ea-F_local*gamma)/kT(ii,jj))*tstep/t0;
            
            % Generate Vo's anywhere prob_gen is high enough relative to test_prob 
            if test_prob(ii,jj)<prob_gen(ii,jj)  %%% Vo generation
                vmatrix(ii,jj)=1;            % Insert Vo in vmatrix
                num_v=num_v+1;               % Increment number of Vo's
                vindex(ii,jj) = num_v;       % Add Vo to matrix of Vo's
                vlist(num_v,1:2) = [ii,jj];  % Add Vo to list of Vo's
                
                % Placeholder local migration barriers in each direction
                % Used 10 as an arbitrary large value so probability ~ 0
                % until barrier is later reduced
                Em_right=10;
                Em_left=10;
                Em_up=10;
                Em_down=10;
                
                % For the next 4 blocks:
                % Check each of the 4 directions around us to see if it
                % contains an Oi, a Vo, or neither
                % For each direction, verify we are also not going out of
                % the bounds of the device
                % For migration into a Vo, use the reduced Em1 for
                % migration into GB's
                % For migration not into a Vo, use the regular barrier Em0
                % Each migration barrier in each direction is also further
                % reduced by deltaV*gamma
                
                % Check to the right
                if jj<M && ion_matrix(ii,jj+1)==0 % Check if no Oi
                    if vindex(ii,jj+1)>0
                        % Oi=no and Vo=yes --> use Em1
                        Em_right=Em1-(V(ii,jj+1)-V(ii,jj))*gamma;
                    else
                        % Oi = no and Vo=no --> use Em0
                        Em_right=Em0-(V(ii,jj+1)-V(ii,jj))*gamma;
                    end
                end
                
                % Check to the left
                if jj>1 && ion_matrix(ii,jj-1)==0
                    if vindex(ii,jj-1)>0
                        % Oi=no and Vo=yes --> use Em1
                        Em_left=Em1-(V(ii,jj-1)-V(ii,jj))*gamma;
                    else
                        % Oi = no and Vo=no --> use Em0
                        Em_left=Em0-(V(ii,jj-1)-V(ii,jj))*gamma;
                    end
                end
                
                % Check above
                if ii>1 && ion_matrix(ii-1,jj)==0
                    if vindex(ii-1,jj)>0
                        % Oi=no and Vo=yes --> use Em1
                        Em_up=Em1-(V(ii-1,jj)-V(ii,jj))*gamma;
                    else
                        % Oi = no and Vo=no --> use Em0
                        Em_up=Em0-(V(ii-1,jj)-V(ii,jj))*gamma;
                    end
                end
                
                % Check below
                if ii<N && ion_matrix(ii+1,jj)==0
                    if vindex(ii+1,jj)>0
                        % Oi=no and Vo=yes --> use Em1
                        Em_down=Em1-(V(ii+1,jj)-V(ii,jj))*gamma;
                    else
                        % Oi = no and Vo=no --> use Em0
                        Em_down=Em0-(V(ii+1,jj)-V(ii,jj))*gamma;
                    end
                end
                
                % Calculate the migration rate in each direction
                rate_right=tstep/t0*exp(-Em_right/kT(ii,jj));
                rate_left=tstep/t0*exp(-Em_left/kT(ii,jj));
                rate_up=tstep/t0*exp(-Em_up/kT(ii,jj));
                rate_down=tstep/t0*exp(-Em_down/kT(ii,jj));
                
                % Random number for getting probability of ion motion
                test_num=rand;
                
                % Create new ion
                % NOTE: This section is biased and should be rewritten
                num_ion=num_ion+1; % New ion created, but must also see if it will move
                
                % Check each direction to see if and where the ion moves
                if test_num<rate_right
                    % Move ion rightwards if probability high enough
                    if jj<M
                        ion_list(num_ion,1:2)=[ii,jj+1];
                        ion_matrix(ii,jj+1)=1;
                    end
                elseif test_num<rate_right+rate_left
                    % Else move ion left if probability high enough
                    if jj>1
                        ion_list(num_ion,1:2)=[ii,jj-1];
                        ion_matrix(ii,jj-1)=1;
                    end
                elseif test_num<rate_right+rate_left+rate_up
                    % Else move ion up if probability high enough
                    if ii>1
                        ion_list(num_ion,1:2)=[ii-1,jj];
                        ion_matrix(ii-1,jj)=1;
                    end
                elseif test_num<rate_right+rate_left+rate_up+rate_down
                    % Else move ion down if probability high enough
                    if ii<N
                        ion_list(num_ion,1:2)=[ii+1,jj];
                        ion_matrix(ii+1,jj)=1;
                    end
                else
                    % If no probabilities were high enough, still create
                    % ion, but keep it in the current mesh site
                    ion_list(num_ion,1:2)=[ii,jj];
                    ion_matrix(ii,jj)=1;
                end
                
            end
        end
    end
    
    %%%%% Ion Motion %%%%%%%
    % Iterate over all ions
    % num_ion_initial=num_ion;
    for ss=1:num_ion % ion process
        
        % Get current ion
        xx=ion_list(ss,1);
        yy=ion_list(ss,2);
        
        if xx==0 || yy==0
            % Ion does not exist. Move to next ion
            % In current implementation, deleted ions are replaced with a
            % [0,0] entry instead of the row being deleted...see line
            % 289ish
            continue
        end
        
        % Vo-Oi recombination
        if vindex(xx,yy)>0 % If a Vo is here, calculate Vo-Oi recombinaton
            
            % F is the max field in any direction at this point
            % Obtained by getting potential difference with surrounding points
            if xx==1
                % Field for left edge of device
                F_local=max([abs(V(xx,yy+1)-V(xx,yy)),abs(V(xx,yy-1)-V(xx,yy)),abs(V(xx+1,yy)-V(xx,yy))])/lattice;
            elseif xx==N
                % Field for right edge of device
                F_local=max([abs(V(xx,yy+1)-V(xx,yy)),abs(V(xx,yy-1)-V(xx,yy)),abs(V(xx-1,yy)-V(xx,yy))])/lattice;
            elseif yy==1
                % Field for bottom of device
                F_local=max([abs(V(xx,yy+1)-V(xx,yy)),abs(V(xx+1,yy)-V(xx,yy)),abs(V(xx-1,yy)-V(xx,yy))])/lattice;
            elseif yy==M
                % Field for top of device
                F_local=max([abs(V(xx,yy-1)-V(xx,yy)),abs(V(xx+1,yy)-V(xx,yy)),abs(V(xx-1,yy)-V(xx,yy))])/lattice;
            else
                % Field in interior of device, where all 4 neighbor sites are checked
                F_local=max([abs(V(xx,yy+1)-V(xx,yy)),abs(V(xx,yy-1)-V(xx,yy)),abs(V(xx+1,yy)-V(xx,yy)),abs(V(xx-1,yy)-V(xx,yy))])/lattice;
            end
            
            % Recombination rate rises exponentially with field
            rate_recomb=tstep/t0*exp(-(Er-F_local*gamma)/kT(xx,yy));
            
            if rand<rate_recomb %%% if recombine
                % Delete the VO
                vmatrix(xx,yy)=0;
                num_delete=vindex(xx,yy) ;
                for ll=num_delete:num_v-1 % Delete this row from vlist and entry from vindex (very inefficient - can vectorize)
                    xx_n = vlist(ll+1,1);
                    yy_n = vlist(ll+1,2);
                    vindex(xx_n,yy_n) = vindex(xx_n,yy_n)-1;
                    vlist(ll,1:2) = vlist(ll+1,1:2);
                end
                vlist(num_v,1:2)=[0,0];
                vindex(xx,yy)=-1;
                num_v=num_v-1;
                
                % Delete the oxygen ion (very efficient - can be
                % vectorized)
                ion_matrix(xx,yy)=0;
                for pp=ss:num_ion-1
                    ion_list(pp,1:2) = ion_list(pp+1,1:2);
                end
                ion_list(num_ion,1:2)=[0,0];
                num_ion=num_ion-1;
                
            else % if not recombine, ion moving
                
                % Next blocks of code calculate ion migration probabilities
                % in all directions: left, right, up, and down
                
                % Placeholder directional migration barriers
                % Initially set as just very high numbers so Em ~ Inf in
                % each direction
                Em_right=10;
                Em_left=10;
                Em_up=10;
                Em_down=10;

                % Check ion motion in all directions
                if yy==1 % Ion is at left interface
                    % Need interfacial barriers Em3 and Em4
                    % Em4 is higher and is barrier for migrating into
                    % electrode
                    Em_left=Em4-(V(xx,yy)-V(xx,yy+1))*gamma;
                    if ion_matrix(xx,yy+1)==0
                        % If no ion in the way, use normal barrier Em3 for
                        % migration to the right back into device
                        Em_right=Em3-(V(xx,yy+1)-V(xx,yy))*gamma;
                    end
                elseif yy<M && ion_matrix(xx,yy+1)==0  
                    % Ion not at left interface and no Oi to the right
                    if vindex(xx,yy+1)>0 
                        % There is a Vo to the right
                        % Use reduced Em1 for migration into GB
                        Em_right=Em1-(V(xx,yy+1)-V(xx,yy))*gamma;
                    else
                        % No Vo to the right. Use normal Em0
                        Em_right=Em0-(V(xx,yy+1)-V(xx,yy))*gamma;
                    end
                end
                
                % Ion not at left interface and no Oi in spot to left
                if yy>1 && ion_matrix(xx,yy-1)==0
                    if vindex(xx,yy-1)>0
                        % Vo in next spot --> Use reduced Em1
                        Em_left=Em1-(V(xx,yy-1)-V(xx,yy))*gamma;
                    else
                        % No Vo in next spot --> Use normal Em0
                        Em_left=Em0-(V(xx,yy-1)-V(xx,yy))*gamma;
                    end
                end
                
                % Ion not at top and no Oi above
                if xx>1 && ion_matrix(xx-1,yy)==0
                    if vindex(xx-1,yy)>0
                        % Vo in next spot --> Use reduced Em1
                        Em_up=Em1-(V(xx-1,yy)-V(xx,yy))*gamma;
                    else
                        % No Vo in next spot --> Use normal Em0
                        Em_up=Em0-(V(xx-1,yy)-V(xx,yy))*gamma;
                    end
                end
                
                % Ion not at bottom and no Oi in spot below
                if xx<N && ion_matrix(xx+1,yy)==0
                    if vindex(xx+1,yy)>0
                        % Vo in next spot --> Use reduced Em1
                        Em_down=Em1-(V(xx+1,yy)-V(xx,yy))*gamma;
                    else
                        % No Vo in next spot --> Use normal Em0
                        Em_down=Em0-(V(xx+1,yy)-V(xx,yy))*gamma;
                    end
                end
                
                % Migration probability in each direction 
                % (called rate here but is treated like a probability)
                rate_right=tstep/t0*exp(-Em_right/kT(xx,yy));
                rate_left=tstep/t0*exp(-Em_left/kT(xx,yy));
                rate_up=tstep/t0*exp(-Em_up/kT(xx,yy));
                rate_down=tstep/t0*exp(-Em_down/kT(xx,yy));
                
                % Random chance of migration
                % Lower number = higher migration probability
                test_num=rand;
                
                % Determine if/where the ion moves
                % NOTE: This code is biased for rightward motion
                % NOTE: Code also fails for ions at edge (they will get
                % stuck...need to combine the two if statements in each
                % block)
                if test_num<rate_right
                    if yy<M % Verify we aren't at the right edge
                        % Update ion coordinates
                        ion_list(ss,1:2)=[xx,yy+1];
                        ion_matrix(xx,yy+1)=1;
                        ion_matrix(xx,yy)=0;
                    end
                elseif test_num<rate_right+rate_left
                    if yy>1 % Verify we aren't at left edge
                        % Update ion coordinates leftward
                        ion_list(ss,1:2)=[xx,yy-1];
                        ion_matrix(xx,yy-1)=1;
                        ion_matrix(xx,yy)=0;
                    elseif yy==1
                        % We are at left interface
                        % Ion migrates into electrode and disappears
                        for pp=ss:num_ion-1  %% delete ion
                            ion_list(pp,1:2) = ion_list(pp+1,1:2);
                        end
                        ion_list(num_ion,1:2)=[0,0];
                        num_ion=num_ion-1;
                        ion_matrix(xx,yy)=0;
                    end
                elseif test_num<rate_right+rate_left+rate_up
                    if xx>1 % Verify we aren't at the top
                        % Update ion coordinates upward
                        ion_list(ss,1:2)=[xx-1,yy];
                        ion_matrix(xx-1,yy)=1;
                        ion_matrix(xx,yy)=0;
                    end
                elseif test_num<rate_right+rate_left+rate_up+rate_down
                    if xx<N % Verify we aren't at the bottom
                        % Update ion coordinates downward
                        ion_list(ss,1:2)=[xx+1,yy];
                        ion_matrix(xx+1,yy)=1;
                        ion_matrix(xx,yy)=0;
                    end
                else
                    % Ion doesn't move, but it is still a new ion
                    % Add its entry to the ion list and matrix
                    ion_list(ss,1:2)=[xx,yy];
                    ion_matrix(xx,yy)=1;
                end
                
            end
        else % No Vo is here, but there is an ion here.
             % Calculate ion hopping isntead of Oi-Vo recombination
            
            % Set up placeholder directional ion migration barriers
            % Just using arbitrary large numbers
            Em_right=10;
            Em_left=10;
            Em_up=10;
            Em_down=10;
            
            if yy==1 %%%%%%%%%%%%%%%%%%%%%%%%% interface
                % At left interface, ion can migrate into interface but
                % faces a larger barrier Em4
                Em_left=Em4-(V(xx,yy)-V(xx,yy+1))*gamma;
                if ion_matrix(xx,yy+1)==0
                    % Migration away from interface uses barrier Em3
                    % (typically same barrier as in bulk)
                    Em_right=Em3-(V(xx,yy+1)-V(xx,yy))*gamma;
                end
            elseif yy<M && ion_matrix(xx,yy+1)==0
                % Not at interface and spot to right is free
                if vindex(xx,yy+1)>0
                    % Vo in spot; use lower barrier Em1 for migration into GB's 
                    Em_right=Em1-(V(xx,yy+1)-V(xx,yy))*gamma;
                else
                    % No Vo, use normal barrier Em0
                    Em_right=Em0-(V(xx,yy+1)-V(xx,yy))*gamma;
                end
            end
            
            if yy>1 && ion_matrix(xx,yy-1)==0
                % No Oi in spot to left
                if vindex(xx,yy-1)>0
                    % Vo in spot; use lower barrier Em1 for migration into GB's
                    Em_left=Em1-(V(xx,yy-1)-V(xx,yy))*gamma;
                else
                    % No Vo, use normal barrier Em0
                    Em_left=Em0-(V(xx,yy-1)-V(xx,yy))*gamma;
                end
            end
            
            if xx>1 && ion_matrix(xx-1,yy)==0
                % No Oi in spot above
                if vindex(xx-1,yy)>0
                    % Vo in spot; use lower barrier Em1 for migration into GB's
                    Em_up=Em1-(V(xx-1,yy)-V(xx,yy))*gamma;
                else
                    % No Vo, use normal barrier Em0
                    Em_up=Em0-(V(xx-1,yy)-V(xx,yy))*gamma;
                end
            end
            
            if xx<N && ion_matrix(xx+1,yy)==0
                % No Oi in spot below
                if vindex(xx+1,yy)>0
                    % Vo in spot; use lower barrier Em1 for migration into GB's
                    Em_down=Em1-(V(xx+1,yy)-V(xx,yy))*gamma;
                else
                    % No Vo, use normal barrier Em0
                    Em_down=Em0-(V(xx+1,yy)-V(xx,yy))*gamma;
                end
            end
            
            % Calculate migration probabilities
            % (Called rates, but treated like probabilities)
            rate_right=tstep/t0*exp(-Em_right/kT(xx,yy));
            rate_left=tstep/t0*exp(-Em_left/kT(xx,yy));
            rate_up=tstep/t0*exp(-Em_up/kT(xx,yy));
            rate_down=tstep/t0*exp(-Em_down/kT(xx,yy));
            
            % Random odds of migration
            test_num=rand;
            
            if test_num<rate_right
                % Move right if not at right interface
                if yy<M 
                    ion_list(ss,1:2)=[xx,yy+1];
                    ion_matrix(xx,yy+1)=1;
                    ion_matrix(xx,yy)=0;
                end
            elseif test_num<rate_right+rate_left
                % Move left
                if yy>1
                    % Move ion if not at left interface
                    ion_list(ss,1:2)=[xx,yy-1];
                    ion_matrix(xx,yy-1)=1;
                    ion_matrix(xx,yy)=0;
                elseif yy==1
                    % If at left interface, ion is absorbed into electrode
                    % and disappears
                    for pp=ss:num_ion-1  %% delete ion
                        ion_list(pp,1:2) = ion_list(pp+1,1:2);
                    end
                    ion_list(num_ion,1:2)=[0,0];
                    num_ion=num_ion-1;
                    ion_matrix(xx,yy)=0;
                end
            elseif test_num<rate_right+rate_left+rate_up
                % Move up if not at top of device
                if xx>1
                    ion_list(ss,1:2)=[xx-1,yy];
                    ion_matrix(xx-1,yy)=1;
                    ion_matrix(xx,yy)=0;
                end
            elseif test_num<rate_right+rate_left+rate_up+rate_down
                % Move down if not at bottom of device
                if xx<N
                    ion_list(ss,1:2)=[xx+1,yy];
                    ion_matrix(xx+1,yy)=1;
                    ion_matrix(xx,yy)=0;
                end
            else
                % No motion happens, but this is still a new ion, so
                % add the ion to the list and matrix
                ion_list(ss,1:2)=[xx,yy];
                ion_matrix(xx,yy)=1;
            end
            
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%
        
        % Vo self-migration (substitutional diffusion)
        % Iterate over all Vo's
        for ff=1:num_v % calculate the Vo self-migration
            
            % Get Vo coordinates
            xx=vlist(ff,1);
            yy=vlist(ff,2);
            if xx==0 || yy==0
                % Vo was deleted at some point. Move to next one
                continue
            end
            
            % Placeholder directional Vo migration barriers
            % These are just arbitrary high numbers that will be updated
            Em_right=10;
            Em_left=10;
            Em_up=10;
            Em_down=10;
            
            if yy<M-1
                % If not at right edge and spot to right has no Vo
                if vindex(xx,yy+1)<0
                    Em_right=Em2+(V(xx,yy+1)-V(xx,yy))*gamma;
                end
            end
            
            if yy>round((M-1)/10)+1 %%%%%%%%%%%%%%%%%%%%%%%%% 2
                % If not in the left interfacial region where Oi's gather
                % And if there is no Vo in the spot to the left
                if vindex(xx,yy-1)<0
                    Em_left=Em2+(V(xx,yy-1)-V(xx,yy))*gamma;
                end
            end
            
            if xx>1
                % If not at top of device and no Vo above
                if vindex(xx-1,yy)<0
                    Em_up=Em2+(V(xx-1,yy)-V(xx,yy))*gamma;
                end
            end
            
            if xx<N
                % If not at bottom of device and no Vo below
                if vindex(xx+1,yy)<0
                    Em_down=Em2+(V(xx+1,yy)-V(xx,yy))*gamma;
                end
            end
            
            % Extract migration probabilities
            % (Called rates here, but treated like probabiltiies)
            rate_right=tstep/t0*exp(-Em_right/kT(xx,yy));
            rate_left=tstep/t0*exp(-Em_left/kT(xx,yy));
            rate_up=tstep/t0*exp(-Em_up/kT(xx,yy));
            rate_down=tstep/t0*exp(-Em_down/kT(xx,yy));
            
            % Generate random migration probability
            test_num=rand;
            
            % NOTE: This code is biased rightward
            if  test_num<rate_right
                % Migrate right if not at right interface
                if yy<M-1
                    vlist(ff,1:2)=[xx,yy+1];
                    vindex(xx,yy)=0;
                    vindex(xx,yy+1)=ff;
                    vmatrix(xx,yy)=0;
                    vmatrix(xx,yy+1)=1;
                end
            elseif test_num<rate_right+rate_left
                % Migrate left if not at left interface
                if yy>2
                    vlist(ff,1:2)=[xx,yy-1];
                    vindex(xx,yy)=0;
                    vindex(xx,yy-1)=ff;
                    vmatrix(xx,yy)=0;
                    vmatrix(xx,yy-1)=1;
                end
            elseif test_num<rate_right+rate_left+rate_up
                % Migrate up if not at top edge
                if xx>1
                    vlist(ff,1:2)=[xx-1,yy];
                    vindex(xx,yy)=0;
                    vindex(xx-1,yy)=ff;
                    vmatrix(xx,yy)=0;
                    vmatrix(xx-1,yy)=1;
                end
            elseif test_num<rate_right+rate_left+rate_up+rate_down
                % Migrate down if not at bottom edge
                if xx<N
                    vlist(ff,1:2)=[xx+1,yy];
                    vindex(xx,yy)=0;
                    vindex(xx+1,yy)=ff;
                    vmatrix(xx,yy)=0;
                    vmatrix(xx+1,yy)=1;
                end
            % No need to do an else here, since Vo is not new and is
            % already in the vindex
            end
        end
    end
    
end

end



