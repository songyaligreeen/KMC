%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to create random initial configuration
%
% COMPUTATION:
% Creates two separate weak spots in the device. 
% The spots are ((N/2)-diameter) apart
% Each spot is num_diameter wide
% Each spot has 5 Vo's (this is a magic number that is fixed by the
% function)
% The 5 Vo's are located randomly within the spot
%
% INPUT VARIABLES:
% N = grid height = device width
% M = grid width = device height
% num_diameter = Width of weak spot (in mesh points)
%
% RETURN VARIABLES:
% vindex = NxM array; == -1 if there is no Vo at (i,j)
%                     == n if there is a Vo, where n is the index of the Vo
% vlist  = num_Vx2 array where num_V is the total number of Vo's
%          vlist(n,:) = [i j] coordinates of vacancy n
% num_v  = total number ov Vo's (scalar)
%
% TO DO:
% Add error checking - does not account for edge case where 2 Vo's are 
% randomly initialized at the same mesh point
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vindex,vlist,vmatrix,num_v]=creat_filament_random(N,M,num_diameter)

% Initialize variables
vindex = ones(N,M)*(-1);    % NxM array; == -1 where there is no vacancy
                            %            == n if there is a Vo, where n is the Vo's number
vmatrix = zeros(N,M);       % NxM array; == 0 where there is no vacancy
                            %            == 1 if there is a Vo
num_v=0;                    % Scalar; total number of vacancies

% gap=round(0/0.25); % 0 for LRS, 2 nm for HRS


% for ii=1:N
%     
% %         if ii>(N-num_diameter)/2 && ii<=(N+num_diameter)/2
% %         % add the vacancy
% %         vmatrix(ii,3)=1;
% %         num_v=num_v+1;
% %         vindex(ii,3) = num_v;
% %         vlist(num_v,1:2) = [ii,3]; 
% %         end
%    
%     for jj=2+gap:M-1;
%         if ii>(N-num_diameter)/2 && ii<=(N+num_diameter)/2
%         % add the vacancy
%         vmatrix(ii,jj)=1;
%         num_v=num_v+1;
%         vindex(ii,jj) = num_v;
%         vlist(num_v,1:2) = [ii,jj];
%         end
%     end
% end

% Initialize helper variables
random_num=5;             % Number of Vo's to insert at each weak spot
min=N/4-num_diameter/4;   % Creates min and max
max=N*3/4-num_diameter/4; % where max - min = N/2 and offset from center by -n_d/4

% Insert random_num Vo's in the first weak spot
for kk=1:random_num
   
   ii=ceil(num_diameter*rand()+min); % Random number from min to min+num_diameter
   jj=ceil((M-2)*rand())+1;          % Random number from 1 to M-1
   
   %  if (ii>(N/3-num_diameter/2) && ii<=(N/3+num_diameter/2)) || (ii>(N*2/3-num_diameter/2) && ii<=(N*2/3+num_diameter/2))
   vmatrix(ii,jj)=1;        % Insert Vo at (ii,jj)
   num_v=num_v+1;           % Increment number of Vo's
   vindex(ii,jj) = num_v;   % Insert Vo index in 2D array at Vo's location
   vlist(num_v,1:2) = [ii,jj]; % Insert Vo coordinates into list
   
end

% Insert random_num Vo's in the second weak spot
for kk=1:random_num
   
   ii=ceil(num_diameter*rand()+max); % Random number from max to max+num_diameter
   jj=ceil((M-2)*rand())+1;          % Random number from 1 to M-1
   
   %  if (ii>(N/3-num_diameter/2) && ii<=(N/3+num_diameter/2)) || (ii>(N*2/3-num_diameter/2) && ii<=(N*2/3+num_diameter/2))
   vmatrix(ii,jj)=1;        % Insert Vo at (ii,jj)
   num_v=num_v+1;           % Increment number of Vo's
   vindex(ii,jj) = num_v;   % Insert Vo index in 2D array at Vo's location
   vlist(num_v,1:2) = [ii,jj]; % Insert Vo coordinates into list

end

end
