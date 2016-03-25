%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code creates a simple RRAM cell in the same format as the RRAM KMC
% model. It uses the Comsol Livelink interface.
%
% It allows the user to change oxide dimensions, electrode thickness, and
% number of vacancies, and localization of vacancies.
%
% In COMSOL, Vo's and Oi's are approximated as square regions.
% It is necessary for Vo's to be approximated this way in order to allow
% them to form continuous filaments.
% Oi's are represented this way for consistency.
%
% Current device layout:
%
%      Electrode  Oxide  Electrode
%   N |---------|-------|---------|
%     |         |       |         |
%     |         |       |         |
%     |         |       |         |
%     |         |       |         |
%     |         |       |         |
%   0 |---------|-------|---------|
%    -tL        0       M        M+tR
%
%
%
% Eventual device layout:
%
% N+tair |---------------------------|
%        |                           |
%        |   Ambient (SiO2 or air)   |
%        |                           |
%        |                           |
%      N |---------|-------|---------|
%        |         |       |         |
%        |  pt     |       |   TiN   |
%        |Electrode| Oxide |Electrode|
%        |         |       |         |
%        |         |       |         |
%      0 |---------|-------|---------|
%        |                           |
%        |                           |
%        |   Ambient (SiO2 or air)   |
%        |                           |
%  -tair |---------|-------|---------|
%       -tL        0       M        M+tR
%
% Major variables:
% dev is a structure containing all major domains for ocide, electrodes,
% and defects
% dev.ox{1} is the oxide layer (note using curly brackets)
% dev.met is a cell array containing the metal electrodes
%   dev.met{1} = left electrode
%   dev.met{2} = right electrode
% dev.vo is an array of oxygen vacancies, each represented as a square
%   dev.vo{n} will access the nth vacancy (1 <= n <= numVo)
%
% 
% TO DO:
% (1) Figure out how to change a domain's material to one of the materials
% included in comsol, such as copper or aluminum
% (2) Write scripts generating user-defined materials HfO2, Hf, TiN, and Pt
% (3) Add SiO2 or O2 ambient to outside of device
% (4) Calculate Ohmic current
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%function [V,I,T]=COMSOL_heat(VL,N,M,vlist,T0,operation)

% Initialize workspace:
% clear all
clear model* geom* r* t* M N T
clear Vo_idx air_idx BE_idx TE_idx
close all
%addpath(genpath('../'))
home

%load('vlist.mat'); % CF not connected
%load('vlist_yali.mat'); % CF connected
num_v=length(vlist);
numVo = num_v;
% Initialize device geometry
N = 80;         % Grid height = device width (since electrodes go from left to right)
M = 42;         % Grid width  = device height
a = 0.25;       % Width/height of grid point (nm)
tL = 10;        % Left electrode thickness (nm)
tR = 10;        % Right electrode thickness (nm)

%numVo = 170; % Number of vacancies to introduce
%ww    = 5;    % Width of weak spot, in mesh points
              % All vacancies created within ww of middle of cell
useSeed = 1; % 1 to use a seed (holds random number generator constant)


% Get unscaled values of device dimensions
nm = 1e-9;      % Conversion from m to nm
dim.M = M*a*nm;
dim.N = N*a*nm;
dim.tL = tL*nm;
dim.tR = tR*nm;
dim.W = dim.M+dim.tL+dim.tR; % width of ambient SiO2 
dim.tair = dim.N/2;% length of ambient SiO2

VL = 3;
T0 = 300;
% Create oxide geometry
model = ModelUtil.create('Model');  % 5.94s (0.02s for future calls)

%set boundary value
model.param.set('T0', T0);
model.param.set('VL', VL);

%%%%%%%% set material parameters  Pt as top electrode %%%%%%%%
model.param.set('sigma_Pt', '8.9e6[S/m]');
model.param.set('Cv_Pt', '133[J/(kg*K)]');
model.param.set('dens_Pt', '21450[kg/m^3]');
model.param.set('kappa_Pt', '71.6[W/(m*K)]');

% Vo's
model.param.set('sigma_Vo', '1.55e3[S/m]');
model.param.set('Cv_Vo', '500[J/(kg*K)]');
model.param.set('dens_Vo', '9680[kg/m^3]');
model.param.set('kappa_Vo', '4[W/(m*K)]');
model.param.set('perm_Vo', '1');

% HfO2 dielectric
model.param.set('sigma_HfO2', '1e-5[S/m]');
model.param.set('Cv_HfO2', '500[J/(kg*K)]');
model.param.set('dens_HfO2', '9680[kg/m^3]');
model.param.set('kappa_HfO2', '1[W/(m*K)]');

% Cu as bottom electrode 
model.param.set('sigma_Cu', '58.1e6[S/m]');
model.param.set('Cv_Cu', '384[J/(kg*K)]');
model.param.set('dens_Cu', '8960[kg/m^3]');
model.param.set('kappa_Cu', '400[W/(m*K)]');

% SiO2 as ambient
model.param.set('sigma_SiO2', '1e-5[S/m]');
model.param.set('Cv_SiO2', '730[J/(kg*K)]');
model.param.set('dens_SiO2', '2200[kg/m^3]');
model.param.set('kappa_SiO2', '1.4[W/(m*K)]');
model.param.set('perm_SiO2', '4.2');

% Initialize 2D geometry named geom1
geom1 = model.geom.create('geom1',2); % 1.61s (0.04s for future calls)

%%%%%%%%%%%%%%%%%% CREATE MAJOR DOMAINS %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Rectangles for oxide(s) and electrodes %%%%%%%%%%%%
%%%%%%%%%% Eventually include SiO2 or air ambient %%%%%%%%%%%%

% Add a rectangle for oxide which stretches from (0,0) to (M,N)
dev.ox{1} = geom1.feature.create('r1','Rectangle'); % 0.17s (0.04 for future calls)
dev.ox{1}.set('size', [dim.M dim.N]);  % 0.02s (width and length, respectively)
dev.ox{1}.set('pos', [0 0]);   % 0.01s (position is bottom left corner of rectangle)

% Add a rectangle for the left electrode
% Stretches from (-tL,0) to (0,N)
dev.met{1} = geom1.feature.create('r2','Rectangle'); % 0.02s
dev.met{1}.set('size', [dim.tL dim.N]);  % 0.02s (width and length, respectively)
dev.met{1}.set('pos', [-dim.tL 0]);  % 0.01s (position is bottom left corner of rectangle)

% Add a rectangle for the left electrode
% Stretches from (M,0) to (M+tR,N)
dev.met{2} = geom1.feature.create('r3','Rectangle'); % 0.03s
dev.met{2}.set('size', [dim.tR dim.N]);  % 0.01s (width and length, respectively)
dev.met{2}.set('pos', [dim.M 0]);    % 0.01s (position is bottom left corner of rectangle)

% Add two rectangle SiO2 or O2 ambient to outside of devic
dev.amb{1} = geom1.feature.create('r4','Rectangle'); % 0.03s
dev.amb{1}.set('size', [dim.W dim.tair]);  % 0.01s (width and length, respectively)
dev.amb{1}.set('pos', [-dim.tL -dim.tair]);    % 0.01s (position is bottom left corner of rectangle)

dev.amb{2} = geom1.feature.create('r5','Rectangle'); % 0.03s
dev.amb{2}.set('size', [dim.W dim.tair]);  % 0.01s (width and length, respectively)
dev.amb{2}.set('pos', [-dim.tL dim.N]);    % 0.01s (position is bottom left corner of rectangle)

%%%%%%%%%%%%%%%%%%%%%% CREATE VACANCIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate vacancies randomly
% Vacancies occur anywhere from 0 to M with equal probability
% Vacancies evenly distributed within weak region width

% Apply seed if desired
% This makes outputs deterministic instead of random
if useSeed
    seed = 4321; % Arbitrary integer
    rng(seed);   % Sets all random number generators to same starting point
end


% Add all Vo's to model (0.025s per Vo)
dev.vo = {};

for i = 1:numVo
    mc = (vlist(i,2)-1)*a*nm; nc = (vlist(i,1)-1)*a*nm;
    dev.vo{i} = geom1.feature.create(['v',num2str(i)],'Square');
    dev.vo{i}.set('size', a*nm); % Squares are size of 1 grid point
    dev.vo{i}.set('pos', [mc nc]);
end



% Plot device geometry
 figure(1)
 mphgeom(model) % 0.16s (1.4s if many Vo's)

DN = geom1.getNDomains; % get total domain number 
Ox_domain = 1:DN; %initial total domains vector, after remove the Vo, TE, BE domain number,
%the left  number is HfO2 domain 
vo = zeros(numVo,2);
vo(:,1) = vlist(:,2);
vo(:,2) = vlist(:,1);
%%%%%%%%%%%%%%%%%%%%%% find  Vo domain index %%%%%%%%%%%%%%%%% 3.3209s
%%%%%%%%%%%%%%%%%%%%%% depend on the number of Vo
% use a rectangular to enclose one Vo, then get the domain number of this Vo

for i = 1:numVo
 vx1 = (vo(i,1)-1.2)*a*nm; vy1 = (vo(i,2)-1.2)*a*nm;
    vx2 = (vo(i,1)+0.2)*a*nm; vy2 = (vo(i,2)+0.2)*a*nm;
    Vo_idx(i) = mphselectbox(model,'geom1',[vx1 vx2; vy1 vy2], 'domain');
    Ox_domain = Ox_domain(find(Ox_domain~=Vo_idx(i)));%if the domain is Vo,remove the domain number
     
end

%%%%%%%%%%%%%%%%%%%%%% find the top electrode domain index%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%0.0230s
 
 TEx1 = -dim.tL*1.1; TEy1 = 0;
 TEx2 = 0; TEy2 = dim.N*1.1;
 TE_idx = mphselectbox(model,'geom1',[TEx1 TEx2; TEy1 TEy2], 'domain');
 
 Ox_domain = Ox_domain(find(Ox_domain~=TE_idx));%if the domain is TE,remove the domain number  
       
%%%%%%%%%%%%%%%%%%%%%% find the bottom electrode domain index %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% 0.0186s
 
 BEx1 = dim.M; BEy1 = 0;
 BEx2 = dim.M+dim.tR*1.1; BEy2 = dim.N*1.1;
 BE_idx = mphselectbox(model,'geom1',[BEx1 BEx2; BEy1 BEy2], 'domain');
 
 Ox_domain = Ox_domain(find(Ox_domain~=BE_idx));%if the domain is BE,remove the domain number  
 
 %%%%%%%%%%%%%%%%%%%%%% find the air/SiO2 domain index %%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%% 0.0415s

 for i = 1:2
 if i == 1
    Ax1 = -dim.tL*1.1; Ay1 = -dim.tair;
    Ax2 = dim.M+dim.tR*1.1; Ay2 = 0;
    air_idx(i) = mphselectbox(model,'geom1',[Ax1 Ax2; Ay1 Ay2], 'domain');
 else
    Ax3 = -dim.tL*1.1; Ay3 = dim.N;
    Ax4 = dim.M+dim.tR*1.1; Ay4 = dim.N+dim.tair;
    air_idx(i) = mphselectbox(model,'geom1',[Ax3 Ax4; Ay3 Ay4], 'domain');
 end
 Ox_domain = Ox_domain(find(Ox_domain~=air_idx(i)));%if the domain is air,
 %remove the domain number,the left  number is the HfO2 domain 
end

%%%%%%%%%%%%%%%%%% find the boundary conditions %%%%%%%%%%%%%%%%%%%%%% 0.1348s 
 TE_bc = mphselectbox(model,'geom1',[-dim.tL*1.1,-dim.tL;0,dim.N],'boundary');
 BE_bc = mphselectbox(model,'geom1',[dim.M+dim.tR,dim.M+dim.tR*1.1;0,dim.N],'boundary');
 

 Temp_bc1 = mphselectcoords(model,'geom1',[-dim.tL,-dim.tair],...
  'boundary','radius',dim.tair); 
 Temp_bc2 = mphselectcoords(model,'geom1',[-dim.tL,dim.N+dim.tair],...
  'boundary','radius',dim.tair); 
 Temp_bc3 = mphselectbox(model,'geom1',[dim.M+dim.tR,dim.M+dim.tR*1.1;-dim.tair*1.1,0],'boundary'); 
 Temp_bc4 = mphselectbox(model,'geom1',[dim.M+dim.tR,dim.M+dim.tR*1.1;dim.N,dim.N+dim.tair*1.1],'boundary'); 
 Temp_bc = [Temp_bc1 Temp_bc2 Temp_bc3 Temp_bc4];

%%%%%%%%%%%%%%%%%% ADD MATERIALS %%%%%%%%%%%%%%%%%%%%%%
% First, create the material label, then set the parameters. 
% The index of each domain is dynamic, depending on its position. 
% Generally, the domains are numberd from left to right. we use a rectangular 
% to enclose each Vo, then get the domain index give it to the material     

%creat and define parameters of metal electrode Pt as top electrode,

mat{1} = model.material.create('mat1','Common', 'mod1');
mat{1}.label('Pt');
mat{1}.propertyGroup('def').set('electricconductivity', 'sigma_Pt');
mat{1}.propertyGroup('def').set('heatcapacity', 'Cv_Pt');
mat{1}.propertyGroup('def').set('density', 'dens_Pt');
mat{1}.propertyGroup('def').set('thermalconductivity', 'kappa_Pt');
mat{1}.propertyGroup('def').set('relpermittivity', '1');
%choose the domain of Pt
mat{1}.selection.set(TE_idx); %0.1828s

%creat and define parameters of HfO2
mat{2} = model.material.create('mat2','Common', 'mod1');
mat{2}.label('HfO2');
mat{2}.propertyGroup('def').set('electricconductivity', 'sigma_HfO2');
mat{2}.propertyGroup('def').set('thermalconductivity', 'kappa_HfO2');
mat{2}.propertyGroup('def').set('density', 'dens_HfO2');
mat{2}.propertyGroup('def').set('heatcapacity', 'Cv_HfO2');
mat{2}.propertyGroup('def').set('relpermittivity', '25');

%choose the domain of HfO2
mat{2}.selection.set(Ox_domain);


%creat and define parameters of filament elements 
mat{3} = model.material.create('mat3', 'Common', 'mod1');
mat{3}.label('CF');
mat{3}.propertyGroup('def').set('electricconductivity', 'sigma_Vo');
mat{3}.propertyGroup('def').set('heatcapacity', 'Cv_Vo');
mat{3}.propertyGroup('def').set('density', 'dens_Vo');
mat{3}.propertyGroup('def').set('thermalconductivity', 'kappa_Vo');
mat{3}.propertyGroup('def').set('relpermittivity', 'perm_Vo');
%choose the domain of Vo, the domain is a seris number included in Vo_idx 
mat{3}.selection.set(Vo_idx);

%creat and define parameters of metal electrode Cu as bottom electrode
mat{4} = model.material.create('mat4', 'Common', 'mod1');
mat{4}.label('Cu');
mat{4}.propertyGroup('def').set('electricconductivity', 'sigma_Cu');
mat{4}.propertyGroup('def').set('heatcapacity', 'Cv_Cu');
mat{4}.propertyGroup('def').set('density', 'dens_Cu');
mat{4}.propertyGroup('def').set('thermalconductivity', 'kappa_Cu');
mat{4}.propertyGroup('def').set('relpermittivity', '1');

%choose the domain of Cu
mat{4}.selection.set(BE_idx);

%creat and define parameters of ambient SiO2
mat{5} = model.material.create('mat5', 'Common', 'mod1');
mat{5}.label('SiO2');
mat{5}.propertyGroup('def').set('electricconductivity', 'sigma_SiO2');
mat{5}.propertyGroup('def').set('thermalconductivity', 'kappa_SiO2');
mat{5}.propertyGroup('def').set('density', 'dens_SiO2');
mat{5}.propertyGroup('def').set('relpermittivity', 'perm_SiO2');
mat{5}.propertyGroup('def').set('heatcapacity', 'Cv_SiO2');
%choose the domain of SiO2.

mat{5}.selection.set(air_idx);


%%%%%%%%%%%%%%%%%% ADD PHYSICS %%%%%%%%%%%%%%%%%%%%%%
% here we choose the Joule heating, the physics including the Electric
% Current(ec) and Heat Transfer in Solids(ht) and the multiphysics 

ec = model.physics.create('ec', 'ConductiveMedia', 'geom1');
    ec.prop('d').set('d', '1[m]');
ht = model.physics.create('ht', 'HeatTransfer', 'geom1');
    ht.prop('PhysicalModelProperty').set('dz', '1[m]');

emh1 = model.multiphysics.create('emh1', 'ElectromagneticHeatSource', 'geom1', 2);
    emh1.selection.all;
    emh1.set('EMHeat_physics', 'ec');
    emh1.set('Heat_physics', 'ht');
bemh1 = model.multiphysics.create('bemh1', 'BoundaryElectromagneticHeatSource', 'geom1', 1);
    bemh1.selection.all;
    bemh1.set('EMHeatBoundary_physics', 'ec');
    bemh1.set('Heat_physics', 'ht');
tc1 = model.multiphysics.create('tc1', 'TemperatureCoupling', 'geom1');
    tc1.set('TemperatureSource_physics', 'ht');
    tc1.set('TemperatureDestination_physics', 'ec');

% above 1.2792s
% set the boundary conditions 

pot1 = ec.feature.create('pot1', 'ElectricPotential', 1);
    pot1.selection.set(TE_bc);   
    pot1.set('V0', 'VL');
    
   % if operation==1
   %    pot1.selection.set(BE_bc);   
   %    pot1.set('V0', 'VL'); 
   % else
   %    pot1.selection.set(TE_bc);   
   %    pot1.set('V0', 'VL'); 
   % end
   
gnd1 = ec.feature.create('gnd1', 'Ground', 1);
    gnd1.selection.set(BE_bc);  

   % if operation==1
   % gnd1.selection.set(TE_bc);
   % else
   % gnd1.selection.set(BE_bc);
   % end
ht.feature('init1').set('Tinit', 'T0'); % set initial T as T0    
temp1 = ht.feature.create('temp1', 'TemperatureBoundary', 1); %temperature 
    temp1.selection.set(Temp_bc);   
    temp1.set('T0', 'T0');


hs1 = ht.feature.create('hs1', 'HeatSource', 2);% heat source
% hs1.selection.all;
    hs1.set('Q_src', 'root.mod1.ec.Qh'); % choose the general source-total power dissipation density
    hs1.selection.set(Vo_idx);% choose the oxygen vacancies as heat source

    
    
    %%%%%%%%%%%%%%%%%% ADD MESH %%%%%%%%%%%%%%%%%%%%%% 0.2422s

    mesh1 = model.mesh.create('mesh1', 'geom1');
%fq1 is used to set the oxide mesh
fq1 = mesh1.create('fq1', 'FreeQuad');
    fq1.selection.geom('geom1', 2);
    fq1.selection.set([Vo_idx Ox_domain]);
size1 = fq1.create('size1', 'Size');
    size1.set('custom', 'on');
    size1.set('hmaxactive', 'on');
    size1.set('hminactive', 'on');
    size1.set('hmax', '0.25E-9'); %here the mesh we set as max 0.25nm,min0.15nm
    size1.set('hmin', '0.15E-9'); %because larger mesh size can cause problem if the Vo is very few 
    mesh1.run('fq1');
    % to disable fq1
    %fq1.active(false);
  
% other domains are using map, calculation time about 9~10s  
map1 = mesh1.create('map1', 'Map');
    map1.selection.remaining;
    mesh1.run('map1');
    % map1.active(false);% to disable map1

 % plot the mesh   
 figure (2)
    mphmesh(model); %0.1198s


%%%%%%%%%%%%%%%%%% ADD SDUTY %%%%%%%%%%%%%%%%%%%%%% 9.9746s
t = tic;
std1 = model.study.create('std1');
stat = std1.create('stat', 'Stationary'); % set the study type
    stat.activate('ec', true);
    stat.activate('ht', true);
    stat.set('notlistsolnum', 1);
    stat.set('notsolnum', '1');
    stat.set('listsolnum', 1);
    stat.set('solnum', '1');
% set the solver
sol1 = model.sol.create('sol1');
sol1.study('std1');

st1 = sol1.create('st1', 'StudyStep');
    st1.set('study', 'std1');
    st1.set('studystep', 'stat');
v1 = sol1.create('v1', 'Variables');
    v1.set('control', 'stat');
s1 = sol1.create('s1', 'Stationary');
    s1.create('seDef', 'Segregated');
    fc1 = s1.create('fc1', 'FullyCoupled');
        fc1.set('dtech', 'auto');
        fc1.set('initstep', 0.01);
        fc1.set('minstep', 1.0E-6);
        fc1.set('maxiter', 50);
    d1 = s1.create('d1', 'Direct');
        d1.set('linsolver', 'mumps');
        fc1.set('linsolver', 'd1');
        fc1.set('dtech', 'auto');
        fc1.set('initstep', 0.01);
        fc1.set('minstep', 1.0E-6);
        fc1.set('maxiter', 50);
s1.feature.remove('fcDef');
s1.feature.remove('seDef');

model.sol('sol1').attach('std1');
model.sol('sol1').runAll;
y10 = toc(t);
%%%%%%%%%%%%%%%%%% plot result %%%%%%%%%%%%%%%%%%%%%% 0.4259s

pg1 = model.result.create('pg1', 'PlotGroup2D');
    pg1.label('Electric Potential (ec)');
    pg1.set('oldanalysistype', 'noneavailable');
    pg1.set('frametype', 'spatial');
    pg1.set('data', 'dset1');
surf1 = pg1.feature.create('surf1', 'Surface');
	surf1.set('oldanalysistype', 'noneavailable');
    surf1.set('data', 'parent');
pg2 = model.result.create('pg2', 'PlotGroup2D');
    pg2.label('Temperature (ht)');
    pg2.set('oldanalysistype', 'noneavailable');
    pg2.set('data', 'dset1');
surf1 = pg2.feature.create('surf1', 'Surface');
    surf1.set('oldanalysistype', 'noneavailable');
    surf1.set('expr', 'T');
    surf1.set('colortable', 'ThermalLight');
    surf1.set('data', 'parent');  % choose the data you want to show

model.result('pg1').run;
model.result('pg2').run;

% save the model as 'RRAM_Vo'
mphsave(model,'KMC_RRAM_Vo'); % total 11.9355s

%find the node coordinate in KMC, the mesh node number is N*M, 
step = a*nm;
[x,y] = meshgrid(0:step:(dim.M-step), 0:step:(dim.N-step));
X = reshape(x,1,M*N);
Y = reshape(y,1,M*N);
data.coord = [X;Y];
%get the data for each point from COMSOL
data.T = mphinterp(model,'T','coord',data.coord);%0.16s

%give the data to array
T = zeros(N,M);
for ii = 1:length(data.T)
    yy = ceil(ii/N);
    if mod(ii,N)== 0
     xx = N;
    else
     xx = mod(ii,N);
    end
T(xx,yy)= data.T(ii);
end

figure(3)
    temperature5=fliplr(T);
    surf(temperature5')
    colorbar;
    axis('equal')
    view(0,90)
    
%pd_J = mpheval(model, 'ec.Jx');
Jx = mphinterp(model,'ec.Jx','coord',data.coord(:,1:80));%0.17s
I = sum(Jx)*dim.N*3*nm; %suppose the thickness is 5nm
data.V = mphinterp(model,'V','coord',data.coord);%0.17s

V = zeros(N,M);
for ii = 1:length(data.V)
    yy = ceil(ii/N);
    if mod(ii,N)== 0
     xx = N;
    else
     xx = mod(ii,N);
    end
V(xx,yy)= data.V(ii);
end

figure(4)
    V1=fliplr(V);
    surf(V1')
    colorbar;
    axis('equal')
    view(0,90)


