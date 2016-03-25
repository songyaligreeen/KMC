%clear all
close all
clc

%%%%%%%%% INITIAL SETUP     %%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Runtime = 0.01 ms %%%%%%%%%%%%%%%%%%%%
N=80;           % Grid height = device width (since electrodes go from left to right)
M=42;           % Grid width  = device height
diameter=4;     % weak spot region width (in mesh points)
tstep=1E-4;     % s; Size of timestep
% tstep0=tstep/10; % reduce tstep
cycle=1; 
operation=0; %0 forming, 1 reset, 2 set;

Vstop=8;
Vstep=0.5; %0.02
nstep=round(abs(Vstop/Vstep));

V_sat=0.4;
tao=1E-7;
time_step=tao/40;

if operation == 1
    Icomp=1E-4;% for RESET
else
Icomp=1E-6; % for FORM and SET
end


if operation==0
[vindex,vlist,vmatrix,num_v]=creat_filament_random(N,M,diameter);    
[ion_list,ion_matrix,num_ion]=creat_ion_initial(N,M);
% load('initial.mat','vindex','vlist','vmatrix','num_v','ion_list','ion_matrix','num_ion')
elseif operation==1
load('LRS_new3.mat','vindex','vlist','vmatrix','num_v','ion_list','ion_matrix','num_ion')
elseif operation==2
load('HRS.mat','vindex','vlist','vmatrix','num_v','ion_list','ion_matrix','num_ion')
end

num_v_before=num_v;
num_ion_before=num_ion;
%%%%%%%%%%%%%%% electron current parameters for TAT
dx = 0.25; % nm
dy = 0.25; % nm
a0 = 0.33; %nm
P0 = 1e3; % GHz 
tunnel_rate0 = 1e4; % GHz edit here! 
meff = 0.1; % m0
NtL = 1; NtR = 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% parameters for FN
Lz = 0.25; % nm
me = 0.10; mh = 0.10; mt = 0.10; % m0
EFL0 = -1.9; % -1.9eV, Fermi level of the zero-biased left electrode
EFR0 = -1.9; % eV, Fermi level of the zero-biased right electrode
Eg = 5.8; % eV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%% plot before reset
figure(1)
vmatrix1=fliplr(vmatrix);
cla
hold on
p=bar3(vmatrix1',1,'r');
set(p,'facecolor',[1,0.412,0.706])
bar3(zeros(N,M)',1,'w');
axis('equal')
axis([1 N 0.5 M+0.5 0 1.001]); 

% ion_matrix = zeros(N,M);
% for kk=1:num_ion
%     if ion_list(kk,1)>0 && ion_list(kk,2)>0 
%        ion_matrix(ion_list(kk,1),ion_list(kk,2))=1;
%     end
% end

figure(1)
ion_matrix1=fliplr(ion_matrix);
hold on
bar3(ion_matrix1',1,'b');
bar3(zeros(N,M)',1,'w');
axis('equal')
axis([1 N 0.5 M+0.5 0 1.001]); 
%%%%%%%%%%%%%%

for ii=1:nstep
T0=273+25;
V_apply(ii)=Vstop/nstep*ii;
% V_apply(ii)=Vstop;
% time(ii)=ii*tstep;

%%%%%%%%%%%%% calculate the current
VL = V_apply(ii); % V 
VR = 0; % V

if ii==1
    V_initial = ones(N,1)*linspace(VL,VR,M); % linear potential 
else
    V_initial=V;
end

%%%%%%%% Ohmic current %%%%% need to change to COMSOL
t = tic;
[V,I_ohmic(ii),temperature]=COMSOL_heat(VL,N,M,vlist,T0,operation);
%[V,I_ohmic(ii),temperature,R_lateral,R_vertical]=KVL2_old(VL,VR,N,M,V_initial,vmatrix,operation);
x1(ii) = toc(t);
%%%%%%%%%%%%% TAT current
t = tic;
f = zeros(num_v,1);
trial = 0;
finit = ones(num_v,1)*0.5;
[I_TAT(ii), f, flag, restart_num] = tatCurrent_uptocb (finit, num_v, N, M, vindex, vlist, ...
    dx, dy, a0, P0, tunnel_rate0, meff, NtL, NtR, EFL0, EFR0, V, mean(temperature(:,1)));   
x2(ii) = toc(t);
%%%%%%%%%%%% FN current

EFR = EFR0; EFL = EFL0 - VL; 
TL = mean(temperature(:,1)); TR = mean(temperature(:,M)); %K
I_FN(ii) = FNcurrent (N, ...
dx, dy, Lz, me, mh, mt, EFL, EFR, Eg, V, TL, TR);    

I_total(ii)=I_TAT(ii)+I_ohmic(ii)+I_FN(ii); 

%%%%%%%%%%%%% overshoot
% while I_ohmic(ii)>=Icomp 

%  [vindex,vlist,vmatrix,num_v,ion_list,num_ion]=DC_MC(vindex,vlist,vmatrix,num_v,N,M,tstep0,cycle,ion_list,num_ion,V,temperature,operation);   
%  [V,I_ohmic(ii),temperature,R_lateral,R_vertical]=KVL2_old(VL,VR,N,M,V_initial,vmatrix,operation);
 
 if I_ohmic(ii)>=Icomp 
     
    figure(4)
    vmatrix4=fliplr(vmatrix);
    cla
    hold on
    p=bar3(vmatrix4',1,'r');
    set(p,'facecolor',[1,0.412,0.706])
    bar3(zeros(N,M)',1,'w');
    axis('equal')
    axis([1 N 0.5 M+0.5 0 1.001]); 

%     ion_matrix = zeros(N,M);
%     for kk=1:num_ion
%         if ion_list(kk,1)>0 && ion_list(kk,2)>0 
%            ion_matrix(ion_list(kk,1),ion_list(kk,2))=1;
%         end
%     end

    figure(4)
    ion_matrix4=fliplr(ion_matrix);
    hold on
    bar3(ion_matrix4',1,'b');
    bar3(zeros(N,M)',1,'w');
    axis('equal')
    axis([1 N 0.5 M+0.5 0 1.001]); 

    figure(6)
    temperature6=fliplr(temperature);
    surf(temperature6')
    colorbar;
    axis('equal')
    view(0,90)
    %%%%%%%%%%%%% overshoot period  %%%% KVL2_old need to change to COMSOL
    for gg=1:200
        V_drop(gg)=V_sat+(VL-V_sat)*exp(-time_step*gg/tao);
        V_initial=V;
        [V,I_transient(gg),temperature]=COMSOL_heat(VL,N,M,vlist,T0,operation);
        %[V,I_transient(gg),temperature,R_lateral,R_vertical]=KVL2_old(V_drop(gg),VR,N,M,V_initial,vmatrix,operation);
        [vindex,vlist,vmatrix,num_v,ion_list,ion_matrix,num_ion]=DC_MC(vindex,vlist,vmatrix,num_v,N,M,time_step,cycle,ion_list,ion_matrix,num_ion,V,temperature,operation);        
        if num_v>N*M/3
            disp('big bang');
            break
        end
    end

    figure(5)
    plot(I_transient);
    break
    end
% end
% 
% if I_ohmic(ii)>=Icomp
%     break
% end
% %%%%%%%%%%% ion process

[vindex,vlist,vmatrix,num_v,ion_list,ion_matrix,num_ion]=DC_MC(vindex,vlist,vmatrix,num_v,N,M,tstep,cycle,ion_list,ion_matrix,num_ion,V,temperature,operation);
fprintf(1,'%5d/%5d finished.',ii,nstep);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if operation==0 || operation==2
index=100;
for ii=1:nstep

V_back(ii)=ii*Vstop/nstep;

if ii>=index
    I_ohmic_back(ii)=Icomp;
    continue
end

VL = V_back(ii); % V 
VR = 0; % V
%%%%%% need to change to COMSOL
if ii==1
    V_initial = ones(N,1)*linspace(VL,VR,M); % linear potential 
    [V,I_ohmic_back(1),temperature]=COMSOL_heat(VL,N,M,vlist,T0,operation);
    %[V,I_ohmic_back(1),temperature,R_lateral,R_vertical]=KVL2_old(VL,VR,N,M,V_initial,vmatrix,operation); 
    index=ceil(Icomp/I_ohmic_back(1));
else
    V_initial=V;
end

%%%%%%%% Ohmic current 0.1061s
[V,I_ohmic_back(ii),temperature]=COMSOL_heat(VL,N,M,vlist,T0,operation);
%[V,I_ohmic_back(ii),temperature,R_lateral,R_vertical]=KVL2_old(VL,VR,N,M,V_initial,vmatrix,operation); 

end

else
for ii=nstep:-1:1
%%%%%%%%%%%%% calculate the current
V_back(ii)=ii*Vstop/nstep;

VL = V_back(ii); % V 
VR = 0; % V

if ii==1
    V_initial = ones(N,1)*linspace(VL,VR,M); % linear potential 
else
    V_initial=V;
end

%%%%%%%% Ohmic current
t = tic;
[V,I_ohmic_back(ii),temperature]=COMSOL_heat(VL,N,M,vlist,T0,operation);
%[V,I_ohmic_back(ii),temperature,R_lateral,R_vertical]=KVL2_old(VL,VR,N,M,V_initial,vmatrix,operation); 
x6(ii) = toc(t);

if I_ohmic_back(ii)>Icomp
    I_ohmic_back(ii)=Icomp;
end

%%%%%%%%%%%%% TAT current
f = zeros(num_v,1);
trial = 0;
finit = ones(num_v,1)*0.5;
[I_TAT_back(ii), f, flag, restart_num] = tatCurrent_uptocb (finit, num_v, N, M, vindex, vlist, ...
    dx, dy, a0, P0, tunnel_rate0, meff, NtL, NtR, EFL0, EFR0, V, mean(temperature(:,1)));   

%%%%%%%%%%%% FN current
EFR = EFR0; EFL = EFL0 - VL; 
TL = mean(temperature(:,1)); TR = mean(temperature(:,M)); %K
I_FN_back(ii) = FNcurrent (N, ...
dx, dy, Lz, me, mh, mt, EFL, EFR, Eg, V, TL, TR);    

I_total_back(ii)=I_TAT_back(ii)+I_ohmic_back(ii)+I_FN_back(ii); 
t = tic;
[vindex,vlist,vmatrix,num_v,ion_list,ion_matrix,num_ion]=DC_MC(vindex,vlist,vmatrix,num_v,N,M,tstep,cycle,ion_list,ion_matrix,num_ion,V,temperature,operation);    
x7(ii) = toc(t);
end
end

figure(2) %0.4789s
semilogy(V_apply,abs(I_ohmic),'go--');
hold on
semilogy(V_apply,abs(I_TAT),'ro--');
semilogy(V_apply,abs(I_total),'ko--');
semilogy(V_apply,abs(I_FN),'bo--');
%semilogy(V_apply,abs(I_ohmic_back),'go--');
legend('I_ohmic','I_TAT','I_total','I_FN','location','southeast');


if operation==1
semilogy(V_back,abs(I_TAT_back),'ro--');
semilogy(V_back,abs(I_ohmic_back),'go--');
semilogy(V_back,abs(I_total_back),'ko--');
semilogy(V_back,abs(I_FN_back),'bo--');
axis([-3,0,1e-12,1e-5]);
else
axis([0,10,1e-12,1e-5]);
end

%%%%%%%%%%%%%%plot after reset

figure(3)
vmatrix3=fliplr(vmatrix);
cla
hold on
p=bar3(vmatrix3',1,'r');
set(p,'facecolor',[1,0.412,0.706])
bar3(zeros(N,M)',1,'w');
axis('equal')
axis([1 N 0.5 M+0.5 0 1.001]); %0.4236s

% ion_matrix = zeros(N,M);
% for kk=1:num_ion
%     if ion_list(kk,1)>0 && ion_list(kk,2)>0 
%        ion_matrix(ion_list(kk,1),ion_list(kk,2))=1;
%     end
% end

figure(3)
ion_matrix3=fliplr(ion_matrix);
hold on
bar3(ion_matrix3',1,'b');
bar3(zeros(N,M)',1,'w');
axis('equal')
axis([1 N 0.5 M+0.5 0 1.001]); 

% figure(5) 
% semilogy(time,I_total,'ro-');

figure(7)
    temperature6=fliplr(temperature);
    surf(temperature6')
    colorbar;
    axis('equal')
    view(0,90)

if operation==1
save('HRS.mat','vindex','vlist','vmatrix','num_v','ion_list','ion_matrix','num_ion')
else
save('LRS.mat','vindex','vlist','vmatrix','num_v','ion_list','ion_matrix','num_ion')
end 