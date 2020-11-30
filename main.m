clc ; clear all ; close all ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Simulation of the Cardiovascular System %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Complete the code in places where the "..." string appears %%


%% Flags:
Plot_flag = 1; % 0 = off , 1 = on
 
%% system parameters:
% Time parameters
HR           = 60 + 6 ; % [BPM] % 60 + sum of last digits from all members
dt           = 5e-4        ; % [sec]
Heart_cycles = 20          ; % total heart cycles
N_per_cycle      = 1/((HR/60)*dt)     ; % number of steps per heart cycle
t = 1:dt:Heart_cycles;
  
% Heart parameters
V0                          = 15  ;                      % [ml] V0 for Plv calculation
Emax                        = 2.0 ;                      % contractility
N_Systole                   = round(N_per_cycle/3) ;     % number of points per ventricle Systole
N_Diastole                  = N_per_cycle - N_Systole ;  % number of points per ventricle Diastole
E_dia                       = 10/(120-V0);               %
En(1:N_Systole)             = 0.5*(1+sin(((2*pi*(1:N_Systole))/N_Systole)-(pi/2)));
En(N_Systole+1:N_per_cycle) = 0;
E            = max(E_dia,Emax*En) ;                      % Heart Elasticity: combine Systole (elasticity function) & Diastole (uniform value) 
 
% Vascular constants:
Ra = 0.1;  % arterial resistance 
Rp = 1.0;  % peripheral resistance
Rv = 0.01; % venous filling resistance
Ca = 2.0;  % arterial compliance
Cv = 300.0; % venous compliance 
 
% Initiate variables:
%Volume [ml]
Vlv(1)  = 120;  % left ventricle
Va(1)   = 270;  % arteries
Vv(1)   = 2700; % veins 
%Pressure [mmHg]
Plv(1)  = 0;    % left ventricle
Pa(1)   = 70;   % arterial capacitor
Pv(1)   = 9;    % venous filling 
Pao(1)  = 100;   % aorta
%Flow [ml/sec]
Qlv(1)  = 0;    % left ventricle (outflow)
Qp(1)   = 0;    % peripheral resistance
Qv(1)   = 0;    % ventricle filling (inflow)
 
%% Main Program
for CycleIdx = 1 : Heart_cycles % main loop for each heart cycle
    for StepInCycle = 2 : N_per_cycle 
	%calculating all variables for each cycle at N points:
        %Volumes [ml] ...
		Vlv(StepInCycle)  = Vlv(StepInCycle-1) + (Qv(StepInCycle-1)-Qlv(StepInCycle-1))*dt;  % left ventricle
        Va(StepInCycle)   = Va(StepInCycle-1) + (Qlv(StepInCycle-1)-Qp(StepInCycle-1))*dt; % arteries
        Vv(StepInCycle)   = Vv(StepInCycle-1) + (Qp(StepInCycle-1)-Qv(StepInCycle-1))*dt; % veins 
        %Pressures [mmHg] ...
        Plv(1)  = 0;    % left ventricle
        Pa(1)   = 70;   % arterial capacitor
        Pv(1)   = 9;    % venous filling 
        Pao(1)  = 100;   % aorta
        %Flows [ml/sec] ...
        Qlv(1)  = 0;    % left ventricle (outflow)
        Qp(1)   = 0;    % peripheral resistance
        Qv(1)   = 0;    % ventricle filling (inflow)
    end
        % Save each variable from the current cycle to a continuos variable:
        ...
        %VolumeS [ml] ...
		...
        %PressureS [mmHg] ...
        ...
        %FlowS [ml/sec] ...
        ...


    % Update the initial variables before the next cycle:  ...
	  %Volume [ml] ...
      ...
      %Pressure [mmHg] ...
      ...
      %Flow [ml/sec] ...
	  ...

end
 
%% Plot
if Plot_flag  
    ... ; % Plotting commands:
end