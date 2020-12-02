clc ; clear all ; close all ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Simulation of the Cardiovascular System %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Complete the code in places where the "..." string appears %%


%% Flags:
Plot_flag = 1; % 0 = off , 1 = on
 
%% system parameters:
% Time parameters
HR            = 60 + 6      ;   % [BPM] % 60 + sum of last digits from all members
dt            = 5e-4        ;   % [sec]
Heart_cycles  = 20          ;   % total heart cycles
N_per_cycle   = 1/((HR/60)*dt); % number of steps per heart cycle
t             = 0:dt:(Heart_cycles)*(60/HR);
  
% Heart parameters
V0                          = 15  ;                      % [ml] V0 for Plv calculation
Emax                        = 2.0 ;                      % contractility
N_Systole                   = round(N_per_cycle/3) ;     % number of points per ventricle Systole
N_Diastole                  = N_per_cycle - N_Systole ;  % number of points per ventricle Diastole
E_dia                       = 10/(120-V0);               %
En(1:N_Systole)             = 0.5*(1+sin(((2*pi*(1:N_Systole))/N_Systole)-(pi/2)));
En(N_Systole+1:N_per_cycle) = 0;
E            = max(E_dia,Emax*En) ;                      % Heart Elasticity: combine Systole (elasticity function) & Diastole (uniform value) 
len_E = length(E);
for i = 1:(Heart_cycles-1)
    
    E(len_E*i+1:len_E*(i+1)) = E(1:len_E); % ffffffff
end

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
Vlv_C = zeros(1,round(N_per_cycle*Heart_cycles)); % left ventricle
Va_C  = zeros(1,round(N_per_cycle*Heart_cycles)); % arteries
Vv_C  = zeros(1,round(N_per_cycle*Heart_cycles)); % veins 

Plv_C   = zeros(1,round(N_per_cycle*Heart_cycles));    % left ventricle
Pa_C    = zeros(1,round(N_per_cycle*Heart_cycles));    % arterial capacitor
Pv_C    = zeros(1,round(N_per_cycle*Heart_cycles));    % venous filling 
Pao_C   = zeros(1,round(N_per_cycle*Heart_cycles));

Qlv_C     = zeros(1,round(N_per_cycle*Heart_cycles));   % left ventricle (outflow)
Qp_C      = zeros(1,round(N_per_cycle*Heart_cycles));   % peripheral resistance
Qv_C      = zeros(1,round(N_per_cycle*Heart_cycles));   % ventricle filling (inflow)
        
for CycleIdx = 1 : Heart_cycles % main loop for each heart cycle
    for StepInCycle = 2 : N_per_cycle 
	%calculating all variables for each cycle at N points:
        %Volumes [ml] ...
		Vlv(StepInCycle)  = Vlv(StepInCycle-1) + (Qv(StepInCycle-1)-Qlv(StepInCycle-1))*dt;  % left ventricle
        Va(StepInCycle)   = Va(StepInCycle-1) + (Qlv(StepInCycle-1)-Qp(StepInCycle-1))*dt; % arteries
        Vv(StepInCycle)   = Vv(StepInCycle-1) + (Qp(StepInCycle-1)-Qv(StepInCycle-1))*dt; % veins 
        %Pressures [mmHg] ...
        Plv(StepInCycle)  = E(StepInCycle)*(Vlv(StepInCycle)-V0);    % left ventricle
        Pa(StepInCycle)   = Va(StepInCycle)/Ca;   % arterial capacitor
        Pv(StepInCycle)   = Vv(StepInCycle)/Cv;    % venous filling 
        Pao(StepInCycle)  = max(Pa(StepInCycle),Plv(StepInCycle));   % aorta
        %Flows [ml/sec] ...
        Qlv(StepInCycle)  = max(0,((Pao(StepInCycle)-Pa(StepInCycle))/Ra));    % left ventricle (outflow)
        Qp(StepInCycle)   = (Pa(StepInCycle)-Pv(StepInCycle))/Rp;    % peripheral resistance
        Qv(StepInCycle)  = max(0,((Pv(StepInCycle)-Plv(StepInCycle))/Rv)) ;    % ventricle filling (inflow)
    end
        % Save each variable from the current cycle to a continuos variable:

        %VolumeS [ml] ...
		Vlv_C(N_per_cycle*(CycleIdx-1)+1:N_per_cycle*CycleIdx) = Vlv  ; % left ventricle
        Va_C(N_per_cycle*(CycleIdx-1)+1:N_per_cycle*CycleIdx)  = Va  ; % arteries
        Vv_C(N_per_cycle*(CycleIdx-1)+1:N_per_cycle*CycleIdx)  = Vv  ; % veins 
        %PressureS [mmHg] ...
        Plv_C(N_per_cycle*(CycleIdx-1)+1:N_per_cycle*CycleIdx)   = Plv;    % left ventricle
        Pa_C(N_per_cycle*(CycleIdx-1)+1:N_per_cycle*CycleIdx)    = Pa;    % arterial capacitor
        Pv_C(N_per_cycle*(CycleIdx-1)+1:N_per_cycle*CycleIdx)    = Pv;    % venous filling 
        Pao_C(N_per_cycle*(CycleIdx-1)+1:N_per_cycle*CycleIdx)   = Pao;
        %FlowS [ml/sec] ...
        Qlv_C(N_per_cycle*(CycleIdx-1)+1:N_per_cycle*CycleIdx)     = Qlv;   % left ventricle (outflow)
        Qp_C(N_per_cycle*(CycleIdx-1)+1:N_per_cycle*CycleIdx)      = Qp;   % peripheral resistance
        Qv_C(N_per_cycle*(CycleIdx-1)+1:N_per_cycle*CycleIdx)      = Qv;   % ventricle filling (inflow)

    % Find max and min of each cycle
      [Pao_max(CycleIdx),Pao_max_ind(CycleIdx)] = max(Pao);
      [Pao_min(CycleIdx),Pao_min_ind(CycleIdx)] = min(Pao);
      Pao_min_ind(CycleIdx) = Pao_min_ind(CycleIdx)+(CycleIdx-1)*N_per_cycle;
      Pao_max_ind(CycleIdx) = Pao_max_ind(CycleIdx)+(CycleIdx-1)*N_per_cycle;
      
    % Update the initial variables before the next cycle:  ...
	  %Volume [ml] ...
      Vlv(1)  = Vlv(end) + (Qv(end)-Qlv(end))*dt; % left ventricle
      Va(1)   = Va(end) + (Qlv(end)-Qp(end))*dt;  % arteries
      Vv(1)   = Vv(end) + (Qp(end)-Qv(end))*dt;   % veins
      %Pressure [mmHg] ...
      Plv(1)  = Plv(end) + E(1)*(Vlv(1)-V0);      % left ventricle
      Pa(1)   = Va(1)/Ca;                         % arterial capacitor
      Pv(1)   = Vv(1)/Cv;                         % venous filling 
      Pao(1)  = max(Pa(1),Plv(1));                % aorta
      %Flow [ml/sec] ...
	  Qlv(1)  = max(0,((Pao(1)-Pa(1))/Ra));       % left ventricle (outflow)
      Qp(1)   = (Pa(1)-Pv(1))/Rp;                 % peripheral resistance
      Qv(1)  = max(0,((Pv(1)-Plv(1))/Rv)) ;       % ventricle filling (inflow)
      
end
 
%% Plot
if Plot_flag  
    ... ; % Plotting commands:
   plot(t(1:length(Pao_C)),Pao_C)
   title('Aorta pressure as a funtion of time')
   xlabel('Time (Sec)')
   ylabel('Pressure (mmHG)')
    
   hold on
   scatter(Pao_max_ind*dt,Pao_max)
   scatter(Pao_min_ind*dt,Pao_min)
    
    
end