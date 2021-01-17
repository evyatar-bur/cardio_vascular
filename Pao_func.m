function [mean_Pao,Vlv_1,Va_1,Vv_1,Plv_1,Pa_1,Pv_1,Pao_1,Qlv_1,Qp_1,Qv_1] = Pao_func(HR,Emax,Cv,Rp,Vlv_1,Va_1,Vv_1,Plv_1,Pa_1,Pv_1,Pao_1,Qlv_1,Qp_1,Qv_1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Simulation of the Cardiovascular System %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% system parameters:
% Time parameters

dt            = 5e-4        ;   % [sec]
Heart_cycles  = 20          ;   % total heart cycles
N_per_cycle   = round(1/((HR/60)*dt)); % number of steps per heart cycle
t             = 0:dt:(Heart_cycles)*(60/HR);
  
% Heart parameters
V0                          = 15  ;                      % [ml] V0 for Plv calculation
N_Systole                   = round(N_per_cycle/3) ;     % number of points per ventricle Systole
E_dia                       = 10/(120-V0);               % diastolic contractility
En(1:N_Systole)             = 0.5*(1+sin(((2*pi*(1:N_Systole))/N_Systole)-(pi/2))); % systolic contractility
En(N_Systole+1:N_per_cycle) = 0;                         
E                           = max(E_dia,Emax*En) ;       % Heart Elasticity: combine Systole (elasticity function) & Diastole (uniform value) 
len_E                       = length(E);                 

% Expanding E to be in the length of all the cycles combined
for i = 1:(Heart_cycles-1)
    E(len_E*i+1:len_E*(i+1)) = E(1:len_E);
end

% Vascular constants:
Ra = 0.1;  % arterial resistance 
Rv = 0.01; % venous filling resistance
Ca = 2.0;  % arterial compliance
 
 
% Initiate variables:
%Volume [ml]
Vlv(1)  = Vlv_1;  % left ventricle
Va(1)   = Va_1;  % arteries
Vv(1)   = Vv_1; % veins 
%Pressure [mmHg]
Plv(1)  = Plv_1;    % left ventricle
Pa(1)   = Pa_1;   % arterial capacitor
Pv(1)   = Pv_1;    % venous filling 
Pao(1)  = Pao_1;  % aorta
%Flow [ml/sec]
Qlv(1)  = Qlv_1;    % left ventricle (outflow)
Qp(1)   = Qp_1;    % peripheral resistance
Qv(1)   = Qv_1;    % ventricle filling (inflow)
 
%% Main Program
        

    for StepInCycle = 2 : N_per_cycle 
	%calculating all variables for each cycle at N points:
        %Volumes [ml] 
		Vlv(StepInCycle)  = Vlv(StepInCycle-1) + (Qv(StepInCycle-1)-Qlv(StepInCycle-1))*dt; % left ventricle
        Va(StepInCycle)   = Va(StepInCycle-1) + (Qlv(StepInCycle-1)-Qp(StepInCycle-1))*dt;  % arteries
        Vv(StepInCycle)   = Vv(StepInCycle-1) + (Qp(StepInCycle-1)-Qv(StepInCycle-1))*dt;   % veins 
        %Pressures [mmHg] 
        Plv(StepInCycle)  = E(StepInCycle)*(Vlv(StepInCycle)-V0);    % left ventricle
        Pa(StepInCycle)   = Va(StepInCycle)/Ca;                      % arterial capacitor
        Pv(StepInCycle)   = Vv(StepInCycle)/Cv;                      % venous filling 
        Pao(StepInCycle)  = max(Pa(StepInCycle),Plv(StepInCycle));   % aorta
        %Flows [ml/sec] 
        Qlv(StepInCycle)  = max(0,((Pao(StepInCycle)-Pa(StepInCycle))/Ra));    % left ventricle (outflow)
        Qp(StepInCycle)   = (Pa(StepInCycle)-Pv(StepInCycle))/Rp;              % peripheral resistance
        Qv(StepInCycle)   = max(0,((Pv(StepInCycle)-Plv(StepInCycle))/Rv)) ;   % ventricle filling (inflow)
    end

      
    % Update the initial variables before the next cycle: 
	  %Volume [ml] 
      Vlv_1  = Vlv(end) + (Qv(end)-Qlv(end))*dt; % left ventricle
      Va_1   = Va(end) + (Qlv(end)-Qp(end))*dt;  % arteries
      Vv_1   = Vv(end) + (Qp(end)-Qv(end))*dt;   % veins
      %Pressure [mmHg]
      Plv_1  = Plv(end) + E(1)*(Vlv(1)-V0);      % left ventricle
      Pa_1   = Va(1)/Ca;                         % arterial capacitor
      Pv_1   = Vv(1)/Cv;                         % venous filling 
      Pao_1  = max(Pa(1),Plv(1));                % aorta
      %Flow [ml/sec] ...
	  Qlv_1  = max(0,((Pao(1)-Pa(1))/Ra));       % left ventricle (outflow)
      Qp_1   = (Pa(1)-Pv(1))/Rp;                 % peripheral resistance
      Qv_1  = max(0,((Pv(1)-Plv(1))/Rv)) ;       % ventricle filling (inflow)

% Calculating mean pao      
mean_Pao = mean(Pao);

end
