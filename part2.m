%% Part 2 - PID controller

clear
clc

% Setting the wanted mean pao to be the stable value without controller and
% without interference
wanted_pao = 72.2803 ;

HR     = 66          ;      % [BPM] % 60 + sum of last digits from all members

% Setting initial values
Emax   = 2           ;      % max contractility
Cv     = 300.0       ;      % venous compliance 
Rp     = 1.0         ;      % peripheral resistance
 

% Initiate variables:
%Volume [ml]
Vlv_1  = 120;  % left ventricle
Va_1   = 270;  % arteries
Vv_1   = 2700; % veins 
%Pressure [mmHg]
Plv_1  = 0;    % left ventricle
Pa_1   = 70;   % arterial capacitor
Pv_1   = 9;    % venous filling 
Pao_1  = 100;  % aorta
%Flow [ml/sec]
Qlv_1  = 0;    % left ventricle (outflow)
Qp_1   = 0;    % peripheral resistance
Qv_1   = 0;    % ventricle filling (inflow)


% Setting controller constants
Kp = 0.005  ;
Ki = 0.003  ;
Kd = 0.0002 ;


% Setting initial error to be zero
error        = 0 ;
last_error   = 0 ;
error_sum    = 0 ;

% setting total number of cycles, and interference index
Heart_cycles = 250;
interferenceIdx = 150;


% Preallocating mean pao vector for speed
mean_pao = zeros(1,Heart_cycles);

for CycleIdx = 1 : Heart_cycles % main loop for each heart cycle
    
    % After 20 cycles we add an interference - blood loss of 5% 
    if CycleIdx == interferenceIdx
        
        % Setting interference (blood loss of 5%)
        Vlv_1 = 0.95*Vlv_1  ;
        Va_1  = 0.95*Va_1   ;
        Vv_1  = 0.95*Vv_1   ;
        
    end

    % Calculating the PID 
    P = Kp*error                  ; % Proportional controller 
    I = Ki*error_sum              ; % Intergrator controller
    D = Kd*(error - last_error)   ; % Derivative controller
    
    PID = P + I + D + 1           ; % Setting PID

    
    % Parameters update
    Emax = max(0,Emax*PID)        ;
    Cv   = max(0.001,Cv*PID)      ;
    Rp   = max(0.001,Rp*PID)      ;
    
    % Computing mean pao
    [mean_pao(CycleIdx),Vlv_1,Va_1,Vv_1,Plv_1,Pa_1,Pv_1,Pao_1,Qlv_1,Qp_1,Qv_1] = Pao_func(HR,Emax,Cv,Rp,Vlv_1,Va_1,Vv_1,Plv_1,Pa_1,Pv_1,Pao_1,Qlv_1,Qp_1,Qv_1);
    
    % Updating error,last error and error sum
    last_error = error;
    error = wanted_pao - mean_pao(CycleIdx);
    error_sum = error_sum + error;
    
end


%% Computing mean pao without PID controller for reference


% Allocating for speed
mean_Pao_without = zeros(1,Heart_cycles);


% Retrieving parameters to initial values
Emax   = 2     ;        % max contractility
Cv     = 300.0 ;        % venous compliance 
Rp     = 1.0   ;        % peripheral resistance


% Initiate variables:
%Volume [ml]
Vlv_1  = 120;  % left ventricle
Va_1   = 270;  % arteries
Vv_1   = 2700; % veins 
%Pressure [mmHg]
Plv_1  = 0;    % left ventricle
Pa_1   = 70;   % arterial capacitor
Pv_1   = 9;    % venous filling 
Pao_1  = 100;  % aorta
%Flow [ml/sec]
Qlv_1  = 0;    % left ventricle (outflow)
Qp_1   = 0;    % peripheral resistance
Qv_1   = 0;    % ventricle filling (inflow)


% Running cycles without PID controller
for CycleIdx = 1 : Heart_cycles % main loop for each heart cycle
    
    % After 20 cycles we add an interference - blood loss of 5% 
    if CycleIdx == interferenceIdx
        
        % Setting interference (blood loss of 5%)
        Vlv_1 = 0.95*Vlv_1  ;
        Va_1  = 0.95*Va_1   ;
        Vv_1  = 0.95*Vv_1   ;
        
    end
   
    % Computing mean pao
    [mean_Pao_without(CycleIdx),Vlv_1,Va_1,Vv_1,Plv_1,Pa_1,Pv_1,Pao_1,Qlv_1,Qp_1,Qv_1] = Pao_func(HR,Emax,Cv,Rp,Vlv_1,Va_1,Vv_1,Plv_1,Pa_1,Pv_1,Pao_1,Qlv_1,Qp_1,Qv_1);
   
end



%% Plotting the computed mean paos together


plot(mean_pao)
hold on
plot(mean_Pao_without)

title('Mean pao as a function of heart cycles')
xlabel('Heart Cycle [au]')
ylabel('mean pao [mmHg]')

legend('With PID controller','Without PID controller')
