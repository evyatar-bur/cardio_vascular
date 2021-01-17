%% Use this code to find the wanted k parameters
clear
clc

% These k values will be checked
kp_group = [0.01 0.1 0.4 1];
ki_group = 1:10:100;
kd_group = 1:10:100;

% This cell will be the output, containing the possible parameters
possible_k = {[]};

best_error = 9999;

for Kp = kp_group
    
    disp(Kp)
    
    for Ki = ki_group
        
        disp('  ' + string(Ki))
        
        for Kd = kd_group
            
            error = check(Kp,Ki,Kd);
            
            
            if abs(error) < best_error

            possible_k = [Kp Ki Kd error];
            
            best_error = error;
            
            end
        end
    end
end

%%
k1 = [0.0001 0.0006 0.0001] ;
