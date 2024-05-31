%% CMIP6 Arctic Amplification (AA) Parameterisation

rc_save = [];  % Save regression coefficients for CMIP6 AA parameterisation (Chapter 2 Section 2.3.2)
AAT_emulation_anomaly = [];  % Save Arcitc annual mean temperature anomaly
AAT_emulation_absolute = []; % Save Arcitc annual mean temperature absolute
for n = 1:3
    for j = 1:12
        
        % CMIP6 Data
        rw = tas_global{n}{j};              % global mean surface air temperature
        raa = tas_arctic_annual{n}{j};      % Arctic annual mean surface air temperature


        % Calculate regression coefficients to get AA
        coefficients = polyfit(rw, raa, 1);
        sl = coefficients(1);
        rc_save{j,n} = sl;
        
        
        % Calculate the Arctic annual temperature anomaly 
        AAT_emulation_anomaly{j,n} = rw .* sl;
        

        % Calculate the Arctic annual temperature absolute  
        CMIP6_AAT = mean(tas_store{n}{j}, 2);
        PI_cmip6 = mean(CMIP6_AAT(1:50)); % preindustrial temperature
        AAT_emulation_absolute{j,n} = (rw .* sl) + PI_cmip6;

    end
end

rc_save = cell2mat(rc_save);