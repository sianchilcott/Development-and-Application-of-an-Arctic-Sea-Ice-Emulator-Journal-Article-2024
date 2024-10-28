%% Observationally Constrained Emulator

%% ___________ Step i: Use the observationally constrained Arctic Amplification to calculate the Arctic annual mean temperature ______________________________________________________________________________________________________________________________________________________________________________________________________________________ %% 


new_emulation_global_annual = [];
pi_new_emulation_global_annual = [];
for n = 1:3
    emulation_global_to_annual = save_GMST{n} .* AA_MCMC_3{n};   % Use the AA that is calculated from randomly sampling the global mean temperature, 'p' and 's'
    new_emulation_global_annual{n} = emulation_global_to_annual;


    % Normalise to absolute observations
    avg_obs = mean(Arctic_annual_temp_obs(1:50));


    % Convert anomaly data to absolute temperature by adding the mean of the observed absolute Arctic annual mean temperature
    pi_new_emulation_global_annual{n} = emulation_global_to_annual + avg_obs;
end

new_emulation_global_annual_OC = new_emulation_global_annual;
pi_new_emulation_global_annual_OC = pi_new_emulation_global_annual;





%% ___________ STEP ii: Arctic Seasonal Temperature Parameterisation with Observational Constraints ______________________________________________________________________________________________________________________________________________________________________________________________________________________ %% 

% 1st calculate the Arctic Seasonal Temperature from CMIP6 calibration parameters before adding bias corrections (explained in Chapter 2, Section 2.4.2.2)
% Initialise
x = linspace(0,2*pi,13); 
final_final_save_OC = [];
final_final_save_hold_OC = [];

close
for n = 1:3                         % SSPs
    disp('__________')
        
    OC_AAT = pi_new_emulation_global_annual_OC{n};

    for j = 1:12                    % MODELS
        disp('NEW MODEL')
        for ii = 1:600


    %       Emulated AAT
            siza_AAT = size(OC_AAT);
            AAMST = OC_AAT;
            AAMST = AAMST(ii,:)';



            f1 = parameters_store_new(j,1);
            f2 = parameters_store_new(j,2);
            g1 = parameters_store_new(j,3);
            g2 = parameters_store_new(j,4);
            a1 = parameters_store_new(j,5);
            a2 = parameters_store_new(j,6);
            a3 = parameters_store_new(j,7);


            e = 0.3;
            f = (f1.*AAMST) + f2;
            g = (g1.*AAMST) + g2;
            a = cos((AAMST.*-a1)-a2)+a3;


            starter_ind_temp = -8;                      % Keep parameter 'f' the same until -8 degrees C is reached to allow July and August temperatures to rise at the same rate as observations suggest
            indf = AAMST <= starter_ind_temp;
            f(indf) = (f1.*starter_ind_temp) + f2;


            ind_temp = 6;                               % Keep all parameters the same after 6 degrees C is reached as the Arctic seasonal temperature doesn't change much after this time, this ensures the parameters don't keep changing
            ind = find(AAMST >= ind_temp);
            if isempty(ind) ~= 1
                f(ind) = (f1.*ind_temp) + f2;
                g(ind) = (g1.*ind_temp) + g2;
            end


            AMST_OC = f.*(cos(x .* g - e .* exp(cos(x.^(a))))+((AAMST - mean(f.*(cos(x .* g - e .* exp(cos(x.^(a))) )), 2, 'omitnan')) ./ f));


            final_final_save_hold_OC{ii,j} = AMST_OC(:,1:12);

        end           
    end
    final_final_save_OC{n} = final_final_save_hold_OC;
end

disp('___END RUN___')

OC_tas_no_BC = final_final_save_OC; 




%% Step ii: 2nd part of this section adds the bias corrections

clc
residuals_obs = [-6.0000   -3.0000   -4.0849   -3.8929   -1.4497   -0.0843    1.0000     1.8000    2.0787    0.4880   -2.0996   -2.2037]; % Add residuals (calculated from difference between obs and CMIP6) to emulated Arctic seasonal temperature, 1 value for each month

AMST_OC_bias_correction_hold = [];
AMST_OC_bias_correction = []; 

for n = 1:3
    disp('__________')
    for j = 1:12
        disp('NEW MODEL')
        for ii = 1:600
            AMST_OC_bias_correction_hold{ii,j} = OC_tas_no_BC{n}{ii,j} + residuals_obs;
        end
    end
    AMST_OC_bias_correction{n} = AMST_OC_bias_correction_hold;
end
disp('___END RUN___')





%% Step 2: Reshape Arctic seasonal temperature for plotting

% SAME INDEX RUN AS ABOVE
clc; tic
index_run = [1:12];
amt_mon_rearrange = [];
amt_mon_rearrange_final_OC = [];
amt_mon_rearrange_hold = [];
amt_mon_ensem = [];
for n = 1:3
    disp('__________')
    amt_mon_rearrange_hold = [];
    for j = index_run
        disp('NEW MODEL')

        amt_hold = AMST_OC_bias_correction{n}(:,j);


        amt_mon_ensem = cellfun(@(mm) num2cell(mm,1), amt_hold, 'UniformOutput', false);
        amt_mon_ensem = vertcat(amt_mon_ensem{:});
        amt_mon_ensem = cellfun(@transpose, amt_mon_ensem, 'UniformOutput', false);

        for i = 1:12
            amt_mon_rearrange{j,i} = cell2mat(amt_mon_ensem(:,i));
        end                

    end
    amt_mon_rearrange_final_OC{n} = amt_mon_rearrange;
end
toc
disp('___END RUN___')







%% Step ii: 3rd part of this section finds the difference between the CMIP6 calibrated AMST and the AMST where 'f' has remains the same to -8 degrees C 
% This allows us to inform our 'b' (sensitivity parameter) in our SIA parameterisation



% Initialisation
x = linspace(0,2*pi,13);
AMST_OC_normed_f = [];
AMST_OC_normed_f_hold = []; 
% Calculate the AMST without using the bias correction 
close
for n = 1:3
    disp('__________')
        
    OC_AAT = pi_new_emulation_global_annual_OC{n};

    for j = 1:12
        disp('NEW MODEL')
        for ii = 1:600


    %       Emulated AAT
            siza_AAT = size(OC_AAT);
            AAMST = OC_AAT;
            AAMST = AAMST(ii,:)';



            f1 = parameters_store_new(j,1);
            f2 = parameters_store_new(j,2);
            g1 = parameters_store_new(j,3);
            g2 = parameters_store_new(j,4);
            a1 = parameters_store_new(j,5);
            a2 = parameters_store_new(j,6);
            a3 = parameters_store_new(j,7);


            e = 0.3;
            f = (f1.*AAMST) + f2;
            g = (g1.*AAMST) + g2;
            a = cos((AAMST.*-a1)-a2)+a3;



            ind_temp = 6;
            ind = find(AAMST >= ind_temp);
            if isempty(ind) ~= 1
                f(ind) = (f1.*ind_temp) + f2;
                g(ind) = (g1.*ind_temp) + g2;
            end


            AMST_OC = f.*(cos(x .* g - e .* exp(cos(x.^(a))))+((AAMST - mean(f.*(cos(x .* g - e .* exp(cos(x.^(a))))), 2, 'omitnan')) ./ f));


            AMST_OC_normed_f_hold{ii,j} = AMST_OC(:,1:12);

        end           
    end
    AMST_OC_normed_f{n} = AMST_OC_normed_f_hold;
end

disp('___END RUN___')





% Step ii: add observational constant bias correction
clc
residuals_obs = [-6.0000   -3.0000   -4.0849   -3.8929   -1.4497   -0.0843    1.0000     1.8000    2.0787    0.4880   -2.0996   -2.2037];

AMST_OC_bias_correction_hold_normed_f = [];
AMST_OC_bias_correction_normed_f  = []; 

for n = 1:3
    disp('__________')
    for j = 1:12
        disp('NEW MODEL')
        for ii = 1:600
            AMST_OC_bias_correction_hold_normed_f{ii,j} = AMST_OC_normed_f{n}{ii,j} + residuals_obs;
        end
    end
    AMST_OC_bias_correction_normed_f{n} = AMST_OC_bias_correction_hold_normed_f;
end
disp('___END RUN___')





%% Find difference between AMST with bias corrected 'f' and AMST where 'f' has not been bias corrected


% Ensure both preindustrial (1850-1900) temperatures are the same for comparison
norm_f_pi_diff = [];
for n = 1   % SSP
    for j = 1:12 % models
        for ii = 1:600  % ensembles
            norm_f_diff_hold = [];
            for i = 1:12    % month
                meana = mean(AMST_OC_bias_correction_normed_f{n}{ii,j}(1:50,i)) - mean(AMST_OC_bias_correction{n}{ii,j}(1:50,i));       % find the mean of 
                monthly_corrects = AMST_OC_bias_correction_normed_f{n}{ii,j}(:,i) - meana;
                norm_f_diff_hold = cat(2, norm_f_diff_hold, monthly_corrects);
            end
            norm_f_pi_diff{ii,j} = norm_f_diff_hold;
        end
    end
end


% Calculate the bias correction to add to parameter 'b' in the SIA parameterisation below
norm_f_diff = [];
save_sept_diff = [];
for n = 1
    for j = 1:12
        for ii = 1:600
            norm_f_diff_hold = [];
            for i = 1:12
                fiida = AMST_OC_bias_correction{n}{ii,j}(:,i) - norm_f_pi_diff{ii,j}(:,i);
                meana = mean( fiida(130:205) );

                norm_f_diff_hold = cat(2, norm_f_diff_hold, meana);
            end
            norm_f_diff{ii,j} = norm_f_diff_hold;
        end
    end
end





%% _________Step iii: SIA Parameterisation _____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________


new_vals = [0.3795    0.3961    0.5408   -1.3240
            0.4006   -0.0411    0.5451   -0.9881
            0.2952    0.0687    0.2140   -4.4366
            0.2958   -0.1261    0.1251   -4.4320
            0.4096   -0.0108    1.0000   -0.7661
            0.4010    0.5436    0.5505    0.6585
            0.3473   -0.1438    0.4485    0.6471
            0.3242    0.4528    0.2728   -0.9961
            0.3240   -0.4721    0.4894   -1.1122
            0.2555    0.2036         0   -1.3647
            0.3889   -0.4409         0   -1.7137
            0.3848    0.0896    0.1904   -1.8352]; % (the calibration parameters for reference [a, sm_shift, w1_store, b])


%% SIA Calculation Initialisation

SIA_bias_corrected_OC = [];
SIA_bias_corrected_OC_store = [];
index_run = 1:12;
sia_resid_diff = [0.0    -0.3758   0.0   0.00    0.1000   0.7000   1.9000    2.0    0.7    0.4885    0.4030    0.5]; % SIA constant offset/ bias correction (explained in Chapter 2, Section 2.4.2.2)

tic
for n = 1:3
    disp('_________')    
    disp('_NEW SSP_')
    for j = index_run
        for ii = 1:600

            hnew3 = [];

        %   Counter for plotting
            tt_store = 1:12;


    %       SIA_max
            SIA_max = SIA_max_param_1850(j,1:12) + sia_resid_diff;          % Add the constant offset/ bias correction to the calibrated SIA in 1850/ SIA_max


    %       Emulated AMT
            x_tas_store = final_final_save_MAGICC_MCMC{n}{ii,j};
            x_tas_store2 = final_final_save_MAGICC_MCMC{n}{ii,j};
            x_tas = x_tas_store;


%           Extract Calibration Parameters
            a = Save_sia_cal_params(j,1);
            sm_shift = Save_sia_cal_params(j,2);
            b = Save_sia_cal_params(j,4);
            w1_store = Save_sia_cal_params(j,3);     

            % Double check no weights are greater than 1 or less than 0 (they shouldn't be from the calibration but just incase)
            if w1_store > 1
                w1_store = 1;
            elseif w1_store < 0
                w1_store = 0;
            end


            b = b + norm_f_diff{ii,j};      % Add the bias correction to 'b' the sensitivity parameter to overcome the increase in sensitivity of SIA to Arctic seasonal temperature after the bias correction to 'f' has been applied



            
            % _________ Calculated Weighted temperature to force the SIA parameterisation _________ %
            tas_hold = x_tas_store(end-length(x_tas_store)+1,:);
            dec_previous_year_tas = tas_hold(12);
            dec_tas_test = x_tas(:,12);
            dec_tas_test = circshift(dec_tas_test,1);
            x_tas_shift = [dec_tas_test, x_tas];
            x_tas_shift(1) = dec_previous_year_tas;



            x_tas_hold = [];
            for tt = 1:12 

                if tt == 1
                    hold_x_tas = ((x_tas(:,tt) .* w1_store) + (x_tas_shift(:,1) .* (1-w1_store))) ./ (w1_store+(1-w1_store));
                else
                    hold_x_tas = ((x_tas(:,tt) .* w1_store) + (x_tas(:,tt-1) .* (1-w1_store))) ./ (w1_store+(1-w1_store));
                end

                x_tas_hold(:,tt) = hold_x_tas; 
            end

            x_tas = x_tas_hold; % weighted temperature

            % _______end of weighting scheme ______________________________________________________________________________ %


            % Run Equation
            SIA_bias_corrected_OC_run = ( ((SIA_max+sm_shift) .* (1 + exp((x_tas(1,:)-b)))) - (a.*(x_tas - x_tas(1,:))) ) ./ ( 1 + exp(x_tas-b) );

            SIA_bias_corrected_OC_store{ii,j} = SIA_bias_corrected_OC_run;


        end
    end
    SIA_bias_corrected_OC{n} = SIA_bias_corrected_OC_store;
end
toc
disp('___END RUN___')








%% Step iii: Reshape SIA for plotting: Plot each month

clc; tic
index_run = [1:12];
sia_mon_rearrange = [];
sia_mon_rearrange_final_OC = [];
sia_mon_ensem = [];

for n = 1:3
    disp('___________')
    disp('__NEW SSP__')
    for j = index_run
        % disp('NEW MODEL')

        sia_hold = SIA_bias_corrected_OC{n}(:,j);


        sia_mon_ensem = cellfun(@(mm) num2cell(mm,1), sia_hold, 'UniformOutput', false);
        sia_mon_ensem = vertcat(sia_mon_ensem{:});
        sia_mon_ensem = cellfun(@transpose, sia_mon_ensem, 'UniformOutput', false);


        for i = 1:12
            sia_mon_rearrange{j,i} = cell2mat(sia_mon_ensem(:,i));
        end                

    end
    sia_mon_rearrange_final_OC{n} = sia_mon_rearrange;
end
toc
disp('___END RUN___')










%% ______________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________

%% ____________________ Force with MAGICC GMST _________________________________________________________________________________________________________________________________________________________________________________________


%% Step i: AAT using CMIP6 AA and apply to MAGICC global mean temperature


new_emulation_global_annual = [];
pi_new_emulation_global_annual = [];

for n = 1:3
    for j = 1:12

%       MAGICC GMST
        rw = GT_ensembles{n};

        % CMIP6 Calibrated AA (from CMIP6_Arctic_Amplification_Parameterisation.m)
        new_emulation_global_annual{j,n} = rw .* rc_save(j,1);  

%       Normalise to observations
        avg_obs = mean(Arctic_annual_temp_obs(1:50));

%       Convert anomaly data to absolute tas by normalising to the 
        pi_new_emulation_global_annual{j,n} = (rw .* rc_save(j,1)) + avg_obs;


    end          
end

AAT_OCCAA_anomaly = new_emulation_global_annual;
AAT_OCCAA_absolute = pi_new_emulation_global_annual;




%% Parameters for step ii:


parameters_store_new = [0.4941   -8.4259   -0.0063    0.9309   -0.0780   -0.5773    0.1017
                        0.5613   -7.1602   -0.0071    0.9319   -0.0333    0.1384    0.1348
                        0.6322   -6.0756   -0.0031    0.9697   -0.0893   -0.6312    0.1086
                        0.5826   -6.2933   -0.0042    0.9557   -0.0964   -0.6432    0.1123
                        0.4154   -7.7789   -0.0007    0.9688   -0.0551   -0.0804    0.1267
                        0.5391   -7.4667   -0.0076    0.9505   -0.0875   -0.3211    0.0084
                        0.3465  -11.3899   -0.0009    1.0264   -0.0147    0.6405    0.3556
                        0.4858   -8.3263   -0.0062    0.9415   -0.0924   -0.5403    0.0691
                        0.4479   -9.0050   -0.0040    0.9766   -0.0973   -0.4534    0.0820
                        0.4378   -8.0920   -0.0057    0.9489   -0.1060   -0.3702    0.0921
                        0.5550   -6.4767   -0.0016    0.9877   -0.0218    0.9385    0.5679
                        0.5196   -7.3501   -0.0036    0.9572   -0.0943   -0.5025    0.0998];        % Reminder of the AMST parameters from Arctic_Seasonal_Temperature_Calibration.m



% Initialisation
x = linspace(0,2*pi,13);
AMST_OCCAA_hold = [];
AMST_OCCAA = [];
    
close
for n = 1:3                     % SSP
    disp('__________')
	for j = 1:12                % model
        disp('NEW MODEL')
        for ii = 1:600          % ensembles
        
        
    
            % Emulated AAT for input into AMST parameterisation
            siza_AAT = size(AAT_OCCAA_absolute);
            AAMST = AAT_OCCAA_absolute{j,n};
            AAMST = AAMST(ii,:)';



            % Get input calibrated parameters for each model
            f1 = parameters_store_new(j,1);
            f2 = parameters_store_new(j,2);
            g1 = parameters_store_new(j,3);
            g2 = parameters_store_new(j,4);
            a1 = parameters_store_new(j,5);
            a2 = parameters_store_new(j,6);
            a3 = parameters_store_new(j,7);


            % Calculate each parameter that is dependent on temperature
            e = 0.3;
            f = (f1.*AAMST) + f2;
            g = (g1.*AAMST) + g2;
            a = cos((AAMST.*-a1)-a2)+a3;

                    
            % Bias correct 'f' so it stays the same until -8degrees C is reached so that it rises at the same rate as observations 
            starter_ind_temp = -8;
            indf = AAMST <= starter_ind_temp;
            f(indf) = (f1.*starter_ind_temp) + f2;


            % Keep all calibrations the same after 6 degrees C is reached as the shape of the seasonal curve doesn't change much after 2100 which is about 6 degrees in each model
            ind_temp = 6;
            ind = find(AAMST >= ind_temp);
            if isempty(ind) ~= 1
                f(ind) = (f1.*ind_temp) + f2;
                g(ind) = (g1.*ind_temp) + g2;
            end
                      

            AMST_OCCAA = f.*(cos(x .* g - e .* exp(cos(x.^(a)))) + ((AAMST - mean(f.*(cos(x .* g - e .* exp(cos(x.^(a))))), 2, 'omitnan')) ./ f));
            
            
            AMST_OCCAA_hold{ii,j} = AMST_OCCAA(:,1:12);
            
        end
    end
    AMST_OCCAA{n} = AMST_OCCAA_hold;
end
disp('___END RUN___')




%% Step ii: 2nd part of this section adds the bias corrections


clc
residuals_obs = [-6.0000   -3.0000   -4.0849   -4.0929   -1.4497   -0.0843    1.0000     1.0000    1.2787    0.4880   -2.0996   -2.2037];       % Observational constant offset/ bias correction 


AMST_OCCAA_bias_correction_hold = [];
AMST_OCCAA_bias_correction = [];
for n = 1:3
    disp('__________')
    for j = 1:12
        disp('NEW MODEL')
        for i = 1:600
            AMST_OCCAA_bias_correction_hold{i,j} = AMST_OCCAA{n}{i,j} + residuals_obs;
        end
    end
    AMST_OCCAA_bias_correction{n} = AMST_OCCAA_bias_correction_hold;
end
disp('___END RUN___')








%% Step ii: Reshape AMT for plotting: Plot each month

clc

tic
index_run = [1:12];
amt_mon_rearrange = [];
amt_mon_rearrange_final = [];
amt_mon_rearrange_hold = [];
amt_mon_ensem = [];
AMST_mon_rearrange_OCCAAn = [];

for n = 1:3
    disp('__________')
    for j = index_run
        disp('NEW MODEL')

        amt_hold = final_final_save_MAGICC_CMIP6_AA{n}(:,j);


        amt_mon_ensem = cellfun(@(mm) num2cell(mm,1), amt_hold, 'UniformOutput', false);
        amt_mon_ensem = vertcat(amt_mon_ensem{:});
        amt_mon_ensem = cellfun(@transpose, amt_mon_ensem, 'UniformOutput', false);


        for i = 1:12
            amt_mon_rearrange{j,i} = cell2mat(amt_mon_ensem(:,i));
        end                

    end
    
    AMST_mon_rearrange_OCCAA{n} = amt_mon_rearrange;
end
toc
disp('___END RUN___')






%% Step ii: 3rd part of this section finds the difference between the CMIP6 calibrated AMST and the AMST where 'f' has remains the same to -8 degrees C 
% This allows us to inform our 'b' (sensitivity parameter) in our SIA parameterisation


% Initialisation
x = linspace(0,2*pi,13);
AMST_OCCAA_hold_normed_f = [];
AMST_OCCAA_normed_f = [];

close
for n = 1:3
    disp('__________')
	for j = 1:12
        disp('NEW MODEL')
        for ii = 1:600
        
        
    %       Emulated AAT
            siza_AAT = size(AAT_OCCAA_absolute);
            AAMST = AAT_OCCAA_absolute{j,n};
            AAMST = AAMST(ii,:)';



            f1 = parameters_store_new(j,1);
            f2 = parameters_store_new(j,2);
            g1 = parameters_store_new(j,3);
            g2 = parameters_store_new(j,4);
            a1 = parameters_store_new(j,5);
            a2 = parameters_store_new(j,6);
            a3 = parameters_store_new(j,7);


            e = 0.3;
            f = (f1.*AAMST) + f2;
            g = (g1.*AAMST) + g2;
            a = cos((AAMST.*-a1)-a2)+a3;


            ind_temp = 6;
            ind = find(AAMST >= ind_temp);
            if isempty(ind) ~= 1
                f(ind) = (f1.*ind_temp) + f2;
                g(ind) = (g1.*ind_temp) + g2;
            end
                             

            AMST_OCCAA = f.*(cos(x .* g - e .* exp(cos(x.^(a))))+((AAMST - mean(f.*(cos(x .* g - e .* exp(cos(x.^(a))))), 2, 'omitnan')) ./ f));
            
            
            AMST_OCCAA_hold_normed_f{ii,j} = AMST_OCCAA(:,1:12);
            
        end
    end
    AMST_OCCAA_normed_f{n} = AMST_OCCAA_hold_normed_f;
end
disp('___END RUN___')



%% Step ii: 2nd part of this section adds the constant offsets to the bias corrected 'f' AMST

clc
residuals_obs = [-6.0000   -3.0000   -4.0849   -4.0929   -1.4497   -0.0843    1.0000     1.0000    1.2787    0.4880   -2.0996   -2.2037];


AMST_OCCAA_normed_f_bias_correction_hold = [];
AMST_OCCAA_normed_f_bias_correction = [];
for n = 1:3
    disp('__________')
    for j = 1:12
        disp('NEW MODEL')
        for i = 1:600
            AMST_OCCAA_normed_f_bias_correction_hold{i,j} = AMST_OCCAA_normed_f{n}{i,j} + residuals_obs;
        end
    end
    AMST_OCCAA_normed_f_bias_correction{n} = AMST_OCCAA_normed_f_bias_correction_hold;
end
disp('___END RUN___')





%% Find difference between AMST with bias corrected 'f' and AMST where 'f' has not been bias corrected

norm_f_pi_diff_OCCAA = [];
for n = 1
    for j = 1:12
        for ii = 1:600
            norm_f_diff_hold = [];
            for i = 1:12
                meana = mean(AMST_OCCAA_normed_f_bias_correction{n}{ii,j}(1:50,i)) - mean(AMST_OCCAA_bias_correction{n}{ii,j}(1:50,i));
                monthly_corrects = AMST_OCCAA_normed_f_bias_correction{n}{ii,j}(:,i) - meana;
                norm_f_diff_hold = cat(2, norm_f_diff_hold, monthly_corrects);
            end
            norm_f_pi_diff_OCCAA{ii,j} = norm_f_diff_hold;
        end
    end
end



%% Find difference between normed f tas and correct f tas


clc
norm_f_diff_cmip6aa = [];
for n = 1
    for j = 1:12
        for ii = 1:600
            norm_f_diff_hold = [];
            for i = 1:12
                fiida = AMST_OCCAA_bias_correction{n}{ii,j}(:,i) - norm_f_pi_diff_OCCAA{ii,j}(:,i);
                meana = mean( fiida(130:200) );

                norm_f_diff_hold = cat(2, norm_f_diff_hold, meana); 
            end
            norm_f_diff_cmip6aa{ii,j} = norm_f_diff_hold;
        end
    end
end










%% Step iii: SIA Parameterisation

close; clc


% Initialisation
new_vals = [0.3795    0.3961    0.5408   -1.3240
            0.4006   -0.0411    0.5451   -0.9881
            0.2952    0.0687    0.2140   -4.4366
            0.2958   -0.1261    0.1251   -4.4320
            0.4096   -0.0108    1.0000   -0.7661
            0.4010    0.5436    0.5505    0.6585
            0.3473   -0.1438    0.4485    0.6471
            0.3242    0.4528    0.2728   -0.9961
            0.3240   -0.4721    0.4894   -1.1122
            0.2555    0.2036         0   -1.3647
            0.3889   -0.4409         0   -1.7137
            0.3848    0.0896    0.1904   -1.8352]; %    Reminder if calibrated vals from SIA_Calibration.m 
        

sia_resid_diff = [0.0    -0.3758   0.0   0.00    0.1000   0.7000   1.9000    2.0    0.7    0.4885    0.4030    0.5];             % observational constant offsets

SIA_OCCAA = [];
SIA_OCCAA_store = [];
index_run = 1:12 ;

tic
for n = 1:3
    disp('__________')
    
    for j = index_run
        disp('NEW MODEL')
        
        for i = 1:600
            
            hnew3 = [];

        %   Counter for plotting
            tt_store = 1:12;
            
            

            % SIA_max
            SIA_max = SIA_max_param_1850(j,1:12) + sia_resid_diff;          % Add constant offset to SIA_max
            


            % Emulated AMST
            x_tas_store = AMST_OCCAA_bias_correction{n}{i,j};
            x_tas_store2 = AMST_OCCAA_bias_correction{n}{i,j};
            x_tas = x_tas_store;
            

            % Extract Calibration Parameters
            a = Save_sia_cal_params(j,1);
            sm_shift = Save_sia_cal_params(j,2);
            b = Save_sia_cal_params(j,4);
            w1_store = Save_sia_cal_params(j,3);


            if w1_store > 1
                w1_store = 1;
            elseif w1_store < 0
                w1_store = 0;
            end
            
            b = b + norm_f_diff_cmip6aa{ii,j};       % Add bias correction to 'b' sensitivity parameter
        
                            
            % ________ Weighting Scheme _________
            tas_hold = x_tas_store(end-length(x_tas_store)+1,:);
            dec_previous_year_tas = tas_hold(12);
            dec_tas_test = x_tas(:,12);
            dec_tas_test = circshift(dec_tas_test,1);
            x_tas_shift = [dec_tas_test, x_tas];
            x_tas_shift(1) = dec_previous_year_tas;
            


            x_tas_hold = [];
            for tt = 1:12 

                if tt == 1
                    hold_x_tas = ((x_tas(:,tt) .* w1_store) + (x_tas_shift(:,1) .* (1-w1_store))) ./ (w1_store+(1-w1_store));
                else
                    hold_x_tas = ((x_tas(:,tt) .* w1_store) + (x_tas(:,tt-1) .* (1-w1_store))) ./ (w1_store+(1-w1_store));
                end

                x_tas_hold(:,tt) = hold_x_tas; 
            end

            x_tas = x_tas_hold;         % Weighted temperature array
            % _____________________________________
            


%           Run Equation
            SIA_OCCAA_calc = ( ((SIA_max+sm_shift) .* (1 + exp((x_tas(1,:)-b)))) - (a.*(x_tas - x_tas(1,:))) ) ./ ( 1 + exp(x_tas-b) );

            SIA_OCCAA_store{i,j} = SIA_OCCAA_calc;
            

        end
    end
    SIA_OCCAA{n} = SIA_OCCAA_store;
end
toc
disp('___END RUN___')





%% Step 3: Reshape SIA for plotting: Plot each month


tic
index_run = [1:12];
sia_mon_rearrange = [];
sia_mon_rearrange_final = [];
sia_mon_rearrange_hold = [];
sia_mon_ensem = [];
sia_mon_rearrange_OCCAA = [];

for n = 1:3
    disp('__________')
    for j = index_run
        disp('NEW MODEL')
    
        sia_hold = SIA_OCCAA{n}(:,j);
    
            
        sia_mon_ensem = cellfun(@(mm) num2cell(mm,1), sia_hold, 'UniformOutput', false);
        sia_mon_ensem = vertcat(sia_mon_ensem{:});
        sia_mon_ensem = cellfun(@transpose, sia_mon_ensem, 'UniformOutput', false);


        for i = 1:12
            sia_mon_rearrange{j,i} = cell2mat(sia_mon_ensem(:,i));
        end                

    end
    sia_mon_rearrange_OCCAA{n} = sia_mon_rearrange;
end
toc
disp('___END RUN___')











%% _____________ Plot probabilistic SIA and CMIP6 AA SIA all SSPS on same graph _____________ %


% Plot Initialisation
SSP_set = {'SSP585', 'SSP245', 'SSP126'};
x = 1850:2300;
threshold = [17 83]; % likely rang
shades_ssps_cmip6_aa = [];
hnew_means = [];    
shades_ssps = [];
hnew_means_mcmc = [];
counter = 0;
index_plot = reshape(1:36, 6, 6).';

% Define colour RBG triplets   
colorss1 = [1 0.8 0.8; 0.8 0.8 1; 0 0.6906 0.5];
colorss2 = [1 0 0; 0 0 1; 0 0.3906 0];
 


close
figure(40)
set(gcf, 'Units', 'Inches', 'Position', [.5 .5 20 11])   
for n = 1:3         % SSPs
    
    % Load SIA Emulations
    emul_OC = sia_mon_rearrange_final_OC{n};
    emul_OCCAA = sia_mon_rearrange_OCCAA{n};
    
    ssp = ssp_char_store{n};

    for i = 1:12    % Months
        
        counter = counter + 1;
        h = subplot(6,6,index_plot(counter));

        

        % Median SSP 
        y1 = median(cell2mat(emul_OC(:,i)), 'omitnan');
        y_cmip6_aa = median(cell2mat(emul_OCCAA(:,i)), 'omitnan');


        % Percentile SSP   
        stdY = prctile(cell2mat(emul_OC(:,i)), threshold);
        std_cmip6_aa = prctile(cell2mat(emul_OCCAA(:,i)), threshold);


        % Create curves: OC Likely Range
        curve1 = stdY(2,:);
        curve2 = stdY(1,:);
        
        % Shading std area: OC Likely Range
        x2 = [x, fliplr(x)];
        inBetween = [curve1, fliplr(curve2)];
        
        % Fill SSP: OC Likely Range
        c = colorss1(n,:);
        f = fill(x2, inBetween, c, 'edgecolor','none');
        if n == 2 || n == 1
            set(f,'facealpha',.7)
        else
            set(f,'facealpha',.2)
        end
        hold on
        if i == 12
            shades_ssps = cat(1, shades_ssps, f);
        end



        % Create curves: OCCAA Likely Range
        curve3 = std_cmip6_aa(2,:);
        curve4 = std_cmip6_aa(1,:);
               
        % Shading std area: OCCAA Likely Range
        x_cmip6_aa = [x, fliplr(x)];
        inBetween_cmip6_aa = [curve3, fliplr(curve4)];
        
        % Fill SSP: OCCAA Likely Range
        c = colorss2(n,:);
        f = fill(x_cmip6_aa, inBetween_cmip6_aa, c, 'edgecolor','none');
        set(f,'facealpha',.3)
        hold on
        if i == 12
            shades_ssps_cmip6_aa = cat(1, shades_ssps_cmip6_aa, f);
        end


        % Plot Means from both emulators
        plt = plot(x, y1, '--', 'color', colorss2(n,:), 'LineWidth', 2);
        plt_cmip6_aa = plot(x, y_cmip6_aa, 'color', colorss2(n,:), 'LineWidth', 2);
        if i == 12
            hnew_means = cat(1, hnew_means, plt_cmip6_aa);
            hnew_means_mcmc = cat(1, hnew_means_mcmc, plt);
        end
              

        % Plot SIA Obs Range
        x_obs = 1850:2019;
        hnew_obs = [];
        for m = 3:6
            obs = monthly_obs_SIA{m}{i};
            obs_plot = plot(x_obs, obs, 'color', [0.3 0.3 0.3], 'LineWidth', 1);
        end



        % Plot SIA Obs Mean        
        hnew_obs_mean = [];
        for m = 3:6
            obs = median(cell2mat(sia_mon_obs_rearrange_final(:,i)), 'omitnan');
            obs = movmean(obs, 5);
            obs_plot_2 = plot(x_obs, obs, 'color', 'k', 'LineWidth', 2);
        end
        hnew_obs_mean = cat(1, hnew_obs_mean, obs_plot_2);







        % Plot edits
        title_label = title([month_label(i)], 'HorizontalAlignment', 'Right', 'VerticalAlignment', 'top', 'fontsize', 22);
        ax = gca;
        ax.TitleHorizontalAlignment = 'Right';
        posnew2 = get(title_label, 'Position');

         
        xticks([1900, 2000, 2100, 2200, 2300])
        if counter ~= [6,12,18,24,30,36]
            ax = gca;
            set(ax,'XTick',[]);
        end
                
        ax = gca;
        ax.FontSize = 16;
        
        if counter ~= [1:6, 13:18, 25:30]
            ax = gca;
            set(ax,'YTick',[]);
        end
        
        xtickangle(30)


        % Y Label
        if ismember(counter, [1])
            y2 = ylabel('SIA (million km^2)', 'Position', [1780 -45.4869 -1], 'Fontsize', 20, 'Fontweight', 'bold');
        end
        
       
        posnew = get(h, 'Position');
        if ismember(i, [6,12])
            posnew(2) = posnew(2) - 0.04;
        elseif ismember(i, [5,11])
            posnew(2) = posnew(2) - 0.03;
        elseif ismember(i, [4,10])
            posnew(2) = posnew(2) - 0.02;
        elseif ismember(i, [3,9])
            posnew(2) = posnew(2) - 0.01;
        elseif ismember(i, [2,8])
            posnew(2) = posnew(2) - 0.000;
        else 
            posnew(2) = posnew(2) + 0.01;
        end
        posnew(1) = posnew(1) - 0.085;
        posnew(3) = posnew(3) + 0.01;
        posnew(4) = posnew(4) + 0.047;
        set(h, 'Position', posnew);


        if n == 3
            ylim([-0.01 19])
        else
            ylim([-0.01 17])
        end
        xlim([1850 2300]);
        yline(1)
        
    end
end
legend_holders = [shades_ssps_cmip6_aa; shades_ssps; hnew_means; hnew_means_mcmc; hnew_obs_mean];
legend_labels = ["OCCAA: SSP585"; "OCCAA: SSP245"; "OCCAA: SSP126";...
                 "OC: SSP585"; "OC: SSP245"; "OC: SSP126";...
                 "OCCAA Median: SSP585";  "OCCAA Median: SSP245";  "OCCAA Median: SSP126";...
                 "OC Median: SSP585";  "OC Median: SSP245";  "OC Median: SSP126";...
                 "Observational Median"]; 

tt2 = legend(legend_holders, legend_labels, 'location', 'eastoutside', 'NumColumns', 1, 'fontsize', 13.5);
tt2.Box = 'off'; 


% Posistion legend in sthe correct place
newPosition = [0.89 0.47 0.05 0.05];
newUnits = 'normalized';
set(tt2,'Position', newPosition,'Units', newUnits);



% % Save
% set(gcf, 'PaperOrientation', 'landscape')
% set(gcf,'PaperSize',[51 30]);
% 
% temp=[' Fig2.17 ', '.pdf']; 
% % saveas(gca,temp);
% print('-dpdf', '-painters', ['alt', temp])









