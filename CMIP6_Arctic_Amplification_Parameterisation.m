%% CMIP6 Arctic Amplification (AA) Parameterisation

rc_save = [];  % Save regression coefficients for CMIP6 AA parameterisation (Section 2.3)
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
        rc_save{j,n} = sl; % sl = beta (ß)
        
        
        % Calculate the Arctic annual temperature anomaly 
        AAT_emulation_anomaly{j,n} = rw .* sl;
        

        % Calculate the Arctic annual temperature absolute  
        CMIP6_AAT = mean(tas_store{n}{j}, 2);
        PI_cmip6 = mean(CMIP6_AAT(1:50)); % preindustrial temperature
        AAT_emulation_absolute{j,n} = (rw .* sl) + PI_cmip6;

    end
end

rc_save = cell2mat(rc_save);




%% Figure 2


% Re-shape variables for plotting
cmip6_anom = [tas_arctic_annual{1}, tas_arctic_annual{2}, tas_arctic_annual{3}]; % CMIP6 AAST (for comparison to my emulation)
cmip6_anom_and_absolute = {cmip6_anom, AAT_ABS_PI};                              % Combine CMIP6 anomaly and absolute temperature
emul_anom_and_absolute = {AAT_emulation_anomaly, AAT_emulation_absolute};        % Emulated anomaly and absolute AAST


% Plot Initialisation   
threshold = [17 83]; % Likely Range
ssp_legend = {'SSP5-8.5', 'SSP2-4.5', 'SSP1-2.6'};
hnew_VL_range_emul = [];
hnew_VL_range_cmip6 = [];
hnew_mean_emul = [];
hnew_mean_cmip6 = [];



% Start plotting
close all; clc
figure(40)
set(gcf, 'Units', 'Inches', 'Position', [.4 .4 19 10])
counter = 0;
for tt = 1:2                                                        % when tt=1 we plot the anomaly temperature and tt=2 we plot the absolute

    % counter = counter + 1;
   

    emul_anomaly_or_abso_data = emul_anom_and_absolute{tt};         % Emulation
    cmip6_anomaly_or_abso_data = cmip6_anom_and_absolute{tt};       % CMIP6 for comparison           

    
    for n = 1:3     % Plot each SSP     
         
        counter = counter + 1;

        h = subplot(2,3,counter);

        Emul = emul_anomaly_or_abso_data(:,n);       
        CMIP6 = cmip6_anomaly_or_abso_data(:,n);         


        % Initialise the mean of emulation and CMIP6
        emul_mean = mean(cell2mat(Emul));
        cmip6_mean = mean(cell2mat(CMIP6));

        
        % Initialise the year to plot the timeseries
        x = 1850:2100;

        
        % Calculate the likely percentile range
        percentile_Emul = prctile(cell2mat(Emul), threshold);      % Emulation
        percentile_CMIP6 = prctile(cell2mat(CMIP6), threshold);    % CMIP6

        
        % Initialise the percentile range to plot the shaded likely range 
        curve1_emul = percentile_Emul(2,:);
        curve2_emul = percentile_Emul(1,:);

        curve3_cmip6 = percentile_CMIP6(2,:);
        curve4_cmip6 = percentile_CMIP6(1,:);


        % Plot the shaded area between the upper and lower likely percentiles
        x2 = [x, fliplr(x)];
        inBetween_emul = [curve1_emul, fliplr(curve2_emul)];
        inBetween_cmip6 = [curve3_cmip6, fliplr(curve4_cmip6)];


        % Fill in the area between the upper and lower likely percentiles
        % Emulation
        c = colmat(n,:);
        f = fill(x2, inBetween_emul, c, 'edgecolor','none');
        set(f,'facealpha',.3)
        hold on

        
        % CMIP6
        c = [0 0 0];
        f2 = fill(x2, inBetween_cmip6, c, 'edgecolor','none'); 
        set(f2,'facealpha',.3)

        
        % Plot the mean:
        d = colmat(n,:);
        emul = plot(x, emul_mean, 'color', d, 'LineWidth', 3);     % Emulation

        c = [0 0 0];
        cmip6 = plot(x, cmip6_mean, 'color', c, 'LineWidth', 3);   % CMIP6                  


        

        % Title for each panel
        if ismember(tt, [2])
            ending = ['Absolute'];
        else
            ending = ['Anomaly'];
        end
        title([ending, ': ' ssp_legend{n}])

        
        % Plot labels
        if tt == 1 && n == 1
            ylab = ylabel( 'Emulated Arctic Annual Temperature (\circC)' , 'Position', [1.82e+03 -5 -1]);
        end
        if ismember(counter, 4:6)
            xlabel('Year')
        end


        % Save plots for legend
        if tt == 1
            hnew_VL_range_emul = cat(1, hnew_VL_range_emul, f);
            hnew_mean_emul = cat(1, hnew_mean_emul, emul);  

            if n == 3
                hnew_mean_cmip6 = cat(1, hnew_mean_cmip6, cmip6);
                hnew_VL_range_cmip6 = cat(1, hnew_VL_range_cmip6, f2);
            end               
        end

        h.Position(1) = h.Position(1) - 0.06;



        % Set fontsize of plot labels
        set(gca,'FontSize', 19)
        xlim([1850 2100])

        grid

    end
end

% Save
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperSize',[49 25.5]);


% temp=[' Fig2.2 ','.pdf']; 
% saveas(gca,temp);



