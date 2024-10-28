%% Sensitivity Emulator


%% ______ Fig 7 (Section 3.2) (only March and September) _______ % 


% Initialise Plot  
% movmean_ind = 1;
movmean_ind = 20;
years = 1850:2300;
index_1974_2014 = [130:165]; % actually 1979-2014
% index_1974_2014 = [130:170]; % actually 1979-2014
emul_ind = [0:600:7200];
emul_ind(1) = 1;
plot_mods = [4.5:0.1:5.6];
plot_mods2 = [6.5:0.1:7.6];
colorss1 = [0.8 0.8 1; 1 0.8 0.8; 0 0.6906 0.5];
colorss2 = [0 0 1; 1 0 0; 0 0.3906 0];

xticks_hold = [0];
% Get xticks ready
for ii = 1:12
    xticks_hold = [xticks_hold, ii-0.2, ii, ii+0.2];
end

% Initialise empty vectors
obs_sens = [];
obs_sens_co2 = [];
emul_cmip6_sens_co2 = [];
emul_cmip6_sens = [];
cmip6_sens = [];
cmip6_sens_co2 = [];
magicc_sens_OCCAA_co2 = [];
magicc_sens_OCCAA = [];
magicc_sens_OC_co2 = [];
magicc_sens_OC = [];

clc
for month_ind = [3,9]

    % ________ OBS ________ %
    % Calculate the observational Sensitivity: GW
    obs_sens_hold = [];
    clear coefficients
    for nn = [1,3:6]    
        for m = [3:6]
            GT = global_mean_redshaped{nn};
            GT = movmean(GT, movmean_ind);
            GT = GT(index_1974_2014);

            SIA_obs = monthly_obs_SIA{m}{month_ind};
            SIA_obs = movmean(SIA_obs, movmean_ind);
            SIA_obs = SIA_obs(index_1974_2014);

            ind = isnan(SIA_obs);
            if isempty(ind) ~= 1
                SIA_obs = SIA_obs(ind==0);
                GT = GT(ind==0);
            else
                SIA_obs = SIA_obs;
                GT = GT;
            end

            coefficients = polyfit(GT, SIA_obs, 1);
            obs_sens_hold{nn,m} = coefficients(1);
        end
    end
    obs_sens{month_ind} = cell2mat(obs_sens_hold);


    % Calculate the observational Sensitivity: co2
    obs_sens_hold = [];
    clear coefficients 
    for m = [3:6]
        GT = cumsum_co2_emissions_int{1}(101:351);
        GT = movmean(GT, movmean_ind);
        GT = GT(index_1974_2014);

        SIA_obs = monthly_obs_SIA{m}{month_ind};
        SIA_obs = movmean(SIA_obs, movmean_ind);
        SIA_obs = SIA_obs(index_1974_2014);

        ind = isnan(SIA_obs);
        if isempty(ind) ~= 1
            SIA_obs = SIA_obs(ind==0);
            GT = GT(ind==0);
        else
            SIA_obs = SIA_obs;
            GT = GT;
        end

        coefficients = polyfit(GT, SIA_obs, 1);
        obs_sens_hold{m} = coefficients(1);
    end
    obs_sens_co2_hold = cell2mat(obs_sens_hold);
    obs_sens_co2{month_ind} = (obs_sens_co2_hold .* 1e12) / 1e9;

    % ________ OBS END ________ %



    % ________ CMIP6 SENS ________ %
    % Calculate the CMIP6 Sensitivity: GW
    clear coefficients
    cmip6_sens_hold = [];
    for n = 1
        for j = 1:12
            SIA_cmip6 = updated_hist_sia_annual_curve_all_models{j,n}(:,month_ind);
            SIA_cmip6 = movmean(SIA_cmip6, movmean_ind);
            SIA_cmip6 = SIA_cmip6(index_1974_2014);

            GT_cmip6 = tas_global{n}{j};
            GT_cmip6 = movmean(GT_cmip6, movmean_ind);
            GT_cmip6 = GT_cmip6(index_1974_2014);

            coefficients = polyfit(GT_cmip6, SIA_cmip6, 1);
            cmip6_sens_hold{j,n} = coefficients(1); 
        end
    end
    cmip6_sens{month_ind} = cell2mat(cmip6_sens_hold);


        
    % Calculate the CMIP6 Sensitivity: CO2
    clear coefficients
    cmip6_sens_hold = [];
    for n = 1
        for j = 1:12
            SIA_cmip6 = updated_hist_sia_annual_curve_all_models{j,n}(:,month_ind);
            SIA_cmip6 = movmean(SIA_cmip6, movmean_ind);
            SIA_cmip6 = SIA_cmip6(index_1974_2014);

            GT_cmip6 = cumsum_co2_emissions_int{n}(101:351);
            GT_cmip6 = movmean(GT_cmip6, movmean_ind);
            GT_cmip6 = GT_cmip6(index_1974_2014);

            coefficients = polyfit(GT_cmip6, SIA_cmip6, 1);
            cmip6_sens_hold{j,n} = coefficients(1); 
        end
    end
    cmip6_sens_co2_hold = cell2mat(cmip6_sens_hold);
    cmip6_sens_co2{month_ind} = (cmip6_sens_co2_hold .* 1e12) / 1e9;
    % ________ CMIP6 SENS END ________ %


    
    % ________ OCCAA SENS ________ %
    % CO2
    clear coefficients
    magicc_sens_hold_co2 = [];
    for n = 1
        for j = 1:12
            for ii = 1:600
                SIA_emul = sia_mon_rearrange_final_CMIP6_AA{n}(:,month_ind);
                SIA_emul = movmean(SIA_emul{j}(ii,:), movmean_ind);
                SIA_emul = SIA_emul(index_1974_2014);

                GT_emul = cumsum_co2_emissions_int{n}(101:351);
                GT_emul = movmean(GT_emul, movmean_ind);
                GT_emul = GT_emul(index_1974_2014);

                coefficients = polyfit(GT_emul, SIA_emul, 1);
                magicc_sens_hold_co2{ii,j} = coefficients(1);
            end  
        end  
    end
    magicc_sens_OCCAA_co2_hold = cell2mat(magicc_sens_hold_co2); 
    magicc_sens_OCCAA_co2{month_ind} = (magicc_sens_OCCAA_co2_hold .* 1e12) / 1e9;

    
    % GW
    clear coefficients
    magicc_sens_hold = [];
    for n = 1
        for j = 1:12
            for ii = 1:600
                SIA_emul = sia_mon_rearrange_final_CMIP6_AA{n}(:,month_ind);
                SIA_emul = movmean(SIA_emul{j}(ii,:), movmean_ind);
                SIA_emul = SIA_emul(index_1974_2014);

                GT = GT_ensembles{n};
                GT_emul = movmean(GT(ii,:), movmean_ind);
                GT_emul = GT_emul(index_1974_2014);

                coefficients = polyfit(GT_emul, SIA_emul, 1);
                magicc_sens_hold{ii,j} = coefficients(1);
            end  
        end
    end
    magicc_sens_OCCAA{month_ind} = cell2mat(magicc_sens_hold);   
    % ________ OCCAA SENS END ________ %



    % ________ OC SENS ________ %
    % CO2
    clear coefficients
    magicc_sens_hold_co2 = [];
    for n = 1
        for j = 1:12
            for ii = 1:600
                SIA_emul = sia_mon_rearrange_final_MCMC{n}(:,month_ind);
                SIA_emul = movmean(SIA_emul{j}(ii,:), movmean_ind);
                SIA_emul = SIA_emul(index_1974_2014);

                GT_emul = cumsum_co2_emissions_int{n}(101:351);
                GT_emul = movmean(GT_emul, movmean_ind);
                GT_emul = GT_emul(index_1974_2014);

                coefficients = polyfit(GT_emul, SIA_emul, 1);
                magicc_sens_hold_co2{ii,j} = coefficients(1);
            end  
        end  
    end
    magicc_sens_OC_co2_hold = cell2mat(magicc_sens_hold_co2); 
    magicc_sens_OC_co2{month_ind} = (magicc_sens_OC_co2_hold .* 1e12) / 1e9;
    

    % GW
    clear coefficients
    magicc_sens_hold = [];
    for n = 1
        for j = 1:12
            for ii = 1:600
                SIA_emul = sia_mon_rearrange_final_MCMC{n}(:,month_ind);
                SIA_emul = movmean(SIA_emul{j}(ii,:), movmean_ind);
                SIA_emul = SIA_emul(index_1974_2014);

                GT = GT_ensembles{n};
                GT_emul = movmean(GT(ii,:), movmean_ind);
                GT_emul = GT_emul(index_1974_2014);

                coefficients = polyfit(GT_emul, SIA_emul, 1);
                magicc_sens_hold{ii,j} = coefficients(1);
            end  
        end
    end
    magicc_sens_OC{month_ind} = cell2mat(magicc_sens_hold);
    % ________ OC SENS END ________ %

end



%% ________ PLOT SENSITIVITIY ________ %%

mods_legend = [];

close; clc
figure(2)
set(gcf, 'Units', 'Inches', 'Position', [.4 .4 19 10])
hold on

%________GW PLOTS________%
for month_ind = [3,9]

    hold on
    h = subplot(1,2,1);
    
    if month_ind == 3
        ind = 0.2;
    else
        ind = 0.55; 
    end

    % Dirk's 
    if month_ind == 9
        boxplot([-3.81, -4.39], 'positions', [ind-0.045-0.045], 'color', [0    0.5000    0.5000], 'Labels', [' '], 'plotstyle','compact'); hold on 
    end

  
    
    
    % OBS GMST
    if month_ind == 9
        % Plausible range
        boxplot([-2.73, -5.28], 'positions', [ind-0.045], 'color', colorss1(2,:), 'plotstyle','compact'); hold on 

        % Actual obs
        obs_sens_bp = boxplot(obs_sens{month_ind}(:), 'positions', [ind], 'color', [0 0 0]+0.7, 'plotstyle', 'compact'); hold on
        hAx = gca;
        lines = hAx.Children;
        lines = lines(1);
        bp_cmip6 = findobj(lines, 'tag', 'Box');
        bp_cmip6.YData(1,1) = quantile(obs_sens{month_ind}(:), 0.25);
        bp_cmip6.YData(1,2) = quantile(obs_sens{month_ind}(:), 0.75);
    else

        % Actual obs
        obs_sens_bp = boxplot(obs_sens{month_ind}(:), 'positions', [ind], 'color', [0 0 0]+0.7, 'plotstyle', 'compact'); hold on
        hAx = gca;
        lines = hAx.Children;
        lines = lines(1);
        bp_cmip6 = findobj(lines, 'tag', 'Box');
        bp_cmip6.YData(1,1) = quantile(obs_sens{month_ind}(:), 1.0);
        bp_cmip6.YData(1,2) = quantile(obs_sens{month_ind}(:), 0.0);
    end

 
    
        
    % CMIP6 GMST 
    boxplot(cmip6_sens{month_ind}, 'positions', [ind+0.045], 'color', [1 0 0], 'plotstyle','compact'); hold on    % CMIP6
    hAx = gca;
    lines = hAx.Children;
    lines = lines(1);
    bp_cmip6 = findobj(lines, 'tag', 'Box');
    bp_cmip6.YData(1,1) = quantile( cmip6_sens{month_ind}, 0.25);
    bp_cmip6.YData(1,2) = quantile( cmip6_sens{month_ind}, 0.75);
    

       
    % % OCCAA GMST 
    % OCCAA_GMST = boxplot(magicc_sens_OCCAA{month_ind}(:), 'positions', [ind+0.045], 'color', [0    0.5000    0.5000], 'plotstyle','compact'); hold on  
    % set(OCCAA_GMST(length(OCCAA_GMST),:),'Visible','off')     
    % hAx = gca;
    % lines = hAx.Children;
    % lines = lines(1);
    % bp_cmip6 = findobj(lines, 'tag', 'Box');
    % bp_cmip6.YData(1,1) = quantile( magicc_sens_OCCAA{month_ind}(:), 0.25);
    % bp_cmip6.YData(1,2) = quantile( magicc_sens_OCCAA{month_ind}(:), 0.75);

 
    % OC GMST 
    OC_GMST = boxplot(magicc_sens_OC{month_ind}(:), 'positions', [ind+0.09], 'color', colorss1(1,:), 'plotstyle','compact'); hold on   
    set(OC_GMST(length(OC_GMST),:),'Visible','off')  
    hAx = gca;
    lines = hAx.Children;
    lines = lines(1);
    bp_cmip6 = findobj(lines, 'tag', 'Box');
    bp_cmip6.YData(1,1) = quantile( magicc_sens_OC{month_ind}(:), 0.25);
    bp_cmip6.YData(1,2) = quantile( magicc_sens_OC{month_ind}(:), 0.75);   

       
   
    
    
    % Set width of boxplots
    a = get(get(gca,'children'),'children');  
    for ii = 1:length(a)
        t = get(a{ii},'tag'); 
        idx=strcmpi(t,'box'); 
        a_new = a{ii};
        boxes=a_new(idx);         
        set(boxes,'linewidth',20); 
    end
    
    % Add plot labels and edits
    title([num2str(years(index_1974_2014(1))), '-', num2str(years(index_1974_2014(end))), ' Sensitivity (dSIA/dGTA)'], 'fontsize', 22, 'interpreter', 'none')
    
    if ismember(month_ind, [1,5,9])
        ylabel('million km^2/ \circC')
    end
    
    
end
h.Position(1) = h.Position(1) + 0.06;
set(gca, 'fontsize', 20)
    
xticks([0.2, 0.5])
xticklabels( [string({month_label{3}}), string({month_label{9}})] );
length_curve_2 = 0.75;
xlim([0 length_curve_2])    
ylim([-6.8 0])  
grid
h.Position(4) = h.Position(4) - 0.1;
h.Position(3) = h.Position(3) - 0.03;
h.Position(3) = h.Position(3) - 0.03;

% hLegend = legend(findall(gca,'Tag','Box'), {'OC Emulator', 'OCCAA Emulator', 'CMIP6', 'Niederdrenk and Notz, (2018)', 'Plausible Range (from Observations)'}, 'location', 'southwest', 'fontsize', 18);
% hLegend.Box = 'off';

hLegend = legend(findall(gca,'Tag','Box'), {'This Study', 'CMIP6', 'Observations', 'Plausible Range', 'Niederdrenk and Notz, (2018)'}, 'location', 'southwest', 'fontsize', 16);
hLegend.Box = 'off';




%________CO2 PLOTS________%
for month_ind = [3,9]   

    h = subplot(1,2,2);

    if month_ind == 3
        ind = 0.2;
    else
        ind = 0.55; 
    end
    
    % OBS CO2
    if month_ind == 9
        boxplot([-1.36, -4.1], 'positions', [ind-0.045-0.045], 'color', colorss1(2,:), 'plotstyle','compact'); hold on

        obs_sens_bp = boxplot(obs_sens_co2{month_ind}(:), 'positions', [ind-0.045], 'color', [0 0 0]+0.7, 'plotstyle', 'compact'); hold on
        hAx = gca;
        lines = hAx.Children;
        lines = lines(1); 
        bp_cmip6 = findobj(lines, 'tag', 'Box');
        bp_cmip6.YData(1,1) = quantile(obs_sens_co2{month_ind}(:), 0.25);
        bp_cmip6.YData(1,2) = quantile(obs_sens_co2{month_ind}(:), 0.75);

    else
        obs_sens_bp = boxplot(obs_sens_co2{month_ind}(:), 'positions', [ind-0.045], 'color', [0 0 0]+0.7, 'plotstyle', 'compact'); hold on
        hAx = gca;
        lines = hAx.Children;
        lines = lines(1); 
        bp_cmip6 = findobj(lines, 'tag', 'Box');
        bp_cmip6.YData(1,1) = quantile(obs_sens_co2{month_ind}(:), 1.0);
        bp_cmip6.YData(1,2) = quantile(obs_sens_co2{month_ind}(:), 0.0);
    end

    
    
    % CMIP6 CO2
    boxplot(cmip6_sens_co2{month_ind}, 'positions', [ind], 'color', [1 0 0], 'plotstyle','compact'); hold on    % CMIP6
    hAx = gca;
    lines = hAx.Children;
    lines = lines(1);
    bp_cmip6 = findobj(lines, 'tag', 'Box');
    bp_cmip6.YData(1,1) = quantile( cmip6_sens_co2{month_ind}, 0.25);
    bp_cmip6.YData(1,2) = quantile( cmip6_sens_co2{month_ind}, 0.75);
    

   
    
    % % OCCAA CO2
    % OCCAA_CO2 = boxplot(magicc_sens_OCCAA_co2{month_ind}(:), 'positions', [ind+0.045], 'color', [0    0.5000    0.5000], 'plotstyle','compact'); hold on
    % set(OCCAA_CO2(length(OCCAA_CO2),:),'Visible','off')     
    % hAx = gca;
    % lines = hAx.Children;
    % lines = lines(1);
    % bp_cmip6 = findobj(lines, 'tag', 'Box');
    % bp_cmip6.YData(1,1) = quantile( magicc_sens_OCCAA_co2{month_ind}(:), 0.25);
    % bp_cmip6.YData(1,2) = quantile( magicc_sens_OCCAA_co2{month_ind}(:), 0.75); 

    
    % OC CO2
    % OC_CO2 = boxplot(magicc_sens_OC_co2{month_ind}(:), 'positions', [ind+0.045], 'color', colorss1(2,:), 'plotstyle','compact'); hold on
    OC_CO2 = boxplot(magicc_sens_OC_co2{month_ind}(:), 'positions', [ind+0.045], 'color', colorss1(1,:), 'plotstyle','compact'); hold on
    set(OC_CO2(length(OC_CO2),:),'Visible','off')     
    hAx = gca;
    lines = hAx.Children;
    lines = lines(1);
    bp_cmip6 = findobj(lines, 'tag', 'Box');
    bp_cmip6.YData(1,1) = quantile( magicc_sens_OC_co2{month_ind}(:), 0.25);
    bp_cmip6.YData(1,2) = quantile( magicc_sens_OC_co2{month_ind}(:), 0.75);
    
    
    
    
    % Set width of boxplots
    a = get(get(gca,'children'),'children');  
    for ii = 1:length(a)
        t = get(a{ii},'tag'); 
        idx=strcmpi(t,'box'); 
        a_new = a{ii};
        boxes=a_new(idx);         
        set(boxes,'linewidth',20); 
    end
    
    
    % Add plot labels and edits
    title([num2str(years(index_1974_2014(1))), '-', num2str(years(index_1974_2014(end))), ' Sensitivity (dSIA/dCO2)'], 'fontsize', 22, 'interpreter', 'none')
    
    if ismember(month_ind, [1,5,9])
        ylabel('m^2/ t')
    end

    
end
h.Position(1) = h.Position(1) + 0.02;
set(gca, 'fontsize', 20)
    
xticks([0.2, 0.5])
xticklabels( [string({month_label{3}}), string({month_label{9}})] );

length_curve_2 = 0.75;
xlim([0 length_curve_2])    
ylim([-4.5 0])    
grid
h.Position(4) = h.Position(4) - 0.1;
h.Position(3) = h.Position(3) - 0.03;
h.Position(1) = h.Position(1) - 0.07;


hLegend = legend(findall(gca,'Tag','Box'), {'This Study', 'CMIP6', 'Observations', 'Plausible Range'}, 'location', 'southwest', 'fontsize', 16);
hLegend.Box = 'off';


% Save Normalised
set(gcf, 'PaperOrientation', 'landscape')
set(gcf,'PaperSize',[51 30]);

% Save Normalised
temp=['SENS_GMST_CO2_2', '.pdf']; 
saveas(gca,temp); 



