%% Carbon Budget Calculations: Section 3.3

%% ____ Loading RCMIP cumulative CO2 data ____ %% 
% Extract MAGICC Ensemble Members for Emulation

M = readtable('rcmip-emissions-annual-means-v5-1-0.csv');

% Convert data from table to cell to matrix
carbon_data = table2cell(M);


% _____ Emissions|CO2 VARIABLE ALONE: Extract CO2 data from SSP585, 245, 126 ______ % 

% Extract "Variable": Emissions|CO2
CO2_ind = cellfun(@string,carbon_data(:,4));
CO2_ind = find(CO2_ind == "Emissions|CO2");


% "Mip_Era": Extract CMIP6 (not CMIP5)
CO2_Mip_Era = cellfun(@string,carbon_data(:,6));
CO2_Mip_Era_ind = find(CO2_Mip_Era == "CMIP6");


% "Region": world
CO2_region = cellfun(@string,carbon_data(:,3));
CO2_region_ind = find(CO2_region == "World");


% Emissions|CO2 & CMIP6 
Mip_in_CO2_ind = ismember(CO2_ind, CO2_Mip_Era_ind);
Mip_in_CO2_ind = find(Mip_in_CO2_ind == 1);
Mip_in_CO2_ind = CO2_ind(Mip_in_CO2_ind);


% Emissions|CO2 & CMIP6 & World
region_CO2_ind = ismember(Mip_in_CO2_ind, CO2_region_ind);
region_CO2_ind = find(region_CO2_ind == 1);
region_CO2_ind = Mip_in_CO2_ind(region_CO2_ind);


% Final CO2 data index
final_data_index = region_CO2_ind;

% Extract the CO2 data "DATA"
CO2_data = cell2mat(carbon_data(final_data_index,8:end));
 


% Extract SSP string data "SCENARIO"
carbon_data_updated = carbon_data(final_data_index,:);
CO2_ensemble_ssp_data = cellfun(@string,carbon_data_updated(:,2));



% Extract SSP Scenarios
scenarios_to_extract = ["ssp585"; "ssp245"; "ssp126"];
ind = [];
CO2_ensembles = [];
CO2_ensembles_GT = [];
for i = 1:length(scenarios_to_extract)
    
%   Find index 
    ind_hold = (CO2_ensemble_ssp_data==scenarios_to_extract(i));
    ind = find(ind_hold==1);
    
    
%   Extract the desired ssp scenarios ensembles
    CO2_ensembles{i} = CO2_data(ind,:);
    CO2_ensembles_GT{i} = (CO2_data(ind,:) .* 1e-3) / 1;

end


% ________________ Interpolate: 2015:2500 (index: 266:751) (as data every 10 years) ________________% 
colorss1 = [1 0.8 0.8; 0.8 0.8 1; 0 0.6906 0.5];
colorss2 = [1 0 0; 0 0 1; 0 0.3906 0];
years = 1750:2500;
co2_emissions_int = [];

% Plot to test how it looks
close all
figure(1)
set(gcf, 'Units', 'Inches', 'Position', [.4 .4 19 10])
for i = 1:3 % SSP
    
    co2_emissions = CO2_ensembles_GT{i};
    co2_emissions_int{i} = fillmissing(co2_emissions, "linear");
    hh(1) = plot(years, co2_emissions_int{i}, 'color', colorss2(i,:), 'linewidth', 3); hold on
      
end
% ___ Plot Edits ____ % 
ylabel(' Gt CO2/yr ')
ax = gca;
ax.YAxis.Exponent = 0;
grid
set(gca,'FontSize', 22)
sgtitle('CO_2 Emissions', 'fontsize', 28)
clear allChildren
allChildren = get(gca, 'Children');
legend_labels = ["SSP5-8.5", "SSP2-4.5", "SSP1-2.6"];
legend(legend_labels)
set(gca,'FontSize', 22)





% ________ Calculate Cumulative CO2 from Emissions Data ________ 5

colorss2 = [1 0 0; 0 0 1; 0 0.3906 0];
cumsum_co2_emissions_int = [];
years = 1750:2500;

close all
figure(2)
set(gcf, 'Units', 'Inches', 'Position', [.4 .4 19 10])
for i = 1:3
    cumsum_co2_emissions_int{i} = cumsum(co2_emissions_int{i});
    hh(1) = plot(years, cumsum_co2_emissions_int{i}, 'color', colorss2(i,:), 'linewidth', 3); hold on
end
% ___ Plot Edits ___ % 
ylabel(' Cumulative CO2 (Gt CO2)  ')
xlabel(' Years ')
ax = gca;
ax.YAxis.Exponent = 0;
sgtitle('CO_2 Cumulative Emissions', 'fontsize', 28)
grid
set(gca,'FontSize', 22)
legend_labels = ["SSP5-8.5", "SSP2-4.5", "SSP1-2.6"];
legend(legend_labels, 'location', 'northwest')













%% _____________ Linear mode months remaining carbon budget calculations _____________ %% 
%  Section 3.4


% Septemeber sea ice as a function of CO2:
% Initialisation 
threshold = [5, 95];
percentiles_0 = [];
percentiles_100 = [];
SIA_at_max_co2 = [];
SIA_store = [];
all_means = [];
all_percentiles = [];
all_percentiles_store_1 = [];
all_percentiles_store_2 = [];
co2_store = [];
co2_store_array = [];
y_store = [];
magicc_sens_co2 = [];
magicc_sens_co2_all = [];
cb_1_ensembles = [];
magicc_coefficients_hold_all = [];
SIA = sia_mon_rearrange_final_MCMC{1}(:,9);
SIA = cell2mat(SIA);
percentiles = prctile(SIA, threshold);
colorss1 = [1 0.8 0.8; 0.8 0.8 1; 0 0.6906 0.5];
colorss2 = [1 0 0; 0 0 1; 0 0.3906 0];
index_1974_2014 = [130:165]; % 1979-2014
clear text


close all
figure(3)
set(gcf, 'Units', 'Inches', 'Position', [.4 .4 19 10])
for n = 1:3
   
    h = subplot(1,2,1);
   


    % Extract emulated MAGICC SIA     
    SIA_mean = mean(cell2mat(sia_mon_rearrange_final_OC{n}(:,9)));
    SIA = cell2mat(sia_mon_rearrange_final_OC{n}(:,9));
    
    % Extract MAGICC cumulative CO2 from 1850-2300
    cum_co2_one = cumsum_co2_emissions_int{n}(101:551);
    cum_co2 = repmat(cum_co2_one, [600*12, 1]);
        
    
    
    
    %__ Plot all of September SIA as a function of Cumulative CO2____%
    stdY = prctile(SIA, threshold);
    
    
    % Years to plot the calibrated timeseries over
    x = cum_co2_one;   
        
    % Create curves: 
    curve1 = stdY(2,:);
    curve2 = stdY(1,:);

    
    % Shading std area: 
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];


    % Fill: 
    if n == 1
        c = [0 0 0]+0.55;
    elseif n == 2
        c = [0 0 0]+0.9;
    elseif n == 3
        c = [0 0 0]+0.9;
    end
    f = fill(x2, inBetween, c, 'edgecolor','none');
    set(f,'facealpha',.2)
    hold on

    
    % Extract only the linear period
    x = cum_co2;                % CO2
    y = SIA;                    % SIA
    x(y<1) = nan;
    y(y<1) = nan; 
    
    pi_SIA = mean(SIA(:,1:51), 2);      % Find pre-industrial to calculate 10% below this
    constant_co2 = max(max(cum_co2));   % Find when CO2 becomes constant in each scenario if this happens
    
    x(y > pi_SIA-(pi_SIA.*0.1)) = nan;  % Remove values that lie outside of the linear period (defined in Section 5.2.1)
    y(y > pi_SIA-(pi_SIA.*0.1)) = nan;
    
    % Only plot 5th-95th percentiles
    y(x == constant_co2) = nan;         % Remove values that are constant in CO2
    x(x == constant_co2) = nan;
     

    % remove SIA in SSP1-2.6 when its at a minumum as SIA recovers in this scenario effecting the linearity of the SIA/ CO2 relationship
    ind_ssp126 = [];
    if n == 3 || n == 2 
        min_SIA = min(y, [], 2);
        for ii = 1:length(SIA)
            ind_ssp126_hold = find(y(ii,:) == min_SIA(ii));
            y(ii,ind_ssp126_hold:end) = nan;
            x(ii,ind_ssp126_hold:end) = nan;
        end
    end
    
    x_linear = x;
    y_linear = y;

    % Only plot 5th-95th percentiles for plotting
    x(y > stdY(2,:)) = nan;
    x(y < stdY(1,:)) = nan;
    
    y(y > stdY(2,:)) = nan;
    y(y < stdY(1,:)) = nan;
    
    % Scatter only the linear period
    scatter(x(:), y(:), 'markerfacecolor', colorss1(1,:), 'markeredgecolor', colorss1(1,:)); hold on  
    
    
  
    
    
    
    % ____________ SENSITIVITIES ____________ %
    % Calculate the Sensitivity over 1979-2014
   
    clear coefficients
    magicc_sens_hold = [];
    for ii = 1:length(SIA)
        SIA_emul = SIA(ii,:);
        SIA_emul = SIA_emul(index_1974_2014);

        GT_emul = cum_co2_one;
        GT_emul = GT_emul(index_1974_2014);

        coefficients = polyfit(GT_emul, SIA_emul, 1);
        magicc_sens_hold{ii} = coefficients(1);
    end  
    magicc_sens_co2{1} = cell2mat(magicc_sens_hold)';   
    
    
    % Calculate the Sensitivity over the linear period
    clear coefficients
    magicc_sens_hold_all = [];
    for ii = 1:length(SIA)
        y_hold = y_linear(ii,:);
        x_hold = x_linear(ii,:);
        
        x_hold(isnan(x_hold) == 1) = [];
        y_hold(isnan(y_hold) == 1) = [];
        
        SIA_emul = y_hold;
        co2_emul = x_hold;
        
        if length(SIA_emul) < 2
            SIA_emul = nan;
            co2_emul = nan;
        end

        coefficients = polyfit(co2_emul, SIA_emul, 1);
        if coefficients == 0
            coefficients = [nan, nan];
        end
        magicc_sens_hold_all{ii} = coefficients(1);
        magicc_coefficients_hold_all{ii,n} = coefficients;
        
    end  
    magicc_sens_co2_all{n} = cell2mat(magicc_sens_hold_all)'; 
    
    
    
    
    
    % _____________________ Extract CO2 when SIA falls below ice-free conditions (<1million km2) to end the linear period _______________________________%
    x = cum_co2;
    y = SIA;

    
    % Calculate the Carbon Budget: OC
    for ii = 1:length(SIA)
        
        y_hold = y(ii,:);
        x_hold = x(ii,:);
        
        SIA_emul = y_hold;
        co2_emul = x_hold;
       
                 
        % Find CO2 cumulative emission when SIA falls below ice-free conditions & if it doesn't then assign a NaN
        ind_1_SIA = find(SIA_emul <= 1);
        if isempty(ind_1_SIA) ~= 1
            ind_1_co2 = co2_emul(ind_1_SIA(1));
        elseif isnan(y_linear(ii,1)) == 1
            ind_1_co2 = nan;
        else
            x_hold = x_linear(ii,:);
            y_hold = y_linear(ii,:);
           
            x_hold(isnan(x_hold) == 1) = [];
            y_hold(isnan(y_hold) == 1) = [];

            ind_1_co2 = interp1(y_hold, x_hold, 1,'linear','extrap');   % Extrapolate if SIA doesn't reach 1
        end
        
        cb_1_ensembles{ii,n} = (ind_1_co2);
    end
    

    
 
    
    % ____ Plot Edits ____ % 
    yline(1)    % to show ice-free conditions
    ax = gca;
    ax.XAxis.Exponent = 0;
    xlabel(' Cumulative CO_2 (Gt CO_2)  ')
    set(gca,'FontSize', 22)
    ylim([0 10])
    
end
text(150, 1.15, 'Ice Free', 'color', [0 0 0]+0.3, 'FontSize',15)
 

% Plot means from OC emulator to they are visible
for n = [2,3,1]
    % Extract emulated MAGICC SIA     
    SIA_mean = median(cell2mat(sia_mon_rearrange_final_OC{n}(:,9)));
    
    % Extract MAGICC cumulative CO2 from 1850-2300
    cum_co2_one = cumsum_co2_emissions_int{n}(101:551)';
    
    
    % Plot mean
    plt = plot(cum_co2_one, SIA_mean, 'color', colorss2(n,:), 'LineWidth', 2); hold on

end



% Plot Observations 
for m = 3:6
    i = 9;
    obs = monthly_obs_SIA{m}{i};
    cum_co2_one = cumsum_co2_emissions_int{1}(101:101+170-1);
    obs_plot = plot(cum_co2_one, obs, 'color', [0.3 0.3 0.3], 'LineWidth', 0.5); hold on
end



% ____ More Plot Edits ____ % 
xlim([cum_co2_one(1), 6000])
hold on
ylim([0, 10])
hlabel = text(-600, 1.8, ' September SIA (million km^2)  ', 'color', [0 0 0], 'FontSize', 27);
set(hlabel,'Rotation',90)




%  ________ Calculate the remaining carbon budget from CMIP6 data ________ % 
cb_1_ensembles_cmip6_months = [];
cb_1_ensembles_cmip6 = [];
magicc_sens_hold_all_cmip6 = []; 
for n = 1:3
    for i = 1:12


        % CMIP6 SIA
        SIA_CMIP6 = cell2mat(updated_hist_sia_annual_curve_all_models(:,n));
        SIA_CMIP6 = SIA_CMIP6(:,i);
        SIA_CMIP6 = reshape(SIA_CMIP6,[], 12);
        SIA_CMIP6 = SIA_CMIP6';
        SIA_CMIP6_hold = [];
        for j = 1:12
            SIA_CMIP6_hold = cat(1, SIA_CMIP6_hold, movmean(SIA_CMIP6(j,:), 20));
        end
        SIA_CMIP6 = SIA_CMIP6_hold;

        % CMIP6 CO2
        CO2_CMIP6 = cumsum_co2_emissions_int{n}(101:351);
        CO2_CMIP6 = repmat(CO2_CMIP6, [12, 1]);

        

            
        %____plot CMIP6 data_____%
        % Extract only the linear period
        x = CO2_CMIP6;
        y = SIA_CMIP6;
        x(y<1) = nan;
        y(y<1) = nan; 

        pi_SIA = mean(SIA_CMIP6(:,1:51), 2);
        constant_co2 = max(max(CO2_CMIP6));
        
        x(y > pi_SIA-(pi_SIA.*0.1)) = nan;
        y(y > pi_SIA-(pi_SIA.*0.1)) = nan;
        
        % Only plot 5th-95th percentiles
        y(x == constant_co2) = nan;
        x(x == constant_co2) = nan;
        
        x_linear = x;
        y_linear = y;




        % _____Calculate the Sensitivity: OVER LINEAR PERIOD______ % 
        clear coefficients
        for j = 1:12
            y_hold = y_linear(j,:);
            x_hold = x_linear(j,:);
            
            x_hold(isnan(x_hold) == 1) = [];
            y_hold(isnan(y_hold) == 1) = [];
            
            SIA_emul = y_hold;
            co2_emul = x_hold;
            
            if length(SIA_emul) < 2
                SIA_emul = nan;
                co2_emul = nan;
            end
    
            coefficients = polyfit(co2_emul, SIA_emul, 1);
            if coefficients == 0
                coefficients = [nan, nan];
            end
            magicc_sens_hold_all_cmip6{j,n} = coefficients(1);            
        end  

        
        
                
        % ____Find cumulative emission when SIA falls below ice-free conditions & if it doesn't then assign a NaN_______ % 
        for j = 1:12
            ind_1_SIA = find(SIA_CMIP6(j,:) <= 1);
            if isempty(ind_1_SIA) ~= 1
                ind_1_co2 = CO2_CMIP6(j,ind_1_SIA(1));
            elseif isnan(y_linear(j,1)) == 1
                ind_1_co2 = nan;
            else
                x_hold = x_linear(j,:);
                y_hold = y_linear(j,:);
    
                x_hold(isnan(x_hold) == 1) = [];
                y_hold(isnan(y_hold) == 1) = [];
    
                ind_1_co2 = interp1(y_hold, x_hold, 1,'linear','extrap');   % Extrapolate if SIA doesn't reach 1
            end
    
            cb_1_ensembles_cmip6{j,i} = (ind_1_co2);
            
        end
    end
    cb_1_ensembles_cmip6_months = cat(1, cb_1_ensembles_cmip6_months, cb_1_ensembles_cmip6); 
end
grid
cb_1_ensembles_cmip6 = cb_1_ensembles_cmip6_months;






% ____ Get info ready for boxplots next to plot above ____ %
% Carbon Budget since 1750
cb_1_ensembles2 = cell2mat(cb_1_ensembles(:));
cb_1_ensembles = cell2mat(cb_1_ensembles);

% Calculate the remaining carbon budget from 2023 based on values above
years = 1850:2023;
ind = find(years == 2023);
co2_2023 = co2_store(:,ind)';
CB_store = cb_1_ensembles - co2_2023;       % Carbon Budget 2023 EMULATOR
co2_2023_cmip6 = CO2_CMIP6(:,ind)';
CB_store_cmip6_2023 = [];
for i = 1:12
    CB_store_cmip6_2023{i} = cell2mat(cb_1_ensembles_cmip6(:,i)) - (co2_2023_cmip6(1));    % Carbon Budget 2023 CMIP6  
end


% Sensitivities Conversions: Convert carbon budget to the amount of SIA lost per ton of warming in m2
% OC emulator
sensitivities_co2_1979_2014 = (magicc_sens_co2{1} .* 1e12) ./ 1e9;
magicc_sens_co2_all_2 = cell2mat(magicc_sens_co2_all(:));
sensitivities_co2_all = (magicc_sens_co2_all_2 .* 1e12) ./ 1e9;

% CMIP6 data
magicc_sens_co2_all_cmip6 = cell2mat(magicc_sens_hold_all_cmip6(:));
sensitivities_co2_all_cmip6 = (magicc_sens_co2_all_cmip6 .* 1e12) ./ 1e9;
sensitivities_co2_each_ssp = [];
for n = 1:3
    sensitivities_co2_each_ssp{n} = (magicc_sens_co2_all{n} .* 1e12) ./ 1e9;
end
summer_sens = sensitivities_co2_all;




% _______ Add Boxplots: September _____ % 
% Calculate the observational Sensitivity (SIA lost per emitted tonne of co2 in m2/ton)
obs_sens_co2 = [];
obs_sens_hold = [];
clear coefficients 
index_1974_2014 = [130:165];
movmean_ind = 20;
month_ind = 9;
for m = [3,4,5]
    CO2_MAGICC = cumsum_co2_emissions_int{1}(101:271);
    CO2_MAGICC = CO2_MAGICC(index_1974_2014);

    SIA_obs = monthly_obs_SIA{m}{month_ind};
    SIA_obs = SIA_obs(index_1974_2014);

    ind = isnan(SIA_obs);
    if isempty(ind) ~= 1
        SIA_obs = SIA_obs(ind==0);
        CO2_MAGICC = CO2_MAGICC(ind==0);
    else
        SIA_obs = SIA_obs;
        CO2_MAGICC = CO2_MAGICC;
    end

    coefficients = polyfit(CO2_MAGICC, SIA_obs, 1);
    obs_sens_hold{m} = coefficients(1);
end
obs_sens_co2 = cell2mat(obs_sens_hold);
sensitivities_co2_obs = (obs_sens_co2 .* 1e12) ./ 1e9;





% Box plot of CMIP6 Sensitivities 
clear coefficients
cmip6_sens = [];
cmip6_sens_hold = [];
index_1974_2014 = [130:165];
movmean_ind = 20;
month_ind = 9;
for n = 1
    for j = 1:12
        SIA_cmip6 = updated_hist_sia_annual_curve_all_models{j,n}(:,month_ind);
        SIA_cmip6 = movmean(SIA_cmip6, movmean_ind);
        SIA_cmip6 = SIA_cmip6(index_1974_2014);

        % CMIP6 CO2
        CO2_CMIP6 = cumsum_co2_emissions_int{n}(101:351);
        CO2_CMIP6 = movmean(CO2_CMIP6, movmean_ind);
        CO2_CMIP6 = CO2_CMIP6(index_1974_2014);

        coefficients = polyfit(CO2_CMIP6, SIA_cmip6, 1);
        cmip6_sens_hold{j,n} = coefficients(1); 
    end
end
cmip6_sens = cell2mat(cmip6_sens_hold);
cmip6_sens = (cmip6_sens * 1e12) / 1e9;





% BOX 1: Sensitivity over 1979-2014 period
axes('Position',[.55 .57 .13 .33])
box on

% Plot CMIP6 1979-2014
hold on
sens_79_14 = boxplot(cmip6_sens(:), 'positions', 0.9, 'color', [1 0 0], 'plotstyle', 'compact');
a = get(get(gca,'children'),'children');   % Get the handles of all the objects
t = get(a,'tag');   % List the names of all the objects 
idx=strcmpi(t,'box');  % Find Box objects
boxes=a(idx);          % Get the children you need
set(boxes,'linewidth',10); % Set width


Median_dot_lines = findobj(sens_79_14, 'type', 'line', 'Tag', 'MedianInner');     % Make median bigger and bolder
set(Median_dot_lines, 'MarkerFaceColor', 'k');
set(Median_dot_lines, 'Color', 'k');
set(Median_dot_lines, 'MarkerSize', 12);

Median_dot_lines = findobj(sens_79_14, 'type', 'line', 'Tag', 'MedianOuter');
Median_dot_lines(1).Visible = 'off';


% ____ Plot observed sensitivity 1979-2014 ____ % 
hold on
sens_79_14 = boxplot(sensitivities_co2_obs, 'positions', 1.2, 'color', [0 0 0]+0.5, 'plotstyle', 'compact');
a = get(get(gca,'children'),'children');   % Get the handles of all the objects
a = a{1};
t = get(a,'tag');   % List the names of all the objects 
idx=strcmpi(t,'box');  % Find Box objects
boxes=a(idx);          % Get the children you need
set(boxes,'linewidth',10); % Set width


Median_dot_lines = findobj(sens_79_14, 'type', 'line', 'Tag', 'MedianInner');     % Make median bigger and bolder
set(Median_dot_lines, 'MarkerFaceColor', 'k');
set(Median_dot_lines, 'Color', 'k');
set(Median_dot_lines, 'MarkerSize', 12);

Median_dot_lines = findobj(sens_79_14, 'type', 'line', 'Tag', 'MedianOuter');
Median_dot_lines(1).Visible = 'off';



% ____ Plot OC emulator sensitivity 1979-2014 ____ % 
sens_79_14_em = boxplot(sensitivities_co2_1979_2014, 'positions', 1.5, 'color', colorss1(1,:), 'plotstyle', 'compact');
hAx = gca;
lines = hAx.Children;
lines = lines(1);
bp_cmip6 = findobj(lines, 'tag', 'Box');
bp_cmip6.YData(1,1) = quantile(sensitivities_co2_1979_2014, 0.75);
bp_cmip6.YData(1,2) = quantile(sensitivities_co2_1979_2014, 0.25);


Median_dot_lines = findobj(sens_79_14_em, 'type', 'line', 'Tag', 'MedianInner');     % Make median bigger and bolder
set(Median_dot_lines, 'MarkerFaceColor', 'k');
set(Median_dot_lines, 'Color', 'k');
set(Median_dot_lines, 'MarkerSize', 12);

Median_dot_lines = findobj(sens_79_14_em, 'type', 'line', 'Tag', 'MedianOuter');
Median_dot_lines(1).Visible = 'off';

a = get(get(gca,'children'),'children');   % Get the handles of all the objects
a = a{1};
t = get(a,'tag');   % List the names of all the objects 
idx=strcmpi(t,'box');  % Find Box objects
boxes=a(idx);          % Get the children you need
set(boxes,'linewidth',10); % Set width





% Box 2: Plot Sensitivities over linear period from OC emulator and other sstudies
axes('Position',[.77 .52 .13 .379])
box on
sens_linear = boxplot(sensitivities_co2_all, 'positions', 1.2, 'color', colorss1(1,:), 'plotstyle', 'compact'); hold on
hAx = gca;
lines = hAx.Children;
bp_cmip6 = findobj(lines, 'tag', 'Box');
bp_cmip6.YData(1,1) = quantile(sensitivities_co2_all, 0.75);
bp_cmip6.YData(1,2) = quantile(sensitivities_co2_all, 0.25);

% Extracting data from above plot
lines = findobj(gcf, 'type', 'line', 'Tag', 'MedianInner');
Median_all_sensitivities = lines.YData;
upper_83_all_sensitivities = bp_cmip6.YData(1,2);
lower_17_all_sensitivities = bp_cmip6.YData(1,1);


Median_dot_lines = findobj(sens_linear, 'type', 'line', 'Tag', 'MedianInner');     % Make median bigger and bolder
set(Median_dot_lines, 'MarkerFaceColor', 'k');
set(Median_dot_lines, 'Color', 'k');
set(Median_dot_lines, 'MarkerSize', 12);

Median_dot_lines = findobj(sens_linear, 'type', 'line', 'Tag', 'MedianOuter');
Median_dot_lines(1).Visible = 'off';


a = get(get(gca,'children'),'children');   % Get the handles of all the objects
t = get(a,'tag');   % List the names of all the objects 
idx=strcmpi(t,'box');  % Find Box objects
boxes=a(idx);          % Get the children you need
set(boxes,'linewidth',10); % Set width


% Plot Edits
grid
set(gca, 'xtick', [0.9, 1.2, 1.5], 'XTickLabels', compose(["CMIP6", "Observations", "This Study"]), 'fontsize', 14)
xtickangle(30)
title([{ ' September Sensitivity ' ; ' 1979-2014 (dSIA/dC02) ' }], 'fontsize', 12, 'interpreter', 'none')
ylabel(' m^2/ t')
xlim([0.5 2])
ylim([-4.7 0])


% __ SIMIP Community, (2020) paper __ %
SIMIP_Community_co2_sens = [-1.36, -2.73, -4.1];  % SIMIP Community (taken from SIMIP Community, 2020)
sens_linear_SIMIP_community = boxplot(SIMIP_Community_co2_sens, 'positions', 1, 'color', colorss1(2,:), 'plotstyle', 'compact');
a = get(get(gca,'children'),'children');   % Get the handles of all the objects
a = a{1};
t = get(a,'tag');   % List the names of all the objects 
idx=strcmpi(t,'box');  % Find Box objects
boxes=a(idx);          % Get the children you need
set(boxes,'linewidth',10); % Set width
Median_dot_lines = findobj(sens_linear_SIMIP_community, 'type', 'line', 'Tag', 'MedianInner');     % Make median bigger and bolder
set(Median_dot_lines, 'MarkerFaceColor', 'k');
set(Median_dot_lines, 'Color', 'k');
set(Median_dot_lines, 'MarkerSize', 12);
Median_dot_lines = findobj(sens_linear_SIMIP_community, 'type', 'line', 'Tag', 'MedianOuter');
Median_dot_lines(1).Visible = 'off';




% __Sensitivity over linear period CMIP6 __ %
sens_linear_cmip6 = boxplot(sensitivities_co2_all_cmip6, 'positions', 0.8, 'color', [1 0 0], 'plotstyle', 'compact');
hAx = gca;
lines = hAx.Children;
clear bp_cmip6
bp_cmip6 = findobj(lines, 'tag', 'Box');
bp_cmip6(1).YData(1,1) = quantile(sensitivities_co2_all_cmip6, 0.75);
bp_cmip6(1).YData(1,2) = quantile(sensitivities_co2_all_cmip6, 0.25);
set(sens_linear_cmip6(length(sens_linear_cmip6),:),'Visible','off')     % Remove the outliers
lines = findobj(gcf, 'type', 'line', 'Tag', 'MedianInner');
Median_all_sensitivities = lines.YData;
upper_83_all_sensitivities = bp_cmip6(1).YData(1,2);
lower_17_all_sensitivities = bp_cmip6(1).YData(1,1);
Median_dot_lines = findobj(sens_linear_cmip6, 'type', 'line', 'Tag', 'MedianInner');     % Make median bigger and bolder
set(Median_dot_lines, 'MarkerFaceColor', 'k');
set(Median_dot_lines, 'Color', 'k');
set(Median_dot_lines, 'MarkerSize', 12);
Median_dot_lines = findobj(sens_linear_cmip6, 'type', 'line', 'Tag', 'MedianOuter');
Median_dot_lines(1).Visible = 'off';


a = get(get(gca,'children'),'children');   % Get the handles of all the objects
a = a{1};
t = get(a,'tag');   % List the names of all the objects 
idx=strcmpi(t,'box');  % Find Box objects
boxes=a(idx);          % Get the children you need
set(boxes,'linewidth',10); % Set width


% Plot Edits
grid
set(gca, 'xtick', [0.8, 1, 1.2], 'XTickLabels', compose(["CMIP6", "SIMIP Community, (2020)", "This Study"]), 'fontsize', 14)
xtickangle(35)
title([' September Sensitivity (dSIA/dC02) '], 'fontsize', 12, 'interpreter', 'none')
title([{ ' September Sensitivity ' ; ' linear period (dSIA/dC02) ' }], 'fontsize', 12, 'interpreter', 'none')
ylabel(' m^2/ t')
xlim([0.5 1.5])
ylim([-4.7 0])




% % Save
% set(gcf, 'PaperOrientation', 'landscape')
% set(gcf,'PaperSize',[53 29.7000]);

% % Save
% temp=['Fig5.2', '.pdf']; 
% saveas(gca,temp);















%% __________ Calculate remaining cabron budget from 2023 and 1850 from all linear mode months __________ % 


close all
figure(3)
set(gcf, 'Units', 'Inches', 'Position', [.4 .4 19 10])
clc
axes('Position',[.3 .2 .5 .4])

cb_1_ensembles_final = [];
cb_1_ensembles_final_months = [];
cb_1_ensembles_final_months_2023 = [];
cb_1_ensembles = [];
index_months = [7:12];
for i = index_months
    co2_store = [];
    co2_2023 = [];
    cb_1_ensembles_final = [];
    for n = 1:3
            
        % Extract emulated MAGICC SIA     
        SIA_mean = mean(cell2mat(sia_mon_rearrange_final_MCMC{n}(:,i)));
        SIA = cell2mat(sia_mon_rearrange_final_MCMC{n}(:,i));
        
        % Extract MAGICC cumulative CO2 from 1850-2300
        cum_co2_one = cumsum_co2_emissions_int{n}(101:551);
        cum_co2 = repmat(cum_co2_one, [600*12, 1]);


        % Calculate the carbon budget based on values above
        co2_store = cat(1, co2_store, cum_co2_one);

           

        % Extract only the linear part
        x = cum_co2;
        y = SIA;
        x(y<1) = nan;
        y(y<1) = nan; 
        
        pi_SIA = mean(SIA(:,1:51), 2);
        constant_co2 = max(max(cum_co2));
        
        x(y > pi_SIA-(pi_SIA.*0.1)) = nan;
        y(y > pi_SIA-(pi_SIA.*0.1)) = nan;
        
        % Only plot 5th-95th percentiles
        y(x == constant_co2) = nan;
        x(x == constant_co2) = nan;

        ind_ssp126 = [];              % remove SIA in SSP1-2.6 when its at a minumum as SIA recovers in this scenario effecting the linearity of the SIA/ CO2 relationship
        if n == 3 || n == 2 
            min_SIA = min(SIA, [], 2);
            for ii = 1:length(SIA)
                ind_ssp126_hold = find(SIA(ii,:) == min_SIA(ii));
                if ind_ssp126_hold ~= 451
                    y(ii,:) = nan;
                    x(ii,:) = nan;
                end
            end
        end
        
        x_linear = x;
        y_linear = y;

    
        
        
        for ii = 1:length(SIA)     % Calculate the Sensitivity: Observationally constrained SIA using CMIP6 calibrated AA   
            
            SIA_emul = SIA(ii,:);
            co2_emul = cum_co2(ii,:);
                                
            ind_1_SIA = find(SIA_emul <= 1);         % Find CO2 cumulative emission when SIA reaches 1 & if it doesn't then assign a NaN
            if isempty(ind_1_SIA) ~= 1
                ind_1_co2 = co2_emul(ind_1_SIA(1));
            elseif isnan(y_linear(ii,1)) == 1
                ind_1_co2 = nan;
            else

                y_hold = y_linear(ii,:);
                x_hold = x_linear(ii,:);

                y_hold(isnan(y_hold) == 1) = [];
                x_hold(isnan(x_hold) == 1) = [];
    
                ind_1_co2 = interp1(y_hold, x_hold, 1,'linear','extrap');   % Extrapolate if SIA doesn't reach 1
            end
            
            cb_1_ensembles{ii,i} = (ind_1_co2);
        end
        cb_1_ensembles_final = cat(1, cb_1_ensembles_final, cb_1_ensembles); 
    end
    cb_1_ensembles_final_months{i} = cb_1_ensembles_final(:,i); 


    years = 1850:2023;
    ind = find(years == 2023);
    co2_2023 = co2_store(:,ind);         % Cumulative CO2 emissions in 2023 EMULATOR
    co2_2023 = mean(co2_2023);
    CB_for_plotting = cell2mat(cb_1_ensembles_final(:,i)) - co2_2023;      % Carbon Budget 2023 EMULATOR
    cb_1_ensembles_final_months_2023{i} = CB_for_plotting;

end
linear_mode_vert_plot = CB_for_plotting;












% ________ Plot CARBON BUDGET Fig 5.3 Chapter 5, Section 5.3.1: linear mode months as horizontal boxplots _______ %
clc
close all
figure(3)
set(gcf, 'Units', 'Inches', 'Position', [.4 .4 19 10])
colorss1 = [1 0.8 0.8; 0.8 0.8 1; 0 0.6906 0.5];
yticks_hold = [];
yticks_hold_2 = [];
h = tiledlayout(1,1);
ax1 = axes(h);
ind_month = [7:11];      % only plot the linear mode months (July-December 7:12)
for i = ind_month  % Keep to double check 1850 matches the 2023 data

    % % if i == 10
    %     % ____ Emulator: CB 1850 ____ %
    %     CBR = boxplot(ax1, cell2mat(cb_1_ensembles_final_months{i})-cumsum_co2_emissions_int{1}(101), 'positions', i+0.3, 'color', colorss1(1,:), 'plotstyle', 'compact', 'orientation', 'horizontal'); hold on % CB from 1750
    %     set(CBR(length(CBR),:),'Visible','off')     % Remove the outliers
    % 
    %     hAx = gca;
    %     lines = hAx.Children;
    %     lines = lines(1);
    %     bp_cmip6 = findobj(lines, 'tag', 'Box');
    %     bp_cmip6.XData(1,1) = quantile(cell2mat(cb_1_ensembles_final_months{i})-cumsum_co2_emissions_int{1}(101), 0.83);
    %     bp_cmip6.XData(1,2) = quantile(cell2mat(cb_1_ensembles_final_months{i})-cumsum_co2_emissions_int{1}(101), 0.17);
    % 
    % 
    %     xline(quantile(cell2mat(cb_1_ensembles_final_months{i})-cumsum_co2_emissions_int{1}(101), 0.83))
    %     xline(quantile(cell2mat(cb_1_ensembles_final_months{i})-cumsum_co2_emissions_int{1}(101), 0.17))
    % % end
    
    % % ____ CMIP6: CB 1750 ____ %
    % cmip6_CB = boxplot(ax1, cell2mat(cb_1_ensembles_cmip6(:,i))-cumsum_co2_emissions_int{1}(101), 'positions', i-0, 'color', [1 0 0], 'plotstyle', 'compact', 'orientation', 'horizontal'); hold on
    % set(cmip6_CB(length(cmip6_CB),:),'Visible','off')     % Remove the outliers
    % 
    % hAx = gca;
    % lines = hAx.Children;
    % lines = lines(1);
    % bp_cmip6 = findobj(lines, 'tag', 'Box');
    % bp_cmip6.XData(1,1) = quantile(cell2mat(cb_1_ensembles_cmip6(:,i))-cumsum_co2_emissions_int{1}(101), 0.83);
    % bp_cmip6.XData(1,2) = quantile(cell2mat(cb_1_ensembles_cmip6(:,i))-cumsum_co2_emissions_int{1}(101), 0.17);

    % xline(quantile(cell2mat(cb_1_ensembles_cmip6(:,i))-cumsum_co2_emissions_int{1}(101), 0.83))
    % xline(quantile(cell2mat(cb_1_ensembles_cmip6(:,i))-cumsum_co2_emissions_int{1}(101), 0.17))
 

    % ax1 axis 
    xlabel(' Carbon Budget 1850 Cumulative CO_2 (Gt CO_2)  ')

end

ax2 = axes(h);
for i = ind_month           % Plot the remaining carbon budget from 2023 in linear mode months only, but change axis from 1850 to match the 2023 carbon budget

    % ____ Emulator: CB 2023 ____ %
    CB23 = boxplot(ax2, cb_1_ensembles_final_months_2023{i}, 'positions', i+0.15, 'color', colorss1(1,:), 'plotstyle', 'compact', 'orientation', 'horizontal'); hold on
    set(CB23(length(CB23),:),'Visible','off')         % Remove the outliers

    hAx = gca;
    lines = hAx.Children;
    lines = lines(1);
    bp_cmip6 = findobj(lines, 'tag', 'Box');
    bp_cmip6.XData(1,1) = quantile(cb_1_ensembles_final_months_2023{i}, 0.25);
    bp_cmip6.XData(1,2) = quantile(cb_1_ensembles_final_months_2023{i}, 0.75);


    % ____ CMIP6: CB 2023 ____ %
    cmip6_CB = boxplot(ax2, CB_store_cmip6_2023{i}, 'positions', i-0.15, 'color', [1 0 0], 'plotstyle', 'compact', 'orientation', 'horizontal'); 
    set(cmip6_CB(length(cmip6_CB),:),'Visible','off')         % Remove the outliers

    hAx = gca;
    lines = hAx.Children;
    lines = lines(1);
    bp_cmip6 = findobj(lines, 'tag', 'Box');
    bp_cmip6.XData(1,1) = quantile(CB_store_cmip6_2023{i}, 0.25);
    bp_cmip6.XData(1,2) = quantile(CB_store_cmip6_2023{i}, 0.75);

    
    a = get(get(gca,'children'),'children');   % Get the handles of all the objects
    for ii = 1:length(a)
        anew = a{ii};
        t = get(anew,'tag');   % List the names of all the objects 
        idx=strcmpi(t,'box');  % Find Box objects
        boxes=anew(idx);          % Get the children you need
        set(boxes,'linewidth',10); % Set width
    end
    
    
    % Plot Edits
    hold on
    yline((i-0.15) - 0.2750, 'color', [0 0 0]+0.7, 'linewidth', 1.5)
    set(ax2, 'Position', get(ax1, 'Position'), 'Color', 'none');
    ax2.XAxisLocation = 'top';   
    ax2.YAxisLocation = 'right';   
    ax1.XLim = [1846 6500];
    ax2.XLim = [-600, 4020];

    % ax2 axis 
    xlabel(' Carbon Budget 2023 Cumulative CO_2 (Gt CO_2)  ')
    yticks_hold = cat(2, yticks_hold, [i-0.15, i+0.15]);
    yticks_hold_2 = cat(2, yticks_hold_2, [i+0.0750]);

    grid
end
% Carbon Budget Shaded Area below 0 Gt CO2
ind_ipcc = -800;
x = ind_ipcc:0;
x2 = [x, fliplr(x)]; 
curve1 = repelem(15, length(x));
curve2 = repelem(-2, length(x));
inBetween = [curve1, fliplr(curve2)];
c = [0 0 0]+0.7;
f2 = fill(x2, inBetween, c, 'edgecolor','none'); 
set(f2,'facealpha',.3)
xline(0, 'color', [0 0 0], 'linewidth', 1.5)
set(ax1, 'fontsize', 14)
set(ax2, 'fontsize', 14)
ax2.YLim = [6.7 11.5];
ax1.YLim = [6.7 11.5];


% ____ axis limits ____ % 
ytick_labels_hold = repmat(["CMIP6", "This Study"], [1,length(ind_month)]);
set(ax1, 'ytick', yticks_hold,'YTickLabels', compose(ytick_labels_hold), 'fontsize', 14)
title('Carbon Budget to Prevent an Ice-free Summer Arctic Ocean', 'fontsize', 18)
h.Position = [0.23 0.3 0.45 0.6];
set(ax2, 'ytick', yticks_hold_2,'YTickLabels', compose(string(month_label(ind_month))), 'fontsize', 14)
grid
clc






% ___ Add IPCC remaining carbon budget to plot to compare (taken from Lamboll et al, 2023) ___ % 
axes('Position',[.148 .07 .608 .12])

% IPCC 1.5 degrees warming
hold on
correction_2020_2023 = 38*3;                % IPCC (from 1st Jan 2020) % they have been corrected for 2023 by assuming 38GtCO2/yr is lost (38*3)
IPCC_CB = boxplot(([80, 496]), 'positions', 1.0, 'color', [0 0 0]+0.5, 'plotstyle', 'compact', 'orientation', 'horizontal');
a = get(get(gca,'children'),'children');   % Get the handles of all the objects
t = get(a,'tag');   % List the names of all the objects 
idx=strcmpi(t,'box');  % Find Box objects
boxes=a(idx);          % Get the children you need
set(boxes,'linewidth',10); % Set width
Median_dot_lines = findobj(IPCC_CB, 'type', 'line', 'Tag', 'MedianOuter');     % Make median bigger and bolder
Median_dot_lines.XData = 275;
Median_dot_lines = findobj(IPCC_CB, 'type', 'line', 'Tag', 'MedianInner');
Median_dot_lines.XData = 275;


% IPCC 2 degrees warming
IPCC_CB_2 = boxplot(([900, 2300]-correction_2020_2023), 'positions', 1.2, 'color', [0 0 0]+0.5, 'plotstyle', 'compact', 'orientation', 'horizontal');
a = get(get(gca,'children'),'children');   % Get the handles of all the objects
a = a{1};
t = get(a,'tag');   % List the names of all the objects 
idx=strcmpi(t,'box');  % Find Box objects
boxes=a(idx);          % Get the children you need
set(boxes,'linewidth',10); % Set width
Median_dot_lines = findobj(IPCC_CB_2, 'type', 'line', 'Tag', 'MedianInner');
Median_dot_lines.XData = 1150; 
Median_dot_lines = findobj(IPCC_CB_2, 'type', 'line', 'Tag', 'MedianOuter');
Median_dot_lines.XData = 1150;
ind = ['IPCC (to limit warming to 2', char(176), 'C)'];
ind_2 = ['IPCC (to limit warming to 1.5', char(176), 'C)'];
set(gca, 'ytick', [1, 1.2],'YTickLabels', compose([string(ind_2), string(ind)]), 'fontsize', 14)
ind_new = [{'Comparison with IPCC Carbon Budget to Prevent Warming Thresholds'}];
text(-200, 1.4, ind_new, 'fontsize', 16, 'fontweight', 'bold')


% Add a Shaded Area below 0 Gt CO2 to indicate we have exceeded the carbon budget (if this occurs in each month)
x = ind_ipcc:0;
x2 = [x, fliplr(x)]; 
curve1 = repelem(15, length(x));
curve2 = repelem(-2, length(x));
inBetween = [curve1, fliplr(curve2)];
c = [0 0 0]+0.7;
f2 = fill(x2, inBetween, c, 'edgecolor','none'); 
set(f2,'facealpha',.3)
xline(0, 'color', [0 0 0], 'linewidth', 1.5)
xlabel(' Cumulative CO_2 (Gt CO_2)  ')
clear ax
ax = gca;
ax.XLim = [-600, 4020];
ax.YLim = [0.9 1.3];
set(ax, 'fontsize', 14)
set(ax, 'YAxisLocation', 'right');
grid


% % Save
% set(gcf, 'PaperOrientation', 'landscape')
% set(gcf,'PaperSize',[53 29.7000]);

% % Save
% temp=['Fig5.3', '.pdf']; 
% saveas(gca,temp); 














%% _______ Non-linear mode month remaining carbon budget _______ %%
% Fig 5.4 (Chapter 5)


threshold = [5, 95];
percentiles_0 = [];
percentiles_100 = [];
SIA_at_max_co2 = [];
co2_store = [];
SIA_store = [];
all_means = [];
all_percentiles = [];
all_percentiles_store_1 = [];
all_percentiles_store_2 = [];
co2_store = [];
co2_store_array = [];
y_store = [];
magicc_march_sens_co2 = [];
magicc_sens_co2_all = [];
cb_1_ensembles = [];
magicc_coefficients_hold_all = [];
magicc_sens_march_79_14 = []; 
colorss1 = [1 0.8 0.8; 0.8 0.8 1; 0 0.6906 0.5];
colorss2 = [1 0 0; 0 0 1; 0 0.3906 0];
index_1974_2014 = [130:165]; % 1979-2014
clear text
month_ind = 3;
close all
hh = figure(3);
set(gcf, 'Units', 'Inches', 'Position', [.4 .4 19 10])


for n = 1:3     % SSPs
   
    h = subplot(1,2,1);
   


    % Extract emulated MAGICC SIA     
    SIA_mean = mean(cell2mat(sia_mon_rearrange_final_OC{n}(:,month_ind)));
    SIA = cell2mat(sia_mon_rearrange_final_OC{n}(:,month_ind));
    
    % Extract MAGICC cumulative CO2 from 1850-2300
    cum_co2_one = cumsum_co2_emissions_int{n}(101:551);
    cum_co2 = repmat(cum_co2_one, [600*12, 1]);
    
    
    
    %__ Plot all of March SIA as a function of Cumulative CO2____%
    stdY = prctile(SIA, threshold);
    
    
    % Years to plot the calibrated timeseries over
    x = cum_co2_one;   
        
    % Create curves: 
    curve1 = stdY(2,:);
    curve2 = stdY(1,:);

    
    % Shading std area: 
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];


    % Fill: 
    if n == 1
        c = [0 0 0]+0.55;
    elseif n == 2
        c = [0 0 0]+0.9;
    elseif n == 3
        c = [0 0 0]+0.9;
    end
    f = fill(x2, inBetween, c, 'edgecolor','none');
    set(f,'facealpha',.2)
    hold on
    
    
    % _________ plot linear period _________ % 
    x = cum_co2;
    y = SIA;
    x(y<1) = nan;
    y(y<1) = nan; 
    
    pi_SIA = mean(SIA(:,1:51), 2);
    constant_co2 = max(max(cum_co2));
    
    x(y > pi_SIA-(pi_SIA.*0.1)) = nan;
    y(y > pi_SIA-(pi_SIA.*0.1)) = nan;
    
    % Only plot 5th-95th percentiles
    y(x == constant_co2) = nan;
    x(x == constant_co2) = nan; 
     

    % remove SIA in SSP1-2.6 when its at a minumum as SIA recovers in this scenario effecting the linearity of the SIA/ CO2 relationship
    ind_ssp126 = [];
    if n == 3 || n == 2 
        min_SIA = min(y, [], 2);
        for ii = 1:length(SIA)
            ind_ssp126_hold = find(y(ii,:) == min_SIA(ii));
            y(ii,ind_ssp126_hold:end) = nan;
            x(ii,ind_ssp126_hold:end) = nan;
        end
    end
    
    x_linear = x;
    y_linear = y;

    % Only plot 5th-95th percentiles for plotting
    x(y > stdY(2,:)) = nan;
    x(y < stdY(1,:)) = nan;
    
    y(y > stdY(2,:)) = nan;
    y(y < stdY(1,:)) = nan;
    
    % Scatter only the linear part
    scatter(x(:), y(:), 'markerfacecolor', colorss1(1,:), 'markeredgecolor', colorss1(1,:)); hold on  
    
    
  
    
    
    
    % _____________ Calculate the Sensitivities _______________%
    % 1979-2014
    clear coefficients
    magicc_sens_hold = [];
    for ii = 1:length(SIA)
        SIA_emul = SIA(ii,:);
        SIA_emul = SIA_emul(index_1974_2014);

        GT_emul = cum_co2_one;
        GT_emul = GT_emul(index_1974_2014);

        coefficients = polyfit(GT_emul, SIA_emul, 1);
        magicc_sens_hold{ii} = coefficients(1);
    end  
    if n == 1
        magicc_sens_march_79_14{n} = cell2mat(magicc_sens_hold)';
    end
    
    
    % linear period
    clear coefficients
    magicc_sens_hold_all = [];
    for ii = 1:length(SIA)
        y_hold = y_linear(ii,:);
        x_hold = x_linear(ii,:);
        
        x_hold(isnan(x_hold) == 1) = [];
        y_hold(isnan(y_hold) == 1) = [];
        
        SIA_emul = y_hold;
        co2_emul = x_hold;
        
        if length(SIA_emul) < 2
            SIA_emul = nan;
            co2_emul = nan;
        end

        coefficients = polyfit(co2_emul, SIA_emul, 1);
        if coefficients == 0
            coefficients = [nan, nan];
        end
        magicc_sens_hold_all{ii} = coefficients(1);
        magicc_coefficients_hold_all{ii,n} = coefficients;
        
    end  
    magicc_sens_co2_all{n} = cell2mat(magicc_sens_hold_all)'; 
    
    
    
    
    
    % _________ Extract CO2 emission SIA falls below ice-free conditions _________ %
    x = cum_co2;
    y = SIA;
    magicc_sens_hold_all = [];
    for ii = 1:length(SIA)
        
        y_hold = y(ii,:);
        x_hold = x(ii,:);
        
        SIA_emul = y_hold;
        co2_emul = x_hold;
       
        
        % Find CO2 cumulative emission when SIA falls below 1million km2 & if it doesn't then assign a NaN
        ind_1_SIA = find(SIA_emul <= 1);
        if isempty(ind_1_SIA) ~= 1
            ind_1_co2 = co2_emul(ind_1_SIA(1));
        elseif isnan(y_linear(ii,1)) == 1
            ind_1_co2 = nan;
        else
            x_hold = x_linear(ii,:);
            y_hold = y_linear(ii,:);
           
            x_hold(isnan(x_hold) == 1) = [];
            y_hold(isnan(y_hold) == 1) = [];

            ind_1_co2 = interp1(y_hold, x_hold, 1,'linear','extrap');   % Extrapolate if SIA doesn't reach 1
        end
        
        cb_1_ensembles{ii,n} = (ind_1_co2);
    end
    

    
 
    
    % Plot Edits
    yline(1)
    ax = gca;
    ax.XAxis.Exponent = 0;
    xlabel(' Cumulative CO_2 (Gt CO_2)  ')
    set(gca,'FontSize', 22)  
end
text(150, 1.25, 'Ice Free', 'color', [0 0 0]+0.3, 'FontSize',15)
 

% Plot medians
for n = [2,3,1]
    % Extract emulated MAGICC SIA     
    SIA_mean = median(cell2mat(sia_mon_rearrange_final_OC{n}(:,month_ind)));
    

    % Extract MAGICC cumulative CO2 from 1850-2300
    cum_co2_one = cumsum_co2_emissions_int{n}(101:551)';
    
    
    % Plot mean
    plt = plot(cum_co2_one, SIA_mean, 'color', colorss2(n,:), 'LineWidth', 2); hold on

    grid
end



% Plot Observations
hnew_obs = [];
for m = 3:6
    i = month_ind;
    obs = monthly_obs_SIA{m}{i};
    cum_co2_one = cumsum_co2_emissions_int{1}(101:101+170-1);
    obs_plot = plot(cum_co2_one, obs, 'color', [0.3 0.3 0.3], 'LineWidth', 0.5); hold on
end




% More Plot Edits
xlim([cum_co2_one(1), max(cum_co2_one)])
ylim([0, 18])
hlabel = text(-2500, 3.8, ' March SIA (million km^2)  ', 'color', [0 0 0], 'FontSize', 27);
set(hlabel,'Rotation',90)
xline( median(CB_at_b(:,1), 'omitnan'), 'LineWidth', 1 )



% CMIP6 data
cb_1_ensembles_cmip6_months = [];
cb_1_ensembles_cmip6 = [];
magicc_sens_hold_all_cmip6 = []; 
for n = 1:3
    for i = 1:12


        % CMIP6 SIA
        SIA_CMIP6 = cell2mat(updated_hist_sia_annual_curve_all_models(:,n));
        SIA_CMIP6 = SIA_CMIP6(:,i);
        SIA_CMIP6 = reshape(SIA_CMIP6,[], 12);
        SIA_CMIP6 = SIA_CMIP6';
        SIA_CMIP6_hold = [];
        for j = 1:12
            SIA_CMIP6_hold = cat(1, SIA_CMIP6_hold, movmean(SIA_CMIP6(j,:), 20));
        end
        SIA_CMIP6 = SIA_CMIP6_hold;

        % CMIP6 CO2
        CO2_CMIP6 = cumsum_co2_emissions_int{n}(101:351);
        CO2_CMIP6 = repmat(CO2_CMIP6, [12, 1]);

        

            
        % _____ Extract linear period _____ %
        x = CO2_CMIP6;
        y = SIA_CMIP6;
        x(y<1) = nan;
        y(y<1) = nan; 

        pi_SIA = mean(SIA_CMIP6(:,1:51), 2);
        constant_co2 = max(max(CO2_CMIP6));
        
        x(y > pi_SIA-(pi_SIA.*0.1)) = nan;
        y(y > pi_SIA-(pi_SIA.*0.1)) = nan;
        
        % Only plot 5th-95th percentiles
        y(x == constant_co2) = nan;
        x(x == constant_co2) = nan;
        
        x_linear = x;
        y_linear = y;




        % ________ Calculate the Sensitivity over linear period ________ % 
        clear coefficients
        for j = 1:12
            y_hold = y_linear(j,:);
            x_hold = x_linear(j,:);
            
            x_hold(isnan(x_hold) == 1) = [];
            y_hold(isnan(y_hold) == 1) = [];
            
            SIA_emul = y_hold;
            co2_emul = x_hold;
            
            if length(SIA_emul) < 2
                SIA_emul = nan;
                co2_emul = nan;
            end
    
            coefficients = polyfit(co2_emul, SIA_emul, 1);
            if coefficients == 0
                coefficients = [nan, nan];
            end
            magicc_sens_hold_all_cmip6{j,n} = coefficients(1);            
        end  

        
        
                
        % ______Find cumulative emission when SIA falls below 1million km2 & if it doesn't then assign a NaN_____ % 
        for j = 1:12
            ind_1_SIA = find(SIA_CMIP6(j,:) <= 1);
            if isempty(ind_1_SIA) ~= 1
                ind_1_co2 = CO2_CMIP6(j,ind_1_SIA(1));
            elseif isnan(y_linear(j,1)) == 1
                ind_1_co2 = nan;
            else
                x_hold = x_linear(j,:);
                y_hold = y_linear(j,:);
    
                x_hold(isnan(x_hold) == 1) = [];
                y_hold(isnan(y_hold) == 1) = [];
    
                ind_1_co2 = interp1(y_hold, x_hold, 1,'linear','extrap');   % Extrapolate if SIA doesn't reach 1
            end
    
            cb_1_ensembles_cmip6{j,i} = (ind_1_co2); 
            
        end
    end
    cb_1_ensembles_cmip6_months = cat(1, cb_1_ensembles_cmip6_months, cb_1_ensembles_cmip6); 
end
cb_1_ensembles_cmip6 = cb_1_ensembles_cmip6_months;





% Calculate the remaining carbon budget based on values above from 1859 and 2023
% Carbon Budget since 1750
cb_1_ensembles2 = cell2mat(cb_1_ensembles(:));
cb_1_ensembles = cell2mat(cb_1_ensembles);
years = 1850:2023;
ind = find(years == 2023);
co2_2023 = co2_store(:,ind)';
CB_store = cb_1_ensembles - co2_2023;       % Carbon Budget 2023 EMULATOR
co2_2023_cmip6 = CO2_CMIP6(:,ind)';
CB_store_cmip6_2023 = [];
for i = 1:12
    CB_store_cmip6_2023{i} = cell2mat(cb_1_ensembles_cmip6(:,i)) - (co2_2023_cmip6(1));    % Carbon Budget 2023 CMIP6  
end







%%  Finding winter non-linearity (the CO2 emission the sensitivity increases) using 'b' calibration parameter

close
clc
coefficients = [];
sens_march_acc = [];
b_final_save = [];
sens_before_TP_updated = [];
CB_at_b = [];
year_of_tp = [];
GT_at_b = [];
month_ind = 3;
for n = 1
    for i = 1:6
        % Extract emulated MAGICC SIA     
        SIA = cell2mat(sia_mon_rearrange_final_OC{n}(:,i));
        amt_AMST = cell2mat(sia_mon_rearrange_final_OC{n}(:,i));
        
        % Extract MAGICC cumulative CO2 from 1850-2300
        cum_co2_one = cumsum_co2_emissions_int{n}(101:551);
    
        % Global-mean temperature anomaly
        GT_anom = repmat(GT_ensembles{n}, [12,1]);
        
    
    
        intervalsa = reshape(1:7200,[600,12]);
        for ii = 1:7200
            [row, col] = find(intervalsa == ii);
            b = new_vals_random(col,4);
            b = b + norm_f_diff{row,col}(i);
            b = b * 2;
            b_final_save{ii,n} = b;
            
            ind_icefree = find(SIA(ii,:) <= 1 );
            [minValue,closestIndex] = min(abs(b-amt_AMST(ii,:)));
            constant_co2 = find(cum_co2_one == max(cum_co2_one));
            constant_co2 = constant_co2(1);
    
            if closestIndex ~= 451
                CB_at_b{ii,i} = cum_co2_one(closestIndex);
    
                GT_at_b{ii,i} = GT_anom(ii, closestIndex);
    
                years_tp = 1850:2300;
                year_of_tp{ii,i} = years_tp(closestIndex);
            else
                CB_at_b{ii,i} = NaN;
    
                GT_at_b{ii,i} = NaN;
    
                year_of_tp{ii,i} = NaN;
            end
    
            sens_before_TP_updated_hold = polyfit( cum_co2_one(1:closestIndex), SIA(ii, 1:closestIndex), 1 );
            sens_before_TP_updated{ii,i} = ( sens_before_TP_updated_hold(1) * 1e12 ) / 1e9;
    
            
            if isempty(ind_icefree) ~= 1
                coefficients{ii, i} = cum_co2_one(ind_icefree(1));
    
                sens_march_acc_hold = polyfit(cum_co2_one(closestIndex:ind_icefree(1)), SIA(ii,closestIndex:ind_icefree(1)), 1);
                sens_march_acc{ii, i} = (sens_march_acc_hold(1) * 1e12) / 1e9;
    
    
            elseif ( amt_AMST(ii,451) < b ) || ( constant_co2 <=  closestIndex ) || ( length(closestIndex:constant_co2(1)) < 51)
                coefficients{ii, i} = nan;
                sens_march_acc{ii, i} = nan; 
    
    
            elseif (amt_AMST(ii,451) > b) && ( constant_co2 >  closestIndex )
                
                polyfita = polyfit(cum_co2_one(closestIndex:constant_co2), SIA(ii,closestIndex:constant_co2) , 1);
                xFit = linspace(min(cum_co2_one(closestIndex:constant_co2)), max(cum_co2_one(closestIndex:constant_co2)), 1000);
                yFit = polyval(polyfita , xFit);
                coefficients{ii, i} = interp1( yFit , xFit, 1, 'linear', 'extrap');
    
                sens_march_acc_hold = polyfit(cum_co2_one(closestIndex:constant_co2), SIA(ii,closestIndex:constant_co2), 1);
                sens_march_acc{ii, i} = (sens_march_acc_hold(1) * 1e12) / 1e9;

            end
        end
    end
end
coefficients = cell2mat(coefficients);
sens_march_acc = cell2mat(sens_march_acc);
CB_MARCH_EMULATOR = coefficients;
b_final_save = cell2mat(b_final_save);
CB_at_b = cell2mat(CB_at_b);
CB_at_b_2023 = CB_at_b;

for i = 1:6
    cum_co2_one = cumsum_co2_emissions_int{1}(101);
    CB_at_b(:,i) = CB_at_b(:,i) + cum_co2_one;

    cum_co2_one = cumsum_co2_emissions_int{1}(274);
    CB_MARCH_EMULATOR_2023(:,i) = CB_MARCH_EMULATOR(:,i) - cum_co2_one;
    CB_at_b_2023(:,i) = CB_at_b_2023(:,i) -  cum_co2_one;
end
cum_co2_one = cumsum_co2_emissions_int{1}(101:551);
med = median( CB_at_b_2023(:,3), 'omitnan' );
[~, ind_of_med_CB_TP] = min(abs( med - cum_co2_one ));





%% ____ Fig 8 ____ %% 

clc
close all
figure(3)
set(gcf, 'Units', 'Inches', 'Position', [.4 .4 19 10])
colorss1 = [1 0.8 0.8; 0.8 0.8 1; 0 0.6906 0.5];


yticks_hold = [];
yticks_hold_2 = [];


h = tiledlayout(1,1);
ax1 = axes(h);
ind_month = [7:12, 1:6];
for i = ind_month

     xlabel(' Cumulative CO_2 emissions since 1850 (GtCO_2)  ')
   
    ax1.XAxisLocation = 'top';

end

ax2 = axes(h);
for i = ind_month



    % ____ Emulator: CB 2023 ____ %
    CB23 = boxplot(ax2, CB_2024_PUB{i}, 'positions', i+0.15, 'color', colorss1(1,:), 'plotstyle', 'compact', 'orientation', 'horizontal'); hold on
    set(CB23(length(CB23),:),'Visible','off')         % Remove the outliers

    hAx = gca;
    lines = hAx.Children;
    lines = lines(1);
    bp_cmip6 = findobj(lines, 'tag', 'Box');
    bp_cmip6.XData(1,1) = quantile(CB_2024_PUB{i}, 0.25);
    bp_cmip6.XData(1,2) = quantile(CB_2024_PUB{i}, 0.75);


    if ismember(i, [7:12])
        % ____ CMIP6: CB 2023 ____ %
        cmip6_CB = boxplot(ax2, CB_store_cmip6_2023{i}, 'positions', i-0.15, 'color', [1 0 0], 'plotstyle', 'compact', 'orientation', 'horizontal'); 
        set(cmip6_CB(length(cmip6_CB),:),'Visible','off')         % Remove the outliers
    
        hAx = gca;
        lines = hAx.Children;
        lines = lines(1);
        bp_cmip6 = findobj(lines, 'tag', 'Box');
        bp_cmip6.XData(1,1) = quantile(CB_store_cmip6_2023{i}, 0.25);
        bp_cmip6.XData(1,2) = quantile(CB_store_cmip6_2023{i}, 0.75);
    
    end
    
    
    a = get(get(gca,'children'),'children');   % Get the handles of all the objects
    for ii = 1:length(a)
        anew = a{ii};
        t = get(anew,'tag');   % List the names of all the objects 
        idx=strcmpi(t,'box');  % Find Box objects
        boxes=anew(idx);          % Get the children you need
        set(boxes,'linewidth',10); % Set width
    end
    
    
    
    hold on
    
    yline((i-0.15) - 0.2750, 'color', [0 0 0]+0.7, 'linewidth', 1.5)


    set(ax2, 'Position', get(ax1, 'Position'), 'Color', 'none');
    ax2.XAxisLocation = 'bottom';   
    ax2.YAxisLocation = 'right';




    ax1.XLim = [1990, 16410];
    ax2.XLim = [-500, 14000];


    xlabpos = xlabel(' Remaining carbon budget (cumulative CO_2 emissions from 2024) (GtCO_2)  ');
    
    yticks_hold = cat(2, yticks_hold, [i-0.15, i+0.15]);
    yticks_hold_2 = cat(2, yticks_hold_2, [i+0.0750]);

    grid
end


set(gca, 'YDir', 'reverse');

% Carbon Budget Shaded Area below 0 Gt CO2
ind_ipcc = -800;
x = ind_ipcc:0;
x2 = [x, fliplr(x)]; 
curve1 = repelem(15, length(x));
curve2 = repelem(-2, length(x));
inBetween = [curve1, fliplr(curve2)];
c = [0 0 0]+0.7;
f2 = fill(x2, inBetween, c, 'edgecolor','none'); 
set(f2,'facealpha',.3)
xline(0, 'color', [0 0 0], 'linewidth', 1.5)


set(ax1, 'fontsize', 14)
set(ax2, 'fontsize', 14)

ax2.YLim = [0.55 12.5];
ax1.YLim = [0.55 12.5]; 


set(ax1, 'ytick', [])


h.Position = [0.23 0.3 0.45 0.6];


set(ax2, 'ytick', [yticks_hold_2(7:12), yticks_hold_2(1:6)], 'YTickLabels', compose(string(month_label)), 'fontsize', 14)

grid
clc


% xline(3000)



% ___ IPCC CB ___ % 
axes('Position',[.2414 .16 .1128 .1])

% IPCC 1.5 degrees warming
hold on
correction_2020_2023 = 38*3;                % IPCC (from 1st Jan 2020) % they have been corrected for 2023 by assuming 38GtCO2/yr is lost (38*3)
IPCC_CB = boxplot(([80, 496]), 'positions', 1.0, 'color', [0 0 0]+0.5, 'plotstyle', 'compact', 'orientation', 'horizontal');

a = get(get(gca,'children'),'children');   % Get the handles of all the objects
t = get(a,'tag');   % List the names of all the objects 
idx=strcmpi(t,'box');  % Find Box objects
boxes=a(idx);          % Get the children you need
set(boxes,'linewidth',10); % Set width

Median_dot_lines = findobj(IPCC_CB, 'type', 'line', 'Tag', 'MedianOuter');     % Make median bigger and bolder
Median_dot_lines.XData = 275;

Median_dot_lines = findobj(IPCC_CB, 'type', 'line', 'Tag', 'MedianInner');
Median_dot_lines.XData = 275;


% IPCC 2 degrees warming
IPCC_CB_2 = boxplot(([900, 2300]-correction_2020_2023), 'positions', 1.2, 'color', [0 0 0]+0.5, 'plotstyle', 'compact', 'orientation', 'horizontal');

a = get(get(gca,'children'),'children');   % Get the handles of all the objects
a = a{1};
t = get(a,'tag');   % List the names of all the objects 
idx=strcmpi(t,'box');  % Find Box objects
boxes=a(idx);          % Get the children you need
set(boxes,'linewidth',10); % Set width

Median_dot_lines = findobj(IPCC_CB_2, 'type', 'line', 'Tag', 'MedianInner');
Median_dot_lines.XData = 1150; 
Median_dot_lines = findobj(IPCC_CB_2, 'type', 'line', 'Tag', 'MedianOuter');
Median_dot_lines.XData = 1150;


ind = ['2', char(176), 'C'];
ind_2 = ['1.5', char(176), 'C'];

set(gca, 'ytick', [1, 1.2],'YTickLabels', compose([string(ind_2), string(ind)]), 'fontsize', 14)


% Carbon Budget Shaded Area below 0 Gt CO2
x = ind_ipcc:0;
x2 = [x, fliplr(x)]; 
curve1 = repelem(15, length(x));
curve2 = repelem(-2, length(x));
inBetween = [curve1, fliplr(curve2)];
c = [0 0 0]+0.7;
f2 = fill(x2, inBetween, c, 'edgecolor','none'); 
set(f2,'facealpha',.3)
xline(0, 'color', [0 0 0], 'linewidth', 1.5)


clear ax
ax = gca;
ax.XLim = [-500, 3000];
ax.YLim = [0.9 1.3];
set(ax, 'fontsize', 14)
set(ax, 'YAxisLocation', 'right');
grid

xline(median(CB_2024_PUB{9}, 'omitnan'), '--')


y = [repelem(0.5975, 2)];
x = [0.16 0.75];
annotation('line', x, y, 'LineStyle', '--', 'Linewidth', 2)

text( -1000, 1.86, 'Linear Mode', 'color', [0 0 0]+0.5, 'Rotation', 90, 'FontSize', 14 )
text( -1000, 3.15, 'Non-linear Mode',  'color', [0 0 0]+0.5, 'Rotation', 90, 'FontSize', 14 )

xlabpos.Position(2) = xlabpos.Position(2) + 3;


% Add legend
legend_labels = ["This Study", "CMIP6", "Paris Agreement Warming Targets"];
tt2 = legend( [ CB23(2), cmip6_CB(2), IPCC_CB(2)], legend_labels, 'location', 'southeast', 'fontsize', 15);
tt2.Box = 'off';
tt2.Position(2) = tt2.Position(2) - 0.2;


% Save
set(gcf, 'PaperOrientation', 'landscape')
set(gcf,'PaperSize',[53 29.7000]);

% % Save Normalised
% temp=['CB_pub_plot', '.pdf']; 
% saveas(gca,temp); 





