%% Probability Calculations: using OC, OCCAA, CE emulators and CMIP6 data: Thesis Chapter 4

%% Calculate the 1st year an ice-free Arctic Ocean occurs from our OC emulator

% Initialisation
closestIndex_save_ALL_YEARS2 = [];
closestIndex_save_ALL_YEARS_cmip6_AA2= [];
closestIndex_save_ALL_YEARS = [];
closestIndex_save_ALL_YEARS_cmip6_AA = [];
data_YEARS = [];
data_YEARS_cmip6_AA = [];
closestIndex_save_CMIP6_AA = [];
closestIndex_save = [];
closestIndex_save_GMST = [];
closestIndex_save_GMST_CMIP6_AA = [];
closestIndex_save_GMST_hold = [];
closestIndex_save_GMST_hold_CMIP6_AA = [];
closestIndex_save_GMST_SSPs_CMIP6_AA = [];
closestIndex_save_GMST_SSPs = [];
new_SIA = repmat((1:451), [600, 1]);

tic
for n = 1:3     % SSP emission scenarios
    disp('__________')
    GMST = GT_ensembles{n};     % MAGICC GMST ensemble    
    for i = [1:12]       % Months (3= March, 9= September)
        
        disp('__NEW MODEL__')
        
        
        % Extract our emulated SIA: when we observationally constrain the Arctic Amplification
        SIA = sia_mon_rearrange_final_OC{n}(:,i);    
        
        % Find closest year SIA falls below 1 million km2 (ice-free conditions)
        data = cellfun(@(v) new_SIA + (sign(inf * (v <= 1)) - 1), SIA, 'un', 0);
        data_YEARS{i} = cellfun(@(v) v+1849, data, 'un', 0);
        data2 = cellfun(@(v) isnan(v), data, 'un', 0);        
        [row col] = cellfun(@(v) min(v, [], 2), data2, 'un', 0);       
        closestIndex = cellfun(@(v) v + (sign(inf * (v ~= 1)) - 1), col, 'un', 0);
        
        
        
            
        % Extract our emulated SIA: when the Arctic Amplification is CMIP6 calibrated
        SIA_CMIP6 = sia_mon_rearrange_OCCAA{n}(:,i);
        
        % Find closest year SIA falls below 1 million km2 (ice-free conditions)
        data_CMIP6_AA = cellfun(@(v) new_SIA + (sign(inf * (v <= 1)) - 1), SIA_CMIP6, 'un', 0);
        data_YEARS_cmip6_AA{i} = cellfun(@(v) v+1849, data_CMIP6_AA, 'un', 0);
        data_CMIP6_AA2 = cellfun(@(v) isnan(v), data_CMIP6_AA, 'un', 0);        
        [row col] = cellfun(@(v) min(v, [], 2), data_CMIP6_AA2, 'un', 0);       
        closestIndex_CMIP6_AA = cellfun(@(v) v + (sign(inf * (v ~= 1)) - 1), col, 'un', 0);


        closestIndex_save_ALL_YEARS2{n,i} = closestIndex;
        closestIndex_save_ALL_YEARS_cmip6_AA2{n,i} = closestIndex_CMIP6_AA;

        
        
        % Extract the GMST the first time SIA falls below 1 million km2 in each model
        closestIndex_save_GMST_hold = [];
        for j = 1:12
            idx = bsxfun(@eq, cumsum(ones(size(GMST)), 2), closestIndex_CMIP6_AA{j});   % when using the CMIP6 calibrated Arctic Amplification
            Result = sum(GMST.*idx, 2);
            Result(Result==0) = nan;
            closestIndex_save_GMST_hold_CMIP6_AA{j} = Result;

            
   
            % when using our observationally constrained Arctic Amplification
            idx = bsxfun(@eq, cumsum(ones(size(GMST)), 2), closestIndex{j});
            Result = sum(GMST.*idx, 2);
            Result(Result==0) = nan;
            closestIndex_save_GMST_hold{j} = Result;
        end
        

        % Save GMST: for each month
        closestIndex_save_GMST{i} = cell2mat(closestIndex_save_GMST_hold(:));        
        closestIndex_save_GMST_CMIP6_AA{i} = cell2mat(closestIndex_save_GMST_hold_CMIP6_AA(:));        
    end
    
    % Save YEARS: for each SSP scenario
    closestIndex_save_ALL_YEARS{n} = data_YEARS;
    closestIndex_save_ALL_YEARS_cmip6_AA{n} = data_YEARS_cmip6_AA;

    
    % Save GMST: for each SSP scenario
    closestIndex_save_GMST_SSPs{n} = closestIndex_save_GMST;
    closestIndex_save_GMST_SSPs_CMIP6_AA{n} = closestIndex_save_GMST_CMIP6_AA;
end
toc








%% Calculate the 1st year of an ice-free Arctic Ocean from CMIP6 data


% Initialisation
closestIndex_save_YEAR_CMIP6 = [];
closestIndex_save_GMST_CMIP6 = [];
closestIndex_save_YEAR_SSPs_CMIP6 = [];
closestIndex_save_GMST_SSPs_CMIP6 = [];
closestIndex_save_GMST_hold_CMIP6_AA_REAL_CMIP6 = [];
new_SIA = repmat((1:251), [1, 1]);

% Rearrange the CMIP6 SIA for plotting
tic
sia_cmip6_rearrange = [];
sia_cmip6_rearrange_final = [];
for n = 1:3
    disp('__________')    
    for j = 1:12
        disp('NEW MODEL')        
        for i = 1:12
            sia_cmip6_rearrange{j,i} = updated_hist_sia_annual_curve_all_models{j,n}(:,i)';
        end
    end
    sia_cmip6_rearrange_final{n} = sia_cmip6_rearrange;
end
toc




tic
for n = 1:3     % SSP emission scenarios
    GMST = GT_ensembles{n};     % MAGICC GMST ensemble    
    for i = [1:12]   % Months (3= March, 9= September)
        
        
        
        % Extract CMIP6 SIA and add a moving mean
        SIA = sia_cmip6_rearrange_final{n}(:,i);
        SIA = cellfun(@(v) movmean(v, 20), SIA,  'un', 0);
        
        % Extract CMIP6 GMST and add a moving mean
        GMST = tas_global{n};
        GMST = cellfun(@(v) movmean(v, 20), GMST,  'un', 0);
       
        
        % Find closest year CMIP6 SIA falls below 1 million km2 (ice-free conditions)
        data = cellfun(@(v) new_SIA + (sign(inf * (v <= 1)) - 1), SIA, 'un', 0);
        data2 = cellfun(@(v) isnan(v), data, 'un', 0);        
        [row col] = cellfun(@(v) min(v, [], 2), data2, 'un', 0);       
        closestIndex = cellfun(@(v) v + (sign(inf * (v ~= 1)) - 1), col, 'un', 0);
        
%         data_YEARS = cellfun(@(v) v+1849, closestIndex, 'un', 0);
        data_YEARS = cellfun(@(v) v+1849, data, 'un', 0);


        closestIndex_save_YEAR_CMIP6{n,i} = closestIndex;
        
        
        % Extract the GMST the first time SIA falls below 1 million km2 in each model
        for j = 1:12
            find_nan = isnan(closestIndex{j});
            if find_nan ~= 1
                idx = bsxfun(@eq, cumsum(ones(size(GMST{j})), 2), closestIndex{j});
                Result = sum(GMST{j}.*idx, 2);
                Result(Result==0) = nan;
                closestIndex_save_GMST_hold_CMIP6_AA_REAL_CMIP6{j} = Result;
            else
                closestIndex_save_GMST_hold_CMIP6_AA_REAL_CMIP6{j} = nan;
            end
        end
        

        
      % Save YEAR: for each month
        closestIndex_save_YEAR_CMIP6{i} = cell2mat(data_YEARS);

        
%       Save GMST: for each month
        closestIndex_save_GMST_CMIP6{i} = cell2mat(closestIndex_save_GMST_hold_CMIP6_AA_REAL_CMIP6);
        
    end

    
%   Save YEARS: for each SSP scenario
    closestIndex_save_YEAR_SSPs_CMIP6{n} = closestIndex_save_YEAR_CMIP6;    
    
%   Save GMST: for each SSP scenario
    closestIndex_save_GMST_SSPs_CMIP6{n} = closestIndex_save_GMST_CMIP6;
end
toc







%% Calculate the year in which GMST reaches IPCC temperature thresholds: 1/ 1.5/ 2 degrees celsius


% Initialisation
years = 1850:2300;
tas_anoms_final = [];
tas_anoms_final_ssps = [];
threshold_tas = [1, 1.5, 2];        % IPCC GMST thresholds
for n = 1:3     % SSP emission scenarios
    for mm = 1:length(threshold_tas)        % each temperature threshold
        for ii = 1:600                      % each MAGICC GMST ensemble 
            tas = GT_ensembles{n}(ii,:);
            ind = find(tas >= threshold_tas(mm));       % find the index nearest to each GMST threshold
            if isempty(ind) ~= 1
                ind = ind(1);
                year_threshold_tas = years(ind);
            else
                year_threshold_tas = nan;           % if the threshold isn't reaches assign a nan
            end
            tas_anoms_final(ii, mm) = year_threshold_tas;
        end
    end
    tas_anoms_final_ssps{n} = tas_anoms_final;      % save the index for each SSP scenario
end



% Find the ensemble median and percentile range of the first year each GMST is reached
year_threshold_tas_median = [];
year_threshold_tas_percentile = [];
for n = 1:3      % SSP emission scenarios
    year_threshold_tas_median = cat(1, year_threshold_tas_median, mean(tas_anoms_final_ssps{n}, 'omitnan'));
    year_threshold_tas_percentile_hold = prctile(tas_anoms_final_ssps{n}, [5 95]);
    year_threshold_tas_percentile_hold = floor(year_threshold_tas_percentile_hold);
    year_threshold_tas_percentile{n} = year_threshold_tas_percentile_hold;
end

% Check if our code suggest the GMST ocurs in 2300 (meaning the the temp never reaches the threshood, so assign it a nan value)
year_threshold_tas_median = floor(year_threshold_tas_median);
[row col] = find(year_threshold_tas_median == 2300);
year_threshold_tas_median(row, col) = nan;








%% _____ The likelihood of an ice-free ocean in each year and at each GMST _______________________________________________________________________________________________________________________________________________________________________%%
% Chapter 4, Section 4.3.1, Fig 4.1
% Only September and March on same graph


% Initialisation
clear text
colorss1 = [1 0.8 0.8; 0.8 0.8 1; 0 0.6906 0.5];
colorss2 = [1 0 0; 0 0 1; 0 0.3906 0];
years = 1850:2300;
P_AA = [];
P_AA_CMIP6 = [];
P_CMIP6 = [];
counter = 0;

close
fig2 = figure(67);
set(gcf, 'Units', 'Inches', 'Position', [.3 .3 19 10]);
for n = 1:3        % each SSP scenario
    counter = counter+1;
    tas = GT_ensembles{n};      % MAGICC GMST
    median_tas = median(tas);   % median of the MAGICC GMST ensemble 
    prctile_tas = prctile(tas, [5 95]); % percentiles of the MAGICC GMST ensemble 
         
    counter = 0;
    for i = [3,9]       % months to plot
        counter = counter + 1;
        h = subplot(2,2,counter);


        % Probability of an ice-free ocean from OC emulator
        data_hist = cell2mat(closestIndex_save_ALL_YEARS{n}{i});
        P = sum(~isnan(data_hist) & data_hist ~= 0)/size(data_hist,1).*100;     % calculate the probability (method explained in section 3.2)
        P_AA{n,i} = P;

        % Extract the median and percentiles of the years each GMST threshod fall at (to add as to the probability plot)
        mediann = year_threshold_tas_median(n,:);
        percentiless_lower = year_threshold_tas_percentile{n}(1,:);
        percentiless_upper = year_threshold_tas_percentile{n}(2,:);

        % Plot the probability
        plota = plot(years, P, 'color', colorss2(n,:), 'LineWidth' , 3); hold on;
        
    



        % CMIP6 AA: Probability of an ice-free ocean from OCCAA emulator
        xprob_CMIP6_AA = closestIndex_save_ALL_YEARS_cmip6_AA{n}{i};
        data_hist = cell2mat(xprob_CMIP6_AA);
        P = sum(~isnan(data_hist) & data_hist ~= 0)/size(data_hist,1).*100;
        P_AA_CMIP6{n,i} = P;

        % Plot the probability
        plota_cmip6_aa = plot(years, P, '-.', 'color', colorss2(n,:), 'LineWidth' , 1.5); hold on; 
        grid




        % CMIP6: Probability of an ice-free ocean when SIA is calculated from the calibrated CMIP6 AA
        if i == 9
            xprob_CMIP6 = closestIndex_save_YEAR_SSPs_CMIP6{n}{i};
            data_hist = xprob_CMIP6;
            P = sum(~isnan(data_hist) & data_hist ~= 0)/size(data_hist,1).*100; 

            P_CMIP6{n,i} = P;
        end





        % ____ Plot Edits ____ %
        % Add labels to plot
        xlabel(' Year ');
        if n == 1
            y_lab_1 = ylabel({'Probability of an'; 'Ice-free Arctic (%)'});
        end
        set(gca, 'fontsize', 20)

        pos_title = title([ month_label{i} ], 'fontsize', 24);  % Title of each column (Month)
    
    


        % Add more labels to plot
        ylim([0 100])
        xlim([1850 2300])
    
    
        % Add lines and labels to indicate the IPCC 'likely', 'very likely', 'unlikely' and 'very unlikely' boundaries
        c = [0 0 0]+0.6;
        yline(10, 'color', c, 'linewidth', 1)
        hold on
        yline(90, 'color', c, 'linewidth', 1)
        yline(33, 'color', c, 'linewidth', 1)
        yline(66, 'color', c, 'linewidth', 1)
    
        text(1855, 5, 'Very Unlikely', 'color', c, 'FontSize',15)
        text(1855, 95, 'Very Likely', 'color', c, 'FontSize',15)
        text(1855, 78, 'Likely', 'color', c, 'FontSize',15)
        text(1855, 23, 'Unlikely', 'color', c, 'FontSize',15)

        if i == 3
            text(1855, 105, 'a)', 'color', [0 0 0], 'FontSize',19, 'Fontweight', 'bold')
        elseif i == 9

            text(1855, 108, 'b)', 'color', [0 0 0], 'FontSize',19, 'Fontweight', 'bold')
        end

    
        % Change the position and size of each panel based on the month
        if n == 3 && i == 3
            posnew = get(h, 'Position');
            posnew(1) = posnew(1) - 0.055;
            posnew(2) = posnew(2) - 0.02;
            posnew(3) = posnew(3) - 0.02;
    	    set(h, 'Position', posnew);
        elseif n == 3 && i == 9
            posnew = get(h, 'Position');
            posnew(1) = posnew(1) - 0.095;
            posnew(2) = posnew(2) - 0.02;
            posnew(3) = posnew(3) - 0.02;
    	    set(h, 'Position', posnew);
        end
      
    end
end


% _________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
% Add GMST probability to the plot

OC_emulator_gmst_likely = [];
OCCAA_emulator_gmst_likely = [];
clear hh
clear tt2
for n = 1:3           % SSP scenario
    counter = 2;
    for i = [3,9]  %  Month to plot     
        counter = counter + 1;
        h = subplot(2,2,counter);   
        hold on


        % Probability of an ice-free ocean from the OC emulator
        data_hist_UR = closestIndex_save_GMST_SSPs{n}{i};
        find_nan = find(isnan(data_hist_UR)==0);
        if isempty(find_nan) ~= 1
            cdf_hold = cdfplot(data_hist_UR); 
            cdf_hold.Color = [1 1 1];
            ind = find(isnan(data_hist_UR)==1); 
            data_hist_hold = data_hist_UR;
            data_hist_hold(ind) = [];
            nbins2 = length(data_hist_hold);
            updated_ydata = cdf_hold.YData .* (nbins2/7200);
            updated_xdata = cdf_hold.XData;
            hh(1) = plot(updated_xdata, updated_ydata, 'color', colorss2(n,:), 'LineWidth', 3);     % Plot the cdf
            uistack(hh(1),'top');  

            delete(cdf_hold)
        end
        OC_emulator_gmst_likely_hold = find(updated_ydata <= 0.66);
        if OC_emulator_gmst_likely_hold(end) ~= length(updated_ydata)
            OC_emulator_gmst_likely{n,i} = updated_xdata(OC_emulator_gmst_likely_hold(end));
        else
            OC_emulator_gmst_likely{n,i} = NaN;
        end
        if n == 3 & i == 9
            clear max
            [OC_emulator_gmst_likely_hold, ind] = max(updated_ydata);
            OC_emulator_gmst_likely{n,i} = updated_xdata(ind);
        end
        if n == 3 & i == 3
            OC_emulator_gmst_likely{n,i} = NaN;
        end


        xline(1.5)
        xline(2)

        % Probability of an ice-free ocean from OCCAA emulator
        data_hist = closestIndex_save_GMST_SSPs_CMIP6_AA{n}{i};            
        find_nan = find(isnan(data_hist)==0);
        if isempty(find_nan) ~= 1
            cdf_hold = cdfplot(data_hist); 
            cdf_hold.Color = [1 1 1];
            ind = find(isnan(data_hist)==1);
            data_hist_hold = data_hist;
            data_hist_hold(ind) = [];
            nbins2 = length(data_hist_hold); 
            updated_ydata = cdf_hold.YData .* (nbins2/7200);
            updated_xdata = cdf_hold.XData;
            hh(2) = plot(updated_xdata, updated_ydata, '-.', 'color', colorss2(n,:), 'LineWidth', 1.5);
            uistack(hh(2),'top');
            delete(cdf_hold)
        end    
        OCCAA_emulator_gmst_likely_hold = find(updated_ydata >= 0.66);
        if isempty(OCCAA_emulator_gmst_likely_hold) ~= 1
            OCCAA_emulator_gmst_likely{n,i} = updated_xdata(OCCAA_emulator_gmst_likely_hold(1));
        elseif isempty(OCCAA_emulator_gmst_likely_hold) == 1
            OCCAA_emulator_gmst_likely{n,i} = NaN;
        end
        if n == 3 & i == 3
            disp('n=1 i=3')
            OCCAA_emulator_gmst_likely{n,i} = NaN;
        end
        if n == 3 & i == 9
            [OCCAA_emulator_gmst_likely_hold, ind] = max(updated_ydata);
            OCCAA_emulator_gmst_likely{n,i} = updated_xdata(ind);
        end



        % Probability of an ice-free ocean from the CMIP6 data
        data_hist = closestIndex_save_GMST_SSPs_CMIP6{n}{i};            
        find_nan = find(isnan(data_hist)==0);
        if isempty(find_nan) ~= 1
            cdf_hold = cdfplot(data_hist); 
            cdf_hold.Color = [1 1 1];
            ind = find(isnan(data_hist)==1);
            data_hist_hold = data_hist;
            data_hist_hold(ind) = [];
            nbins2 = length(data_hist_hold);
            updated_ydata = cdf_hold.YData .* (nbins2/12);
            updated_xdata = cdf_hold.XData;
            hh(2) = plot(updated_xdata, updated_ydata, '-.', 'color', colorss2(n,:), 'LineWidth', 2);
            uistack(hh(2),'top');
            delete(cdf_hold)
        end


        % ____ Plot Edits ____ %
        title(' ')      % Don't add a title as we add one to the top panel      
        if n == 3
            y_lab = ylabel({'CDF of an Ice-free'; 'Arctic Ocean at each GMST'});
        end        
        xlabel('GMST (\circC)')
        set(gca,'FontSize', 20)     
    

        % Add lines and labels to indicate the IPCC 'likely', 'very likely', 'unlikely' and 'very unlikely' boundaries
        c = [0 0 0]+0.6;
        yline(0.10, 'color', c, 'linewidth', 1)
        hold on
        yline(0.90, 'color', c, 'linewidth', 1)
        yline(0.33, 'color', c, 'linewidth', 1)
        yline(0.66, 'color', c, 'linewidth', 1)
    
    
        % Make sure the x and ylabels are correct
        if i == 9
            xticks([0:1:4])
            xlim([0 4])
        elseif i == 3
            xticks([0:1:9])
            xlim([0 9.5])
        end    
        yticks([0:0.2:1]);
        yticklabels({"0", "20", "40", "60", "80", "100"})
        
    
                
        text(0.03, 0.05, 'Very Unlikely', 'color', c, 'FontSize',15)
        text(0.03, 0.95, 'Very Likely', 'color', c, 'FontSize',15)
        text(0.03, 0.75, 'Likely', 'color', c, 'FontSize',15)
        text(0.03, 0.23, 0.05, 'Unlikely', 'color', c, 'FontSize',15)
        if i == 3
            text(0.03, 1.05, 'c)', 'color', [0 0 0], 'FontSize',19, 'Fontweight', 'bold')
        elseif i == 9

            text(0.03, 1.08, 'd)', 'color', [0 0 0], 'FontSize',19, 'Fontweight', 'bold')
        end
        grid

    
    
        % Add errorbars in
        hold on
        if n == 3 & i == 9
            hh(3) = errorbar(1.9500, 0.66, 0.25, 0.25, 'horizontal', 'Linewidth', 3, 'color', [0 0 0], 'capsize', 12);
            uistack(hh(3),'top');
    
            hh(4) = errorbar(2.15, 0.9, 0.65, 0.65, 'horizontal', 'Linewidth', 3, 'color', [0 0 0], 'capsize', 12, 'linestyle', '--');
            hh(4).Bar.LineStyle = 'dashed';
            uistack(hh(4),'top');
    
            hh(5) = errorbar(1.7, 0.68, 0.4, 0.4, 'horizontal', 'Linewidth', 2, 'color', [0 0 0], 'capsize', 12, 'linestyle', ':');
            hh(5).Bar.LineStyle = 'dotted';
            uistack(hh(5),'top');
    
            hh(6) = scatter(2, 0.6, 100, 'k', 'filled');
    
            hh(7) = scatter(1.5, 0.01, 100, 'd', 'k', 'filled');
    
            hh(8) = scatter(2, 0.42, 100, 'd', 'k', 'filled');
    
            hh(9) = scatter(2, 0.39, 100, 's', 'k', 'filled');
    
            hh(10) = scatter(2, 0.80, 120, 'hexagram', 'k', 'filled');
    
            hh(11) = scatter(2, 0.41, 120, 'r', 'filled');
    
        elseif n == 3 && i == 3
            hh(3) = errorbar(8.8, 0.66, 0.25, 0.25, 'horizontal', 'Linewidth', 3, 'color', [0 0 0], 'capsize', 12);
            uistack(hh(3),'top');
        end
    
    
        % Change the position and size of each panel based on the month
        if n == 3 && i == 3
            posnew = get(h, 'Position');
            posnew(1) = posnew(1) - 0.055;
            posnew(2) = posnew(2) - 0.02;
            posnew(3) = posnew(3) - 0.02;
    	    set(h, 'Position', posnew);
        elseif n == 3 && i == 9
            posnew = get(h, 'Position');
            posnew(1) = posnew(1) - 0.095;
            posnew(2) = posnew(2) - 0.02;
            posnew(3) = posnew(3) - 0.02;
    	    set(h, 'Position', posnew);
        end   


        if n == 3 && i == 9
            % Add legend
            clear allChildren
            allChildren = get(gca, 'Children');
            legend_labels = ["OC Emulator: SSP5-8.5", "OCCAA: SSP5-8.5",...
                            "OC Emulator: SSP2-4.5", "OCCAA: SSP5-8.5",...
                            "OC Emulator: SSP1-2.6", "OCCAA: SSP1-2.6",...
                            "CMIP6 SSP5-8.5", "Sigmond et al, (2018)", "Screen & Williamson, (2017)", ...
                            "Ridley & Blockley, (2018): 2°C", "Ridley & Blockley, (2018): 1.5°C",...
                            "Jahn, (2018)", "Notz and Stroeve, (2018)", "Mahlstein and Knutti, (2012)", "Niederdrenk and Notz, (2018)"];
    
            tt2 = legend(allChildren([48,45, 35,32, 22,19, 1:9]), legend_labels, 'location', 'southeast', 'fontsize', 15);
            uistack(tt2,'top'); 

            
        end
    end
end
legend boxoff
grid

% Posistion legend in the correct place
newPosition = [0.89 0.28 0.01 0.01];
newUnits = 'normalized';
set(tt2,'Position', newPosition,'Units', newUnits);


% set(gcf, 'PaperOrientation', 'landscape')
% set(gcf,'PaperSize',[53 29.7000]); 
% 
% % Save
% temp=['Fig4.1', '.pdf']; 
% saveas(gca,temp); 





%% Likelihood info for tables appendix: sept

clc
% _________________________________________________YEAR___________________________________________________%
% ___________________OC_____________________%
% OC:L
OC_emulator_year_likely = [];
year = 1850:2300;
for n = 1:2
    for i = [3,9]
        OC_emulator_year_likely_hold = find(P_AA{n,i} <= 66);
        if OC_emulator_year_likely_hold(end) ~= 451
            OC_emulator_year_likely{n,i} = year(OC_emulator_year_likely_hold(end));
        elseif OC_emulator_year_likely_hold(end) == 451
            OC_emulator_year_likely{n,i} = NaN;
        end
    end
end


% OC: SSP1-2.6
for n = 3
    for i = [3,9]
        [OC_emulator_year_likely_hold, ind] = max(P_AA{n,i});
        OC_emulator_year_likely{n,i} = [max(P_AA{n,i}), year(ind)];
    end
end




% OC: VL
OC_emulator_year_verylikely = [];
year = 1850:2300;
for n = 1:2
    for i = [3,9]
        OC_emulator_year_likely_hold = find(P_AA{n,i} <= 90);
        if OC_emulator_year_likely_hold(end) ~= 451
            OC_emulator_year_verylikely{n,i} = year(OC_emulator_year_likely_hold(end));
        elseif OC_emulator_year_likely_hold(end) == 451
            OC_emulator_year_verylikely{n,i} = NaN;
        end
    end
end





% ___________________OCCAA_____________________%
% OCCAA: L
OCCAA_emulator_year_likely = [];
year = 1850:2300;
for n = 1:2
    for i = [3,9]
        OC_emulator_year_likely_hold = find(P_AA_CMIP6{n,i} <= 66);
        if OC_emulator_year_likely_hold(end) ~= 451
            OCCAA_emulator_year_likely{n,i} = year(OC_emulator_year_likely_hold(end));
        elseif OC_emulator_year_likely_hold(end) == 451
            OCCAA_emulator_year_likely{n,i} = NaN;
        end
    end
end


% OCCAA: SSP1-2.6
for n = 3
    for i = [3,9]
        [OC_emulator_year_likely_hold, ind] = max(P_AA_CMIP6{n,i});
        OCCAA_emulator_year_likely{n,i} = [max(P_AA_CMIP6{n,i}), year(ind)];
    end
end




% OCCAA: VL
OCCAA_emulator_year_verylikely = [];
year = 1850:2300;
for n = 1:2
    for i = [3,9]
        OC_emulator_year_likely_hold = find(P_AA_CMIP6{n,i} <= 90);
        if OC_emulator_year_likely_hold(end) ~= 451
            OCCAA_emulator_year_verylikely{n,i} = year(OC_emulator_year_likely_hold(end));
        elseif OC_emulator_year_likely_hold(end) == 451
            OCCAA_emulator_year_verylikely{n,i} = NaN;
        end
    end
end


% _________________________________________________TABLES___________________________________________________%


% ___________________OC: Likely _____________________%
clc
% SEPT
% OC: likely: Year & GMST
T1_year_gmst_OC_sept = table( cell2mat(OC_emulator_year_likely(1:2,9)), cell2mat(OC_emulator_gmst_likely(1:2,9)) );
T1_year_gmst_OC_sept.Properties.VariableNames{'Var1'} = 'Year';
T1_year_gmst_OC_sept.Properties.VariableNames{'Var2'} = 'GMST';



% March
% OC: likely: Year & GMST
T1_year_gmst_OC_march = table( cell2mat(OC_emulator_year_likely(1:2,3)), cell2mat(OC_emulator_gmst_likely(1:2,3)) );
T1_year_gmst_OC_march.Properties.VariableNames{'Var1'} = 'Year';
T1_year_gmst_OC_march.Properties.VariableNames{'Var2'} = 'GMST';



SSPS = ['SSP5-8.5'; 'SSP2-4.5'];

T1_OC = table(SSPS, T1_year_gmst_OC_march, T1_year_gmst_OC_sept);
T1_OC.Properties.VariableNames{'T1_year_gmst_OC_march'} = 'OC: March';
T1_OC.Properties.VariableNames{'T1_year_gmst_OC_sept'} = 'OC: September';
T1_OC.Properties.VariableNames{'SSPS'} = 'SSP';
T1_OC



% ___________________OC: Very Likely_____________________%

T1_OC_VL = table(cell2mat(OC_emulator_year_verylikely(1:2,9)) );
T1_OC_VL.Properties.VariableNames{'Var1'} = 'Very Likely September: Year';

T1_OCCAA_VL = table(cell2mat(OCCAA_emulator_year_verylikely(1:2,9)) );
T1_OCCAA_VL.Properties.VariableNames{'Var1'} = 'Very Likely September: Year';

T1_VL = table(SSPS, T1_OC_VL, T1_OCCAA_VL);
T1_VL.Properties.VariableNames{'T1_OC_VL'} = 'OC Emulator';
T1_VL.Properties.VariableNames{'SSPS'} = 'SSPs';
T1_VL.Properties.VariableNames{'T1_OCCAA_VL'} = 'OCCAA Emulator';
T1_VL


% ___________________OCCAA_____________________%
% SEPT
% OC: likely: Year & GMST
T1_year_gmst_OCCAA_sept = table( cell2mat(OCCAA_emulator_year_likely(1:2,9)), cell2mat(OCCAA_emulator_gmst_likely(1:2,9)) );
T1_year_gmst_OCCAA_sept.Properties.VariableNames{'Var1'} = 'Year';
T1_year_gmst_OCCAA_sept.Properties.VariableNames{'Var2'} = 'GMST';



% March
% OC: likely: Year & GMST
T1_year_gmst_OCCAA_march = table( cell2mat(OCCAA_emulator_year_likely(1:2,3)), cell2mat(OCCAA_emulator_gmst_likely(1:2,3)) );
T1_year_gmst_OCCAA_march.Properties.VariableNames{'Var1'} = 'Year';
T1_year_gmst_OCCAA_march.Properties.VariableNames{'Var2'} = 'GMST';



SSPS = ['SSP5-8.5'; 'SSP2-4.5'];

T1_OCCAA = table(SSPS, T1_year_gmst_OCCAA_march, T1_year_gmst_OCCAA_sept);
T1_OCCAA.Properties.VariableNames{'T1_year_gmst_OCCAA_march'} = 'OCCAA: March';
T1_OCCAA.Properties.VariableNames{'T1_year_gmst_OCCAA_sept'} = 'OCCAA: September';
T1_OCCAA.Properties.VariableNames{'SSPS'} = 'SSP';
T1_OCCAA



% ___________________OC: SSP1-2.6_____________________%

% SEPT: SSP1-2.6
% OC: likely: Year & GMST
T1_year_gmst_OC_sept_126 = table( [cell2mat(OC_emulator_year_likely(3,9));cell2mat(OCCAA_emulator_year_likely(3,9))], [cell2mat(OC_emulator_gmst_likely(3,9)); cell2mat(OCCAA_emulator_gmst_likely(3,9))] );
T1_year_gmst_OC_sept_126.Properties.VariableNames{'Var1'} = 'Year';
T1_year_gmst_OC_sept_126.Properties.VariableNames{'Var2'} = 'GMST';


% MARCH: SSP1-26
% OC: likely: Year & GMST
T1_year_gmst_OC_march_126 = table( [cell2mat(OC_emulator_year_likely(3,3)); cell2mat(OCCAA_emulator_year_likely(3,3))], [cell2mat(OC_emulator_gmst_likely(3,3)); cell2mat(OCCAA_emulator_gmst_likely(3,3))]  );
T1_year_gmst_OC_march_126.Properties.VariableNames{'Var1'} = 'Year';
T1_year_gmst_OC_march_126.Properties.VariableNames{'Var2'} = 'GMST';



T1_126 = table(["SSP1-2.6: OC"; "SSP1-2.6: OCCAA"], T1_year_gmst_OC_march_126, T1_year_gmst_OC_sept_126);
T1_126.Properties.VariableNames{'T1_year_gmst_OC_march_126'} = 'OC: March';
T1_126.Properties.VariableNames{'T1_year_gmst_OC_sept_126'} = 'OC: September';
T1_126.Properties.VariableNames{'Var1'} = 'SSP';
T1_126








%% Global Temperature at which a 'likekly' ice-free ocean occurs: Fig 9

close; clc
figure(38)
set(gcf, 'Units', 'Inches', 'Position', [.4 .4 19 10])

h = subplot(1,1,1);
for n = 1:3; hleg(n) = plot( temp_likely_IFC(n,:), 'color', colorss2(n,:), 'linewidth', 3 ); hold on; end
for n = 1; plot( temp_likely_IFC_cmip6(n,:), '--', 'color', colorss2(n,:), 'linewidth', 2 ); hold on; end
ax = gca; ax.FontSize = 18;

% Shading 1.5 degrees of global warming
hold on
x = 0:12;
x2 = [x, fliplr(x)];
curve1 = repelem(2, length(x));
curve2 = repelem(0, length(x));
inBetween = [curve1, fliplr(curve2)];
c = [1 0 0];
f2 = fill(x2, inBetween, c, 'edgecolor','none'); 
set(f2,'facealpha',.1)


ylabel( ' Global-mean Temperature Anomaly ')
xlabel( ' Month ' )
grid

yline(1.5, 'linewidth', 1.5)
yline(2, 'linewidth', 1.5)

set(gca, 'xtick', 1:12);
tick_length = 1;
set(gca,'xticklabel', {[ blanks(tick_length) 'Jan'], [ blanks(tick_length) 'Feb'], [ blanks(tick_length) 'Mar'],...
    [ blanks(tick_length) 'Apr'],[ blanks(tick_length) 'May'],[ blanks(tick_length) 'Jun'],[ blanks(tick_length) 'Jul'],...
    [ blanks(tick_length) 'Aug'], [ blanks(tick_length) 'Sep'], [ blanks(tick_length) 'Oct'],[ blanks(tick_length) 'Nov'], ... 
[ blanks(tick_length) 'Dec'], ''});
set(gca, 'Fontsize', 18)
xlim([1 12])
ylim([0 10])

h.Position(3) = h.Position(3)-0.25;
h.Position(4) = h.Position(4)-0.25;

tt = legend( hleg,  'SSP5-8.5', 'SSP2-4.5', 'SSP1-2.6' );
legend boxoff


set(gcf, 'PaperOrientation', 'landscape')
set(gcf,'PaperSize',[53 29.7000]); 

% % Save
% temp=['Fig4.2', '.pdf']; 
% saveas(gca,temp); 








%% _________ The likelihood of an ice-free ocean in each year: All Months ___________________________________________________________________________________________________________________________________________________________________%%


% Initialisation
clear text
colorss1 = [1 0.8 0.8; 0.8 0.8 1; 0 0.6906 0.5];
colorss2 = [1 0 0; 0 0 1; 0 0.3906 0];
SSP_set = {'SSP585', 'SSP245', 'SSP126'};
years = 1850:2300;
tas_2150 = [];
hnew_cal = [];
P_ssp = [];
P_AA = [];
P_AA_CMIP6 = [];
P_CMIP6 = [];
P_AA_CMIP6_OG = [];
percentsa = [];
percentsa_cmip6_aa = [];
hnew = [];
tas_at_100 = [];
P_ssp_cmip6_aa = [];
plota_cmip6_aa = [];
tas_at_100_cmip6_aa = [];
extract_66_585 = [];
extract_2100 = [];
extract_2300 = [];
counter = 0;
marker_size = 80;
month_ind = 3;      % Month to plot

close; clc
fig2 = figure(67);
set(gcf, 'Units', 'Inches', 'Position', [.3 .3 19 10]);
for n = 1:3        % each SSP scenario
    counter = counter+1;
    tas = GT_ensembles{n};      % MAGICC GMST
    median_tas = median(tas);   % median of the MAGICC GMST ensemble 
    prctile_tas = prctile(tas, [5 95]); % percentiles of the MAGICC GMST ensemble 
         
    counter = 0;
    for i = [1:12]       % months to plot
        counter = counter + 1;
        h = subplot(3,4,counter);


        % Probability of an ice-free ocean when SIA is calculated from our observationally constrained AA
        data_hist = cell2mat(closestIndex_save_ALL_YEARS{n}{i});
        P = sum(~isnan(data_hist) & data_hist ~= 0)/size(data_hist,1).*100;     % calculate the probability (method explained in section 3.2)
        P_AA{n,i} = P;


        years = 1850:2300;
        indnew = P(451); 
        extract_2300(n,i) = max(P);
        extract_2100(n,i) = P(251);


        if n == 1
            years = 1850:2300;
            ind = find(P <= 66); 
            prob_ind = P(ind(end));
            indnew = ind(end);
            indoc = years(indnew);
        end

        % Extract the median and percentiles of the years each GMST threshod fall at (to add as to the probability plot)
        mediann = year_threshold_tas_median(n,:);
        percentiless_lower = year_threshold_tas_percentile{n}(1,:);
        percentiless_upper = year_threshold_tas_percentile{n}(2,:);


        % Plot the probability
        plota = plot(years, P, 'color', colorss2(n,:), 'LineWidth' , 3); hold on;
        
    


        % CMIP6 AA: Probability of an ice-free ocean when SIA is calculated from the calibrated CMIP6 AA
        xprob_CMIP6_AA = closestIndex_save_ALL_YEARS_cmip6_AA{n}{i};
        data_hist = cell2mat(xprob_CMIP6_AA);
        P = sum(~isnan(data_hist) & data_hist ~= 0)/size(data_hist,1).*100;
        P_AA_CMIP6{n,i} = P;

        % Plot the probability
        plota_cmip6_aa = plot(years, P, '-.', 'color', colorss2(n,:), 'LineWidth' , 1.5); hold on; 
        grid

        if n == 1
            indoccaa = find(P <= 66);
            indoccaa= indoccaa(end);
            extract_66_585(i) = years(indoccaa) - indoc;
        end


        % CMIP6: Probability of an ice-free ocean when SIA is calculated from the calibrated CMIP6 AA
        if i == 9
            xprob_CMIP6 = closestIndex_save_YEAR_SSPs_CMIP6{n}{i};
            data_hist = xprob_CMIP6;
            P = sum(~isnan(data_hist) & data_hist ~= 0)/size(data_hist,1).*100;

            P_CMIP6{n,i} = P;
        end


        % Add labels to plot
        if ismember(counter, 9:12)
            xlabel(' Year ');
        end
        if ismember(counter, [1,5,9])
            if n == 1
                y_lab_1 = ylabel({'Probability of an'; 'Ice-free Arctic (%)'});
            end
        end
        set(gca, 'fontsize', 16)

        pos_title = title([ month_label{i} ], 'fontsize', 24);  % Title of each column (Month)
    
    


        % Add more labels to plot
        ylim([0 100])
        xlim([1850 2300])
    
    
        % Add lines and labels to indicate the IPCC 'likely', 'very likely', 'unlikely' and 'very unlikely' boundaries
        c = [0 0 0]+0.6;
        yline(10, 'color', c, 'linewidth', 1)
        hold on
        yline(90, 'color', c, 'linewidth', 1)
        yline(33, 'color', c, 'linewidth', 1)
        yline(66, 'color', c, 'linewidth', 1)
    
        text(1855, 5, 'Very Unlikely', 'color', c, 'FontSize',12)
        text(1855, 95, 'Very Likely', 'color', c, 'FontSize',12)
        text(1855, 78, 'Likely', 'color', c, 'FontSize',12)
        text(1855, 23, 'Unlikely', 'color', c, 'FontSize',12)

      
    end
end



% set(gcf, 'PaperOrientation', 'landscape')
% set(gcf,'PaperSize',[53 29.7000]); 

% % Save
% temp=['FigC.1', '.pdf']; 
% saveas(gca,temp);  





%% _________ The likelihood of an ice-free ocean at each GMST: All Months ___________________________________________________________________________________________________________________________________________________________________%%
% Add GMST probability to the plot

temp_likely_IFC_cmip6aa = [];
temp_likely_IFC = [];
temp_likely_IFC_cmip6 = [];
P_OC_GMST = [];
P_X_OC_GMST = [];

close; clc
fig2 = figure(67);
set(gcf, 'Units', 'Inches', 'Position', [.3 .3 19 10]);
clear hh
clear tt2
for n = 1:3           % SSP scenario
    counter = 2;
    for i = [1:12]  %  Month to plot     
        counter = counter + 1;
        h = subplot(3,4,i);
    
        hold on

        % Probability of an ice-free ocean from the OC emulator
        data_hist_UR = closestIndex_save_GMST_SSPs{n}{i};
        find_nan = find(isnan(data_hist_UR)==0);
        if isempty(find_nan) ~= 1
            cdf_hold = cdfplot(data_hist_UR); 
            cdf_hold.Color = [1 1 1];
            ind = find(isnan(data_hist_UR)==1); 
            data_hist_hold = data_hist_UR;
            data_hist_hold(ind) = [];
            nbins2 = length(data_hist_hold);
            updated_ydata = cdf_hold.YData .* (nbins2/7200);
            updated_xdata = cdf_hold.XData;
            hh(1) = plot(updated_xdata, updated_ydata, 'color', colorss2(n,:), 'LineWidth', 3);     % Plot the cdf
            uistack(hh(1),'top');  
            delete(cdf_hold)
        end
        P_OC_GMST{n,i} = updated_ydata;
        P_X_OC_GMST{n,i} = updated_xdata;

        never_rises_above_threshold = all(updated_ydata <= 0.66);
        if never_rises_above_threshold
            temp_likely_IFC(n,i) = nan;
        else
            temp_likely_IFC_hold = find(updated_ydata <= 0.66);
            temp_likely_IFC(n,i) = updated_xdata( temp_likely_IFC_hold(end) );
        end


        xline(1.5)
        xline(2)

        % Probability of an ice-free ocean from OCCAA emulator
        data_hist = closestIndex_save_GMST_SSPs_CMIP6_AA{n}{i};            
        find_nan = find(isnan(data_hist)==0);
        if isempty(find_nan) ~= 1
            cdf_hold = cdfplot(data_hist); 
            cdf_hold.Color = [1 1 1];
            ind = find(isnan(data_hist)==1);
            data_hist_hold = data_hist;
            data_hist_hold(ind) = [];
            nbins2 = length(data_hist_hold); 
            updated_ydata = cdf_hold.YData .* (nbins2/7200);
            updated_xdata = cdf_hold.XData;
            hh(2) = plot(updated_xdata, updated_ydata, '-.', 'color', colorss2(n,:), 'LineWidth', 1.5);
            uistack(hh(2),'top');
            
            delete(cdf_hold)
        end

        never_rises_above_threshold = all(updated_ydata <= 0.66);
        if never_rises_above_threshold
            temp_likely_IFC_cmip6aa(n,i) = nan;
        else
            temp_likely_IFC_cmip6aa_hold = find(updated_ydata <= 0.66);
            temp_likely_IFC_cmip6aa(n,i) = updated_xdata( temp_likely_IFC_cmip6aa_hold(end) );
        end


        % Probability of an ice-free ocean from the CMIP6 data
        data_hist = closestIndex_save_GMST_SSPs_CMIP6{n}{i};            
        find_nan = find(isnan(data_hist)==0);
        if isempty(find_nan) ~= 1
            cdf_hold = cdfplot(data_hist); 
            cdf_hold.Color = [1 1 1];
            ind = find(isnan(data_hist)==1);
            data_hist_hold = data_hist;
            data_hist_hold(ind) = [];
            nbins2 = length(data_hist_hold);
            updated_ydata = cdf_hold.YData .* (nbins2/12);
            updated_xdata = cdf_hold.XData;
            hh(2) = plot(updated_xdata, updated_ydata, '--', 'color', colorss2(n,:), 'LineWidth', 2);
            hh(2) = plot(updated_xdata, updated_ydata, '--', 'color', colorss2(n,:), 'LineWidth', 0.2);
            uistack(hh(2),'top');
            delete(cdf_hold)
        end

        never_rises_above_threshold = all(updated_ydata <= 0.66);
        if never_rises_above_threshold
            temp_likely_IFC_cmip6(n,i) = nan;
        else
            temp_likely_IFC_cmip6_hold = find(updated_ydata <= 0.66);
            temp_likely_IFC_cmip6(n,i) = updated_xdata( temp_likely_IFC_cmip6_hold(end) );
        end


        % ____ Plot Edits ____ %
        pos_title = title([ month_label{i} ], 'fontsize', 24);  % Title of each column (Month)

        if ismember(i, 5)
            if n == 3
                y_lab = ylabel('CDF of an Ice-free Arctic Ocean at each GMST');
            end     
        else
            ylabel(' ')
        end
        if ismember(i, 9:12)
            xlabel('GMST (\circC)')
        else
            xlabel(' ')
        end
        set(gca,'FontSize', 16)     
    

        % Add lines and labels to indicate the IPCC 'likely', 'very likely', 'unlikely' and 'very unlikely' boundaries
        c = [0 0 0]+0.6;
        yline(0.10, 'color', c, 'linewidth', 1)
        hold on
        yline(0.90, 'color', c, 'linewidth', 1)
        yline(0.33, 'color', c, 'linewidth', 1)
        yline(0.66, 'color', c, 'linewidth', 1)
    
    
        xticks([0:1:9])
        xlim([0 9.5])
        yticks([0:0.2:1]);
        yticklabels({"0", "20", "40", "60", "80", "100"})
        
    
                
        text(0.03, 0.05, 'Very Unlikely', 'color', c, 'FontSize',12)
        text(0.03, 0.95, 'Very Likely', 'color', c, 'FontSize',12)
        text(0.03, 0.75, 'Likely', 'color', c, 'FontSize',12)
        text(0.03, 0.23, 0.05, 'Unlikely', 'color', c, 'FontSize',12)

        if n == 1
            grid
        end
    end
end



% set(gcf, 'PaperOrientation', 'landscape')
% set(gcf,'PaperSize',[53 29.7000]); 

% % Save
% temp=['FigC.2', '.pdf']; 
% saveas(gca,temp); 







%% ____ Heat Map for Publication: Year and GMST ____ %% 


close all; clc
figure(40)
set(gcf, 'Units', 'Inches', 'Position', [.4 .4 19 14])  
SSPS = {'SSP5-8.5'; 'SSP2-4.5'; 'SSP1-2.6'};
for n = 1:3

    h = subplot(4,1,n);
    if n == 1
        h.Position(2) = h.Position(2) + 0.04;
    elseif n == 2
        h.Position(2) = h.Position(2) + 0.06;
    else
        h.Position(2) = h.Position(2) + 0.08;
    end
    h.Position(4) = h.Position(4) + 0.025;
    h.Position(3) = h.Position(3) - 0.5;

    
    
    data = P_AA(n,:)';
    data = cell2mat(data);
    
    
    values = [0 1 10 33 66 90 100]; 
    colors = sky(length(values) - 1);
    colormap(colors);
    
    heatmap = imagesc(data);
    
    c = colorbar('Location', 'southoutside');
    caxis([0 100]);
    
    ticks = linspace(0, 100, length(values));
    set(c, 'Ticks', ticks);
    set(c, 'TickLabels', values);
    
    ylabel(c,'%','FontSize',16)
    if n == 3
        xlabel('Year')
    end
    ylabel('Month')
    % title(SSPS{n})
    
    yticks(1:12);
    ylim([1 12])
    yticks(1:12); 
    yticklabels({'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'}); % Label the months

    ax = gca;
    % ax.XTickLabel = cellfun(@(v) str2num(v)+1850, ax.XTick, 'UniformOutput', false);
    % ax.XTickLabel = cellfun(@(v) str2num(v)+1850, ax.XTickLabel, 'UniformOutput', false);
    % ax.XLim([1850 2300])
    % ax.XLim([0 451])
    xticks(0:50:450);
    ax.XTickLabel = string(ax.XTick+1850);
    xtickangle(20)


    % xline(250, 'linewidth', 2)
    % xline(175, 'linewidth', 1)
    xlim([0 451])

    % xticks(0:10:450);
    % ax.XTickLabel = string(ax.XTick+1850);
    % xlim([151 250])
    
    ax.FontSize = 16;

    grid minor

    if n == 1 || n == 2
        % colorbar('off');
        set(ax, 'XTick', [])
    end

    colorbar('off');
end

% c.Position = [0.1301    0.3668    0.2750    0.02];





% Add GMST
cumulative_matrix_interpolated = [];
% for n = 1

    h = subplot(4,1,4);
    if n == 1
        h.Position(2) = h.Position(2) + 0.02;
    elseif n == 2
        h.Position(2) = h.Position(2) + 0.04;
    else
        h.Position(2) = h.Position(2) + 0.06;
    end
    h.Position(4) = h.Position(4) + 0.025;
    h.Position(3) = h.Position(3) - 0.5;


    gmst_hold = GT_ensembles';
    gmst_hold = cell2mat(gmst_hold);
    max_gmst = max(max(gmst_hold));
    ice_free_temperatures_hold = GMST_prob_setup_final;
    ice_free_temperatures = cellfun(@transpose, ice_free_temperatures_hold, 'UniformOutput', false)';

    
    
    % Define the common temperature range for the x-axis (e.g., 0 to 20 degrees)
    x_common = 0:0.26:20;  
    
    % Number of months
    num_months = length(ice_free_temperatures);
    
    % Preallocate a matrix to hold the interpolated cumulative probabilities
    cumulative_matrix_interpolated = NaN(num_months, length(x_common));
    
    % Calculate and interpolate cumulative probability for each month
    for i = 1:num_months
        % Get the temperatures and retain NaNs
        temp_data = ice_free_temperatures{i};
        if all(isnan(temp_data)) == 1
            continue
        end
        if length(unique(temp_data(~isnan(temp_data)))) < 2
            continue
        end
        
        % Identify non-NaN values and their corresponding indices
        valid_indices = ~isnan(temp_data);
        sorted_temps = temp_data(valid_indices);
        
        % If there are valid temperatures, proceed
        if ~isempty(sorted_temps)
            % Sort temperatures and calculate cumulative probability
            sorted_temps = sort(sorted_temps); 
            [sorted_temps, ia] = unique(sorted_temps, 'stable');
            cumulative_probability = (1:length(sorted_temps)) / length(temp_data); % normalize by total data points including NaNs

            max_gmst = max(sorted_temps);
            
            % Interpolate over the common temperature range
            interpolated_probs = interp1(sorted_temps, cumulative_probability, x_common, 'linear', 'extrap');
            
            % Fill missing values after the last known temperature with the last cumulative probability
            last_prob = cumulative_probability(end);
            interpolated_probs(x_common > max(sorted_temps)) = last_prob;

            % if m == 3
            %     max_temp = max(sorted_temps);
            %     interpolated_probs(x_common > max_temp) = NaN;
            % end
            
            % Assign interpolated probabilities to the matrix
            cumulative_matrix_interpolated(i, :) = interpolated_probs;
        else
            % If all NaNs, assign zero probability across the range
            cumulative_matrix_interpolated(i, :) = 0;
        end
    end
    
    % Plot the heatmap using imagesc
    imagesc(x_common, 1:num_months, cumulative_matrix_interpolated.*100);
    
    values = [0 1 10 33 66 90 100]; 
    colors = sky(length(values) - 1);
    colormap(colors);
    
    c = colorbar('Location', 'southoutside');
    caxis([0 100]);
    
    ticks = linspace(0, 100, length(values));
    set(c, 'Ticks', ticks);
    set(c, 'TickLabels', values)

    ylabel(c,'%','FontSize', 16)
    xlabel('Global Mean Temperature Anomaly (\circC)')
    ylabel('Month')
    
    yticks(1:12);
    ylim([1 12])
    yticks(1:12); 
    yticklabels({'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'}); % Label the months

    % new_xticks = 0:0.5:9;
    ax = gca;
    ax.XMinorTick = 'on';
    % set(ax, 'XTick', new_xticks);


    if m == 1 || m == 2
        if n == 1
            xlim([0 9])
        else
            xlim([0 round(max_gmst)])
        end
    else
        xlim([0 round(max_gmst)])
    end

    ax = gca; ax.FontSize = 16; 

    grid minor

% end

text(0.2, -39, 'SSP5-8.5', 'color', [0 0 0]+0.2, 'FontSize', 16)
text(0.2, -27, 'SSP2-4.5', 'color', [0 0 0]+0.2, 'FontSize', 16)
text(0.2, -15, 'SSP1-2.6', 'color', [0 0 0]+0.2, 'FontSize', 16)

x = [0.085 0.43];
y = [repelem(0.332, 2)];
annotation('line', x, y, 'LineStyle', '--', 'Linewidth', 2)
    

c.Position = [0.1301    0.05    0.2750    0.015];

% Save Normalised
set(gcf, 'PaperOrientation', 'portrait')
set(gcf,'PaperSize',[50 40]);

% % Save Normalised
% temp=['heat_map_year_gmst', '.pdf']; 
% saveas(gca,temp); 

