%% Sensitivity from OC and OCCAA Emulators


%% ______ Fig 2.16 (Chapter 2, Section 2.5.4) (only March and September) _______ % 


% Initialise Plot  
movmean_ind = 20;
years = 1850:2300;
index_1974_2014 = [130:165]; % to extract 1979-2014

% Get colour RGB triplets
colorss1 = [0.8 0.8 1; 1 0.8 0.8; 0 0.6906 0.5];
colorss2 = [0 0 1; 1 0 0; 0 0.3906 0];


close
figure(2)
set(gcf, 'Units', 'Inches', 'Position', [.4 .4 19 10])
hold on
h = subplot(1,1,1);
for month_ind = [3,9]

    
    % ________ Calculate the observational Sensitivity _________ % 
    obs_sens = [];
    obs_sens_hold = [];
    clear coefficients
    for nn = [1,3:6]    
        for m = [3:6]
            GT = global_mean_redshaped{nn};
            GT = movmean(GT, movmean_ind);        % Apply a running mean of 20 years
            GT = GT(index_1974_2014);             % Extract 1979-2014 period

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

            coefficients = polyfit(GT, SIA_obs, 1);         % Extract linear regression coefficients that are the sensitivity
            obs_sens_hold{nn,m} = coefficients(1);
        end
    end
    obs_sens = cell2mat(obs_sens_hold);




    % Calculate the Sensitivity: OCCAA Emulator
    clear coefficients
    magicc_sens_hold = [];
    OCCAA_sens = [];
    for n = 1:1
        for j = 1:12
            for ii = 1:600
                SIA_emul = sia_mon_rearrange_OCCAA{n}(:,month_ind);
                SIA_emul = movmean(SIA_emul{j}(ii,:), movmean_ind);
                SIA_emul = SIA_emul(index_1974_2014);

                GT = GT_ensembles{n};
                GT_emul = movmean(GT(ii,:), movmean_ind);
                GT_emul = GT_emul(index_1974_2014);

                coefficients = polyfit(GT_emul, SIA_emul, 1);
                magicc_sens_hold{ii,j} = coefficients(1);
            end  
        end
        OCCAA_sens{n} = cell2mat(magicc_sens_hold);   
    end



    % Calculate the Sensitivity: OC Emulator 
    clear coefficients
    magicc_sens_hold = [];
    OC_sens = [];
    magicc_sens_hold_1 = [];
    for n = 1:1
        SIA = sia_mon_rearrange_final_OC{n}(:, month_ind);
        for j = 1:12
            for ii = 1:600
                SIA_emul = movmean(SIA{j}(ii,:), movmean_ind);
                SIA_emul = SIA_emul(index_1974_2014);

                GT = GMST_random_arrangment{n};
                GT_emul = movmean(GT(ii,:), movmean_ind); 
                GT_emul = GT_emul(index_1974_2014);

                coefficients = polyfit(GT_emul, SIA_emul, 1);
                magicc_sens_hold_1{ii} = coefficients(1);
            end  
            magicc_sens_hold{j} = (cell2mat(magicc_sens_hold_1))';
        end
        OC_sens{n} = cell2mat(magicc_sens_hold(:));   
    end





    % Calculate the CMIP6 Sensitivity
    clear coefficients
    cmip6_sens = [];
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
    cmip6_sens = cell2mat(cmip6_sens_hold);



    % ______ MAKE BOXPLOT _______ %     
    if month_ind == 3
        ind = 0.25;
    elseif month_ind == 9
        ind = 0.55;
    end
    
    % Plot observational sensitivity
    obs_sens_bp = boxplot(obs_sens(:), 'positions', [ind-0.06], 'color', [0 0 0]+0.7, 'plotstyle', 'compact'); hold on
    hAx = gca;
    lines = hAx.Children;
    lines = lines(1);
    bp_cmip6 = findobj(lines, 'tag', 'Box');
    bp_cmip6.YData(1,1) = quantile(obs_sens(:), 1.0);   % Set quartile range to show in boxplot
    bp_cmip6.YData(1,2) = quantile(obs_sens(:), 0.0);

    
    
    % Dirk's Sens Data (SIMIP Community, 2020)
    if month_ind == 3
        boxplot([-1.35, -1.85], 'positions', [ind+0.06], 'color', colorss1(1,:), 'Labels', [' '], 'plotstyle','compact'); hold on % MAGICC using the UPPER BOUNDARY of our observationally constrained AA
    else
        boxplot([-3.81, -4.39], 'positions', [ind+0.06], 'color', colorss1(1,:), 'Labels', [' '], 'plotstyle','compact'); hold on % MAGICC using the UPPER BOUNDARY of our observationally constrained AA
    end


    % Plot CMIP6 sensitivity    
    boxplot(cmip6_sens, 'positions', [ind-0.03], 'color', [1 0 0], 'plotstyle','compact'); hold on  
    hAx = gca;
    lines = hAx.Children;
    lines = lines(1);
    bp_cmip6 = findobj(lines, 'tag', 'Box');
    bp_cmip6.YData(1,1) = quantile( cmip6_sens, 0.25);      % specify to plot interquartile range
    bp_cmip6.YData(1,2) = quantile( cmip6_sens, 0.75);


        
    % Plot sensitivity from OCCAA emulator  
    cmip6_aa_sens = boxplot(OCCAA_sens{1}(:), 'positions', [ind], 'color', [0    0.5000    0.5000], 'Labels', [' '], 'plotstyle','compact'); hold on % MAGICC using the UPPER BOUNDARY of our observationally constrained AA
    set(cmip6_aa_sens(length(cmip6_aa_sens),:),'Visible','off')     % Remove the outliers
    hAx = gca;
    lines = hAx.Children;
    lines = lines(1);
    bp_cmip6 = findobj(lines, 'tag', 'Box');
    bp_cmip6.YData(1,1) = quantile( OCCAA_sens{1}(:), 0.25);      % specify to plot interquartile range
    bp_cmip6.YData(1,2) = quantile( OCCAA_sens{1}(:), 0.75);


    % Plot sensitivity from OC emulator  
    SIE_constrained_bp = boxplot(OC_sens{1}, 'positions', [ind+0.03], 'color', colorss1(2,:), 'Labels', [' '], 'plotstyle','compact'); hold on % MAGICC using the UPPER BOUNDARY of our observationally constrained AA
    set(SIE_constrained_bp(length(SIE_constrained_bp),:),'Visible','off')     % Remove the outliers
    hAx = gca;
    lines = hAx.Children;
    lines = lines(1);
    bp_cmip6 = findobj(lines, 'tag', 'Box');
    bp_cmip6.YData(1,1) = quantile( OC_sens{1}, 0.25);  % specify to plot interquartile range
    bp_cmip6.YData(1,2) = quantile( OC_sens{1}, 0.75);



    

    % _____ Plot Edits ____ % 
    a = get(get(gca,'children'),'children');    % Set width of boxplots
    for ii = 1:length(a)
        t = get(a{ii},'tag'); 
        idx=strcmpi(t,'box'); 
        a_new = a{ii};
        boxes=a_new(idx);         
        set(boxes,'linewidth',20); 
    end

    
    % Add plot labels and edits
    sgtitle([num2str(years(index_1974_2014(1))), '-', num2str(years(index_1974_2014(end))), ' Sensitivity (dSIA/dGTA)'], 'fontsize', 28, 'interpreter', 'none')

    if ismember(month_ind, [1,5,9])
        ylabel('million km^2/ \circC')
    end
    set(gca, 'fontsize', 22)
        
    
    length_curve_2 = 0.7;
    xlim([0 length_curve_2])    
    ylim([-6 0])    
    xticks([0.25,0.55])
    xticklabels(string({month_label{3}, month_label{9}}));  % Add month labels to each cluster of boxplots    

 
    grid
    hold on
    
end

h.Position = [0.3 0.075 0.4 0.85];      % Change size of plot
    
grid 
hLegend = legend(findall(gca,'Tag','Box'), {'OC Emulator', 'OCCAA Emulator',  'CMIP6', 'Niederdrenk and Notz, (2018)', 'Observations'}, 'location', 'southwest', 'fontsize', 18);
hLegend.Box = 'off'; 

% Save
set(gcf, 'PaperOrientation', 'landscape')
set(gcf,'PaperSize',[51 30]);

% % Save Normalised
% temp=['Fig2_16', '.pdf']; 
% saveas(gca,temp); 


























%% ______ Fig A.10 (Chapter 2, Supplementary) (all months) _______ % 


% Initialise Plot  
movmean_ind = 20;
years = 1850:2300;
index_1974_2014 = [130:165]; % to extract 1979-2014

% Get colour RGB triplets
colorss1 = [0.8 0.8 1; 1 0.8 0.8; 0 0.6906 0.5];
colorss2 = [0 0 1; 1 0 0; 0 0.3906 0];




close
figure(2)
set(gcf, 'Units', 'Inches', 'Position', [.4 .4 19 10])
hold on
h = subplot(1,1,1);
for month_ind = [1:12]
    
    
    %  _______ Calculate Sensitivity ______ % 

    % Calculate the observational Sensitivity 
    obs_sens = [];
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
    obs_sens = cell2mat(obs_sens_hold);




    % Calculate the Sensitivity: OCCAA Emulator
    clear coefficients
    magicc_sens_hold = [];
    OCCAA_sens = [];
    for n = 1:1
        for j = 1:12
            for ii = 1:600
                SIA_emul = sia_mon_rearrange_OCCAA{n}(:,month_ind);
                SIA_emul = movmean(SIA_emul{j}(ii,:), movmean_ind);
                SIA_emul = SIA_emul(index_1974_2014);

                GT = GT_ensembles{n};
                GT_emul = movmean(GT(ii,:), movmean_ind);
                GT_emul = GT_emul(index_1974_2014);

                coefficients = polyfit(GT_emul, SIA_emul, 1);
                magicc_sens_hold{ii,j} = coefficients(1);
            end  
        end
        OCCAA_sens{n} = cell2mat(magicc_sens_hold);   
    end




    % Calculate the Sensitivity: OC Emulator 
    clear coefficients
    magicc_sens_hold = [];
    OC_sens = [];
    magicc_sens_hold_1 = [];
    for n = 1:1
        SIA = sia_mon_rearrange_final_OC{n}(:, month_ind);
        for j = 1:12
            for ii = 1:600
                SIA_emul = movmean(SIA{j}(ii,:), movmean_ind);
                SIA_emul = SIA_emul(index_1974_2014);

                GT = GMST_random_arrangment{n};
                GT_emul = movmean(GT(ii,:), movmean_ind); 
                GT_emul = GT_emul(index_1974_2014);

                coefficients = polyfit(GT_emul, SIA_emul, 1);
                magicc_sens_hold_1{ii} = coefficients(1);
            end  
            magicc_sens_hold{j} = (cell2mat(magicc_sens_hold_1))';
        end
        OC_sens{n} = cell2mat(magicc_sens_hold(:));   
    end






    % Calculate the CMIP6 Sensitivity
    clear coefficients
    cmip6_sens = [];
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
    cmip6_sens = cell2mat(cmip6_sens_hold);



    % _____ Plot sensitivity as boxplots ______ %
    
    % Plot Observational Sensitivity
    boxplot(obs_sens(:), 'positions', [month_ind-0.3], 'color', [0 0 0]+0.7, 'plotstyle', 'compact');
    hAx = gca;
    lines = hAx.Children;
    lines = lines(1);
    bp_cmip6 = findobj(lines, 'tag', 'Box');
    bp_cmip6.YData(1,1) = quantile(obs_sens(:), 1.0);
    bp_cmip6.YData(1,2) = quantile(obs_sens(:), 0.0);
    

    % Plot CMIP6 Sensitivity
    boxplot(cmip6_sens, 'positions', [month_ind-0.15], 'color', [1 0 0], 'plotstyle','compact'); hold on  
    

    % Plot OCCAA Sensitivity
    bp_cmip6_AA = boxplot(OCCAA_sens{1}(:), 'positions', [month_ind], 'color', [0    0.5000    0.5000], 'Labels', [' '], 'plotstyle','compact'); hold on
    set(bp_cmip6_AA(length(bp_cmip6_AA),:),'Visible','off')


    % Plot OC Sensitivity
    bp_obs_AA = boxplot(OC_sens{1}(:), 'positions', [month_ind+0.15], 'color', colorss1(2,:), 'Labels', [' '], 'plotstyle','compact'); hold on
    set(bp_obs_AA(length(bp_obs_AA),:),'Visible','off')
  
    


    % _____ Plot Edits ____ % 
    a = get(get(gca,'children'),'children');        % Set width of boxplots
    for ii = 1:length(a)
        t = get(a{ii},'tag'); 
        idx=strcmpi(t,'box'); 
        a_new = a{ii};
        boxes=a_new(idx);         
        set(boxes,'linewidth',10); 
    end

    
    % Add plot labels and edits
    text = [num2str(years(index_1974_2014(1))), '-', num2str(years(index_1974_2014(end))), ' Sensitivity (dSIA/dGTA) '];
    ylabel([text, '(million km^2/ \circC)'])
    set(gca, 'fontsize', 22)
          
    length_curve_2 = 13;
    xlim([0 length_curve_2])    
    ylim([-8 0])    
    xticks([1:13])
	xticklabels(string(month_label));
    xtickangle(25) 
    hold on
    
end
hLegend = legend(findall(gca,'Tag','Box'), {'OC Emulator', 'OCCAA Emulator', 'CMIP6', 'Observations'}, 'location', 'southwest');
legend box off
grid 

% Save
set(gcf, 'PaperOrientation', 'landscape')
set(gcf,'PaperSize',[51 30]);

% % Save Normalised
% temp=['SupFigA.10', '.pdf']; 
% saveas(gca,temp); 
