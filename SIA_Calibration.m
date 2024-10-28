%% SIA Calibration
    

test_input_vals = [0.3505    0.0433    0.3647   -1.4719]; %    Test input values for each parameter to be calibrated

% Opimisation Initiialisation
Save_sia_cal_params = [];
calibrated_SIA_store = [];
hnew3 = [];
tt_store = 1:12;
index_run = [1:12];

close
counter = 0;

for j = index_run

%   Counter for plotting
    counter = counter + 1;


%   CMIP6 SIA
    SIA_store = updated_hist_sia_annual_curve_all_models(j,:);
    SIA_store = cellfun(@(v) movmean(v, 20), SIA_store, 'UniformOutput', false);
    SIA_store = cell2mat(SIA_store(:));


%   Emulated AMST
    x_tas_store = AMST_calibration_bias_corrected(j,:);
    x_tas_store = cellfun(@(v) movmean(v, 20), x_tas_store, 'UniformOutput', false);
    x_tas_store = cell2mat(x_tas_store(:));


%   Rename everything for plotting later
    SIA = SIA_store;
    x_tas = x_tas_store;


%   Use parameterised SIA_max (SIA in 1850)
    SIA_max = SIA_max_param_1850(j,1:12);


    % Starting test parameters
    a = test_input_vals(1);
    sm_shift = test_input_vals(2);
    w1_store = test_input_vals(3);
    b = test_input_vals(4);


    % Testable matrix
    p0 = [a,b,sm_shift,w1_store];
    params_hold = p0; 


    % Extraxt dec value before start date
    tas_hold = x_tas_store(end-length(x_tas_store)+1,:);
    dec_previous_year_tas = tas_hold(12);


    f = @(p)SIA_CMIP6_Calibration_Publication(SIA_max,x_tas,SIA,dec_previous_year_tas,p);       % Calibration function (function defined in a seperate script - SIA_CMIP6_Calibration_Publication.m)


    % Fminsearch
    opt_iter_options = optimset('Display','iter','TolX',1e-6,'MaxIter', 10e9);
    [p,fval] = fminsearch(f, p0, opt_iter_options);


    % Assign calibrated parameterss
    a = p(1);
    b = p(2);
    sm_shift = p(3);
    w1_store = p(4);



    % Weightings Optimisation
    w1_store_boundaries = [0, 1];
    if w1_store > w1_store_boundaries(2)
        [~,index_w1] = max(w1_store_boundaries - w1_store);
        w1_store = w1_store_boundaries(index_w1);
    elseif w1_store < w1_store_boundaries(1)
        [~,index_w1] = min(w1_store_boundaries - w1_store);
        w1_store = w1_store_boundaries(index_w1);
    end

    % Save parameters
    Save_sia_cal_params{j} = [a, sm_shift, w1_store, b];

%   Create empty vectors to store the new tas after weightings added
    x_tas_hold = [];
    dec_tas_test = x_tas(:,12);
    dec_tas_test = circshift(dec_tas_test,1);
    x_tas_shift = [dec_tas_test, x_tas];
    x_tas_shift(1) = dec_previous_year_tas;

    for tt = 1:12 

        if tt == 1
            hold_x_tas = ((x_tas(:,tt) .* w1_store) + (x_tas_shift(:,1) .* (1-w1_store))) ./ (w1_store+(1-w1_store));
        else
            hold_x_tas = ((x_tas(:,tt) .* w1_store) + (x_tas(:,tt-1) .* (1-w1_store))) ./ (w1_store+(1-w1_store));
        end

        x_tas_hold(:,tt) = hold_x_tas;
    end

    x_tas = x_tas_hold;     % Updated AMST (The Dependency of SIA on Temperature - section2.2.3)


    % Run the SIA parameterisation with optimised parameters
    calibrated_SIA = ( ((SIA_max+sm_shift) .* (1 + exp((x_tas(1,:)-b)))) - (a.*(x_tas - x_tas(1,:))) ) ./ ( 1 + exp(x_tas-b) );
    calibrated_SIA_store{j,n} = calibrated_SIA;
end

Save_sia_cal_params = cell2mat(Save_sia_cal_params');







%% Reshape SIA for plotting

tic
sia_cal_rearrange = [];
sia_cal_rearrange_final = [];
sia_cmip6_rearrange = [];
sia_cmip6_rearrange_final = [];
for n = 1:3
    disp('__________')    
    for j = 1:12
        disp('NEW MODEL')
        
        for i = 1:12
            sia_cal_rearrange{j,i} = calibrated_SIA_store{j,n}(:,i)';
            sia_cmip6_rearrange{j,i} = updated_hist_sia_annual_curve_all_models{j,n}(:,i)';
        end
    end
    sia_cal_rearrange_final{n} = sia_cal_rearrange;
    sia_cmip6_rearrange_final{n} = sia_cmip6_rearrange;
end






%% Figure 5: Test SIA Calibration

% Initialise Plot
index_run = [6,8,11,10];    % Example models
counter = 0;

hnew3 = [];

% Plot all models into a subplot
close; clc
figure(101)
set(gcf, 'Units', 'Inches', 'Position', [.5 .5 20 11])

for n = [1:3] % SSP
    for j = index_run % model

    %   Counter for plotting
        counter = counter + 1;
        h = subplot(3,4,counter);

        for nn = 1:12 % month

            % SIA           
            SIA_CMIP6 = sia_cmip6_rearrange_final{n}{j,nn};
            SIA_emulated = sia_cal_rearrange_final{n}{j,nn};

            % AMST
            AMT_CMIP6 = tas_cmip6_rearrange_final{n}{j,nn};
            AMT_emulated = tas_calibration_rearrange_final{n}{j,nn};


            % moving means
            AMT_CMIP6 = movmean(AMT_CMIP6, 20);
            AMT_emulated = movmean(AMT_emulated, 20);

            SIA_CMIP6 = movmean(SIA_CMIP6, 20);
            SIA_emulated = movmean(SIA_emulated, 20);


    %       Plot SIA vs AMST and compare emulation to CMIP6
            ss = plot(AMT_emulated, SIA_emulated, 'color', colorblind(nn,:), 'linewidth', 1.5);    % emulated
            hold on
            tt = plot(AMT_CMIP6, SIA_CMIP6, '-.', 'color', colorblind(nn,:), 'linewidth', 1.5);    % cmip6


            % Store plots for legend
            if counter == 1 && n == 1
                hnew3 = cat(1, hnew3, ss);
            end
            hold on

        end

        % Y Label
        if ismember(counter, [1])
            ylabel('Emulated SIA (million km^{2})', 'Position', [-39 -15.0000 -1])
        end
        % X Label
        if ismember(counter, [9])
            xlabel('Arctic Temperature (\circC)', 'Position', [68.0000 -5.0351 -1])
        end


        % Add SSP labels to each panel
        if ismember(counter, [1:12])
            if n == 1
                x = -28;
                y = 1.3;
                c = colmat(n,:);
                t1 = text(x, y, 'SSP-5.85', 'color', c, 'FontSize', 18, 'FontWeight', 'bold');

            elseif n == 2
                x = -28;
                y = 1.3;
                c = colmat(n,:);
                t2 = text(x, y, 'SSP-2.45', 'color', c, 'FontSize', 18, 'FontWeight', 'bold');

            elseif n == 3
                x = -28;
                y = 1.3;
                c = colmat(n,:);
                t3 = text(x, y, 'SSP-1.26', 'color', c, 'FontSize', 18, 'FontWeight', 'bold');
            end
        end



        % Re-position panels to fit the legend
        if ismember(counter, [1,5,9])
            pos = get(h, 'Position');
            posnew = pos; 
            posnew(1) = posnew(1) - 0.04; 
            set(h, 'Position', posnew);
        elseif ismember(counter, [2,6,10])
            pos = get(h, 'Position');
            posnew = pos; 
            posnew(1) = posnew(1) - 0.05; 
            set(h, 'Position', posnew);
        elseif ismember(counter, [3,7,11])
            pos = get(h, 'Position');
            posnew = pos; 
            posnew(1) = posnew(1) - 0.06; 
            set(h, 'Position', posnew);
        elseif ismember(counter, [4,8,12])
            pos = get(h, 'Position');
            posnew = pos; 
            posnew(1) = posnew(1) - 0.07; 
            set(h, 'Position', posnew);
        end

        % Set size of axes    
        set(gca, 'fontsize', 20)


    %   Add model title to each panel
        if n == 1
            sizza = 22;
            title([num2str(models_alone{j})],  'fontsize', sizza, 'interpreter', 'none')
        end

        % Add a grid to each plot
        grid

        % Axis limits
        xlim([-30 12])

    end
end    

% Get the months for plotting   
index_month_label = month_label(1:12);



% Plot legend outside of subplots
[tt2,hObj] = legend(hnew3, index_month_label, 'location', 'eastoutside');         
hL = findobj(hObj,'type','line'); 
set(hL,'linewidth', 4)


% Posistion legend in the correct place
newPosition = [0.86 0.45 0.1 0.1];
newUnits = 'normalized';
set(tt2,'Position', newPosition,'Units', newUnits);


% % Save
% temp=[' SIA Calibrated Models ', '.png']; 
% saveas(gca,temp);

