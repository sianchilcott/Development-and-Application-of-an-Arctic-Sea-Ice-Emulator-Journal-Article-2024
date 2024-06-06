%% Arctic Seasonal Temperature (AMST) Calibration
% Chapter 2 Section 2.3.3

% Initialise
parameters_store = [];      % Save parameters
AMST_calibration = [];   % Save optimised AMST to plot and compare to CMIP6

AAT_tas_store = [emul_absolute_AAT_585, emul_absolute_AAT_245, emul_absolute_AAT_126];  %  Get the calibrated Arctic annual mean temperature into a vector to put into seasonal temperature parameterisation 
x = repmat(linspace(0,2*pi,13), [length(year_ind).*3, 1]);  % initialise x for the cosine parameterisation
year = 1850:2100;       % Years of calibration
counter = 0;
index_run = 1:12;     % number of models to calibrate
year_ind = 1:251;  % Years to calibrate over


tic
% Calibration routine AMST
for j = index_run


%   Counter for plotting
    counter = counter + 1;


    tas = cell2mat(tas_store(:,j));     % extract the CMIP6 AMST of SSP5-8.5, SSP2-4.5, SSP1-2.6
    add_col = tas(2:end,1);             % make the annual temperature curve Jan-Jan of following year
    add_col = [add_col; nan];
    tas2 = [tas, add_col];
    y = tas2;


    y = tas2([year_ind, year_ind+251, year_ind+501],:);     % extract the CMIP6 AMST of the years to be calibrated (all SSPs are optimised together hence +251)

    AAMST = cell2mat(AAT_emulation_absolute(j,:))';                  % Get the emulated AAST from all SSP scenarios used to force the AMST parameterisation
    AAMST = AAMST([year_ind, year_ind+251, year_ind+501]);
    AAMST = repmat(AAMST, [1 13]);                          % Make it the same size as the CMIP6 annual temperature curve we are calibrating to



    % Calibration routine
    e = 0.3;
    ffunc = @(ppp)AMST_parameterisation_publication2(ppp,x,y,AAMST,e);    % Calibration function (in another script: AMST_parameterisation_publication2.m)


    test_params = [0.5016, -7.1394,   -0.01, 0.9262,  -0.0883, -0.4367, +0.1015];    % Initial inputs to start calibration


    % All the calibration parameters to be calibrated
    f1 = test_params(1);
    f2 = test_params(2);
    g1 = test_params(3);
    g2 = test_params(4);
    a1 = test_params(5);
    a2 = test_params(6);
    a3 = test_params(7);


    % Testable matrix
    p0 = [f1,f2,g1,g2,a1,a2,a3];
    params_hold = p0; 



    % Calibration
    opt_iter_options = optimset('Display','iter','TolX', 1e-6, 'MaxIter', 10e19);
    [ppp,fval] = fminsearch(ffunc, p0, opt_iter_options);



    % Assign calibrated variables to each parameters
    f1 = ppp(1);
    f2 = ppp(2);
    g1 = ppp(3);
    g2 = ppp(4);
    a1 = ppp(5);
    a2 = ppp(6);
    a3 = ppp(7);


    % Multiple the calibrated regression coefficients by the AAST
    f = (f1.*AAMST) + f2;
    g = (g1.*AAMST) + g2;
    a = cos((AAMST.*a1)+a2)+a3; 


    % Cap the linear regressions of the calibration parameters so they don't change after the AAST rises above 8 degrees celsius (explained Chapter 2 Section 2.3.6)
    ind_temp = 8;
    indf = AAMST > ind_temp;
    indf = find(indf==1);
    f(indf) = (f1.*ind_temp) + f2;
    g(indf) = (g1.*ind_temp) + g2;
    a(indf) = cos((ind_temp.* a1)+a2)+a3;



    % Put the calibrated parameters into the Arctic annual temperature cycle parameterisation
    AMST_calibration_hold = f.*(cos(x .* g - e .* exp(cos(x.^(a))))+((AAMST - mean(f.*(cos(x .* g - e .* exp(cos(x.^(a))))), 2, 'omitnan')) ./ f));


    AMST_calibration{j} = AMST_calibration_hold(1:12); % Extract only Jan to Dec of each year


    % Save the calibration parameters    
    parameters_store{j} = [f1,f2,g1,g2,a1,a2,a3];


end
parameters_store = cell2mat(parameters_store');
parameters_store = parameters_store;
toc


% Apply bias correction to correct for the pre-industrial offset
calls_for_mean_calibrated = [];
calls_for_mean_cmip6 = [];
AMST_calibration_bias_corrected = [];

mean_correction_year_index = 50;
for n = 1:3
    for j = 1:12
        BC_cmip6 = tas_store{1}{j}(1:mean_correction_year_index,:);     % Extract the 1850:1900 (pre-industrial) temperature 
        calls_for_mean_cmip6 = cat(1, calls_for_mean_cmip6, BC_cmip6);  % concatenate pre-industrial tas from all models

        calls = AMST_calibration{j,1}(1:mean_correction_year_index,:);  % Extract our emulated 1850:1900 (pre-industrial) temperature
        calls_for_mean_calibrated = cat(1, calls_for_mean_calibrated, calls);                 % concatenate temperatures for comparison with CMIP6

    end
    calls_for_mean_calibrated = mean(calls_for_mean_calibrated);      % Calculate the pre-industrial mean
    calls_for_mean_cmip6 = mean(calls_for_mean_cmip6);
    residuals = calls_for_mean_cmip6 - calls_for_mean_calibrated;     % Subtract the mean CMIP6 pre-industrial from our emulation to generate residuals

    for j = 1:12
        AMST_calibration_bias_corrected{j,n} = AMST_calibration{j,n} + residuals;       % Add residuals to our calibrated AMST
    end
end
tas_cal_new_bias_corrected = AMST_calibration_bias_corrected;



% Re-arrange the shape of our emulated AMST for plotting
tic
tas_calibration_rearrange = [];
tas_cmip6_calibration_rearrange = [];
tas_cmip6_rearrange_final = [];
tas_calibration_rearrange_final = [];
for n = 1:3
    disp('__________')    
    for j = 1:12
        disp('NEW MODEL')

        for i = 1:12
            tas_emul_rearrange{j,i} = AMST_calibration_bias_corrected{j,n}(:,i)';
            tas_cmip6_calibration_rearrange{j,i} = tas_store{n}{j}(:,i)';
        end
    end
    tas_calibration_rearrange_final{n} = tas_emul_rearrange;
    tas_cmip6_rearrange_final{n} = tas_cmip6_calibration_rearrange;
end
toc






%% Figure 2.3: Test calibration

% Plot Initialisation
index_run = [1, 4, 7, 12];  % Example models to plot
years_examples = [1, 151, 201, 251];   % Example years to plot (1850, 200, 2050, 2100)
years_labels = 1850:2100;   

colmat2 = [0, 0, 1.0000; 1.0000, 0, 0; 0, 0.3906, 0; 0.5977, 0.1953, 0.7969];
counter = 0;
hnew3 = [];
hnew3_cmip6 = [];


close
figure(101)
set(gcf, 'Units', 'Inches', 'Position', [.5 .5 20 11])
for n = [1:3]     % n = SSP scenario (n=1 is SSP5-8.5, n=2 is SSP2-4.5, n=3 is SSP1-2.6)
    for j = index_run  % Models 
        
        counter = counter + 1;
        h = subplot(3,4,counter);


        for nn = 1:length(years_examples)

            % Extract years to plot form calibrated temperature
            AMT_emulated = AMST_calibration_bias_corrected{j,n}(years_examples(nn),:);


            % Extract years to plot form CMIP6 temperature
            AMT_CMIP6 = tas_store{n}{j}(years_examples(nn),:);


            %  Plot emulated AMST
            ss = plot(1:12, AMT_emulated, 'color', colmat2(nn,:), 'linewidth', 2);  % Emulation
            hold on
            %  Plot CMIP6 AMST
            tt = plot(1:12, AMT_CMIP6, '-.', 'color', colmat2(nn,:), 'linewidth', 2);     % CMIP6


            % Store plots for legend
            if counter == 1 && n == 1
                hnew3 = cat(1, hnew3, ss);
                hnew3_cmip6 = cat(1, hnew3_cmip6, tt);
            end
            hold on

        end

        % Y Label
        if ismember(counter, [1])
            ylabel('Emulated Arctic Monthly Temperature (\circC)', 'Position', [-2.1947 -85.0000 -1])
        end
        % X Label
        if ismember(counter, [9])
            xlabel('Month', 'Position', [27.0000 -52.5819 -1.0000])
        end



        % Add SSP labels to each panel
        clear text
        if ismember(counter, [1:12])
            if n == 1
                x = 5;
                y = -36;
                c = colmat(n,:);
                t1 = text(x, y, 'SSP-5.85', 'color', c, 'FontSize', 20, 'FontWeight', 'bold');
            elseif n == 2
                x = 5;
                y = -36;
                c = colmat(n,:);
                t2 = text(x, y, 'SSP-2.45', 'color', c, 'FontSize', 20, 'FontWeight', 'bold');
            elseif n == 3
                x = 5;
                y = -36;
                c = colmat(n,:);
                t3 = text(x, y, 'SSP-1.26', 'color', c, 'FontSize', 20, 'FontWeight', 'bold');
            end
        end



        % Re-position each panel to fit the legend
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

        % Title
        if n == 1
            sizza = 22;
            title([num2str(models_alone{j})],  'fontsize', sizza, 'interpreter', 'none')
        end


        % Add a grid to each plot
        grid


        % Set axis limits
        ylim([-42 9])
        xlim([1 12])


        % Label the x-axis with the months of the year
        set(gca, 'xtick', 1:12);
        tick_length = 1;
        set(gca,'xticklabel', {[ blanks(tick_length) 'J'], [ blanks(tick_length) 'F'], [ blanks(tick_length) 'M'],...
            [ blanks(tick_length) 'A'],[ blanks(tick_length) 'M'],[ blanks(tick_length) 'J'],[ blanks(tick_length) 'J'],...
            [ blanks(tick_length) 'A'], [ blanks(tick_length) 'S'], [ blanks(tick_length) 'O'],[ blanks(tick_length) 'N'], ... 
        [ blanks(tick_length) 'D'], ''});
        xtickangle(90)

    end
end  

% Get the month names for legend   
index_month_label = string(years_labels(years_examples));

index_month_label_cmip6 = [];
index_month_label_emulated = [];
for i = 1:length(index_month_label)
    index_month_label_emulated_hold = strcat(index_month_label(i), ': Emulated'); 
    index_month_label_emulated = cat(1, index_month_label_emulated, index_month_label_emulated_hold);
    
    index_month_label_cmip6_hold = strcat(index_month_label(i), ': CMIP6'); 
    index_month_label_cmip6 = cat(1, index_month_label_cmip6, index_month_label_cmip6_hold);
end


% Plot legend outside of subplots
[tt2,hObj] = legend([hnew3; hnew3_cmip6], [index_month_label_emulated; index_month_label_cmip6] , 'location', 'eastoutside');         
hL = findobj(hObj,'type','line'); 
set(hL,'linewidth', 3)
tt2.Box = 'off';



% Posistion legend in the correct place
newPosition = [0.86 0.46 0.1 0.1]; 
newUnits = 'normalized';
set(tt2,'Position', newPosition,'Units', newUnits);


% % Save
% temp=[' AMT Calibrated Models, all SSPs 3 ', '.png']; 
% saveas(gca,temp);


