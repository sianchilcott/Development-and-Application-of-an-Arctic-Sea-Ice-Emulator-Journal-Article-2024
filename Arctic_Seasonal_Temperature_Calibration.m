%% Arctic Seasonal Temperature (AMST) Calibration
% Chapter 2 Section 2.3.3

% Initialise
parameters_store = [];      % Save parameters
calibrated_monthlytemp_save = [];   % Save optimised AMST to plot and compare to CMIP6

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
    calibrated_monthlytemp = f.*(cos(x .* g - e .* exp(cos(x.^(a))))+((AAMST - mean(f.*(cos(x .* g - e .* exp(cos(x.^(a))))), 2, 'omitnan')) ./ f));


    calibrated_monthlytemp_save{j} = calibrated_monthlytemp;


    % Save the calibration parameters    
    parameters_store{j} = [f1,f2,g1,g2,a1,a2,a3];


end
parameters_store = cell2mat(parameters_store');
parameters_store = parameters_store;
toc
