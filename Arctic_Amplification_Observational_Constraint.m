%% Arctic Amplification: Section 2.3

% Initialise Obs

close; clc

% Label CMIP6 Calibrated Arctic Amplification (from CMIP6_Arctic_Amplification_Parameterisation.m)
CMIP6_AA = [];
for j = 1:12
    CMIP6_AA{j} = repelem(rc_save(j,1), 251);
end
CMIP6_AA = CMIP6_AA';

% Obs Mean
y1_CMIP6 = median(cell2mat(CMIP6_AA), 'omitnan');


%_____ Calculate the Observational Arctic Amplification from the observed Arctic annual mean temperature and global temperatures _____% (from xxx.m)
movmean_ind = 5;
global_mean_redshaped{4} = HADCRUT_GMST_obs{2}(1:171)';
AAT_OBS = [BERKELY_monthly_mean_AAT(2); nan; CRUTEM_monthly_mean_AAT(1); HADCRUT_AAT; nan; GISTEMP_monthly_mean_AAT(3)];


% OBS AA GISTEMP
length_years = 21;
x = movmean(GISTEMP_monthly_mean_AAT{3}(1:end-2), movmean_ind);
y = movmean(global_mean_redshaped{6}, movmean_ind);
AAT_slope = movingslope(x, length_years, 1);                        % movingslope is a downloaded added function from the matlab database (not my code)
GMST_slope = movingslope(y, length_years, 1);
obs_ensemble_AA_GISTEMP{1} = AAT_slope ./ GMST_slope;


% OBS AA BERKELY
x = movmean(BERKELY_monthly_mean_AAT{2}(1:171), movmean_ind);
y = movmean(global_mean_redshaped{1}, movmean_ind);
AAT_slope = movingslope(x, length_years, 1);
GMST_slope = movingslope(y, length_years, 1);
obs_ensemble_AA_BERKELY{1} = AAT_slope ./ GMST_slope;



% Obs Mean
AA_obs_final = [obs_ensemble_AA; obs_ensemble_AA_BERKELY; obs_ensemble_AA_GISTEMP];





%% Step 1: Create pdf from prior assumptions



% Set boundaries for pdfs created from prior assumptions
factor_p = [2.5, 3];      % parameter 'p' (Chapter 2, Section 2.4.2.1) - constant value after AA has increased over observational period
factor_s = [4, 6.5];                % parameter 's' (Chapter 2, Section 2.4.2.1) - slope of AA increase over observational period


%____ S (slope over obs period) ____%
uniform_AA_dist_factor_a = unifrnd(factor_s(1), factor_s(2), [1, 600])';
S = uniform_AA_dist_factor_a;


%____ P (constant after obs period) ____%
uniform_AA_dist_random_constant = pearsrnd(3.4,0.5,-0.75,3,600,1); 
P = uniform_AA_dist_random_constant;



% Plot prior assumption pdfs
close all
figure(42)
h = subplot(1,1,1);
set(gcf, 'Units', 'Inches', 'Position', [.5 .5 20 11])
h.Position(3) = h.Position(3) - 0.4;
h.Position(4) = h.Position(4) - 0.1;

% Plot P
[f,xi] = ksdensity(P); plot(xi,f, 'b', 'linewidth', 2); hold on

% Plot S
[f,xi] = ksdensity(S); plot(xi,f, 'r', 'linewidth', 2);


title('Prior Assumptions')
legend('P', 'S')
ylabel('Probability')
xlabel('P (AA factor) & S (AA factor change/ yr)')
grid

ax = gca; ax.FontSize = 18;

% % Save
% temp=[' AA obs constraint prior assump ', '.png']; 
% saveas(gca,temp);            





%% ___________________________________________________________%% 
%% Step 2: Create second pdf from models and obs: P



close all
[P_assumptions, xi_P] = ksdensity(P, 'NumPoints', 600);


% Take a gaussian over CMIP6_AA data to inform 'p'
x_P = cell2mat(CMIP6_AA);
std_P = std(x_P(:,1));      % standard devaition of data extract
mu_P = mean(x_P(:,1));      % mean of data extract
x_P = linspace(1, max(xi_P), 600); 
L_P = normpdf(x_P, mu_P, std_P); % Create a gaussian pdf over CMIP6 calibrated AA data (as it remains constant across simulation period we take the maxiumum, but this will be the same at any index)
prior_model_informed_P = L_P;








%%  Step 2: Create second pdf from models and obs: S



% Take a gaussian over observations to inform 's'
coefficients = [];
for Gi = 1:2
    xx = mean(GT_ensembles{1});
    for i = 1:length(AA_obs_final)
        coefficients_hold = polyfit(xx(121:171), AA_obs_final{i}(121:171), 1);  % Extract observations from 1970 to 2020 (as some datasets have data from 1850)
        coefficients{i, Gi} = coefficients_hold(1);
    end
end
coefficients = cell2mat(coefficients(:));
coefficients = coefficients .* 2;   % Get it into a form that works with our AAobs (Eq: 2.4) parameterisation



% S: Prior assumptions
[S_assumptions, xi_S] = ksdensity(S, 'NumPoints', 600);     % Prior assumption data



% S: Gaussian of observational uncertainty
pd = fitdist(coefficients, 'Kernel', 'Kernel', 'normal');
x_values_sia = linspace(min(coefficients), max(coefficients), length(coefficients));
pd_sia = pdf(pd, x_values_sia);          

std_S = std(coefficients, 'omitnan');
mu_S = mean(coefficients, 'omitnan');
x_S = linspace(1, max(xi_S), 600);
L_S = normpdf(x_S, mu_S, std_S);
prior_model_informed_S = L_S;
          





%% Merge Priors (average and normalise)



merged_prior_pdf_P  = (P_assumptions + prior_model_informed_P) ./ 2;
merged_prior_pdf_P = merged_prior_pdf_P ./ sum(merged_prior_pdf_P);
merged_prior_pdf_P = merged_prior_pdf_P .* 100;


merged_prior_pdf_S = (S_assumptions + prior_model_informed_S) ./ 2;
merged_prior_pdf_S = merged_prior_pdf_S ./ sum(merged_prior_pdf_S);
merged_prior_pdf_S = merged_prior_pdf_S .* 100;





%% Plot prior assumption PDFs and informed PDFs together

close all; clc
figure(42)
set(gcf, 'Units', 'Inches', 'Position', [.5 .5 20 11])
x_factor = 0;

% PLOT P: 2 priors
h = subplot(2,2,1);
h.Position(4) = h.Position(4) - x_factor;

plot(xi_P, P_assumptions, 'b', 'linewidth', 2); hold on            %  Plot prior assumptions pdf
plot(x_P, prior_model_informed_P, 'b--', 'linewidth', 2); hold on         %  Plot prior informed pdf

title('Prior pdf: P')
tt2 = legend('prior assumptions pdf', 'prior informed pdf', 'location', 'northwest');
tt2.Box = 'off';
ylabel('Probability')
xlabel('AA factor')
grid

ax = gca; ax.FontSize = 18;
ylim([0 1])




% PLOT S: 2 priors
h = subplot(2,2,2); 
h.Position(4) = h.Position(4) - x_factor;

plot(xi_S, S_assumptions, 'r', 'linewidth', 2); hold on                       %  Plot prior assumptions pdf
plot(x_S, prior_model_informed_S, 'r--', 'linewidth', 2); hold on           %  Plot prior informed pdf


title('Prior pdf: S')
tt2 = legend('prior assumptions pdf', 'prior informed pdf', 'location', 'northwest');
tt2.Box = 'off';
ylabel('Probability')
xlabel('AA factor change/ yr')
grid
ax = gca; ax.FontSize = 18; 
ylim([0 1])
xlim([1 7])




% PLOT P: Merged prior
h = subplot(2,2,3); 
h.Position(4) = h.Position(4) - x_factor;

plot(xi_P, merged_prior_pdf_P, 'k', 'linewidth', 2); hold on         %  Plot prior informed pdf
title('Final Prior pdf: P')
ylabel('Probability')
xlabel('AA factor')
grid

ax = gca; ax.FontSize = 18;
ylim([0 1])



% PLOT S: Merged prior
h = subplot(2,2,4);
h.Position(4) = h.Position(4) - x_factor;

plot(xi_S, merged_prior_pdf_S, 'k', 'linewidth', 2); hold on      
title('Final Prior pdf: S')
ylabel('Probability')
xlabel('AA factor change/ yr')
grid
ax = gca; ax.FontSize = 18; 
ylim([0 1])
% xlim([1 7])
xlim([1 8])



% % Save
% set(gcf, 'PaperOrientation', 'landscape')
% set(gcf,'PaperSize',[53 29.7000]);

% temp=[' AA merged prior pdf all 4 plots ', '.pdf']; 
% saveas(gca,temp); 








%% Extract values based on probability for input into AAobs (Eq:1)


num_samples = 600; % Number of samples you want to generate

% P
new_x_vector_P = []; 
for i = 1:length(xi_P)
    percentage_of_total = num_samples .* (round(merged_prior_pdf_P(i), 2));

    new_x_vector_hold = repelem(xi_P(i), round(percentage_of_total)); 

    new_x_vector_P = cat(2, new_x_vector_P, new_x_vector_hold);
end

xx = linspace(1, size(new_x_vector_P,2), 600);
xx = round(xx);
new_x_vector_P_new = new_x_vector_P(xx);
pd = fitdist(new_x_vector_P_new', 'Kernel', 'Kernel', 'normal');
x_values_sia = linspace(min(new_x_vector_P_new), max(new_x_vector_P_new), length(new_x_vector_P_new));
pd_sia = pdf(pd, x_values_sia);




% S
new_x_vector_S = []; 
for i = 1:length(xi_S)
    percentage_of_total = num_samples .* (round(merged_prior_pdf_S(i), 2));

    new_x_vector_hold = repelem(xi_S(i), round(percentage_of_total)); 

    new_x_vector_S = cat(2, new_x_vector_S, new_x_vector_hold);
end

xx = linspace(1, size(new_x_vector_S,2), 600);
xx = round(xx);
new_x_vector_S_new = new_x_vector_S(xx);
pd = fitdist(new_x_vector_S_new', 'Kernel', 'Kernel', 'normal');
x_values_sia = linspace(min(new_x_vector_S_new), max(new_x_vector_S_new), length(new_x_vector_S_new));
pd_sia = pdf(pd, x_values_sia);







%% Recalculate the Arctic Amplification with pdfs: 1: sample from pdfs using mean of GMST 


unifrm_extrction_num_P = randperm(600);
unifrm_extrction_num_S = randperm(600);
AA_MCMC = [];
GMST = [];
for n = 1:3
    % Generate random GMST ensembles array
    GMST = mean(GT_ensembles{n});

    % Randomly extract P
    P = new_x_vector_P_new(unifrm_extrction_num_P);

    % Randomly extract S
    S = new_x_vector_S_new(unifrm_extrction_num_S);


    % Use the randomly generate AA parameterisation values and randomly extracted GMST ensemble to calculate the AA 
    AA_MCMC{n} = P' ./ (1 + exp(S' .* (-GMST + 0.5))); 

    AA_MCMC{n}(:,1) = repelem(0, 600)';  % Make sure 1850 AA is 0 instead of NaN
end

           





%% Recalculate the AA with pdfs: 2: mean of p and s sample from GMST ensemble


unifrm_GMST = randperm(600);        % Randomly sample an ensemble from the MAGICC global mean temperature 600-member ensemble
AA_MCMC_2 = [];
GMST = [];
for n = 1:3
    % Generate random GMST ensembles array
    GMST{n} = GT_ensembles{n}(unifrm_GMST,:);

    % Randomly extract P
    P = mean(new_x_vector_P_new);           % Take mean of 'p' distribution

    % Randomly extract S
    S = mean(new_x_vector_S_new);           % Take mean of 's' distribution

    % Use the randomly generate AA parameterisation values and randomly extracted GMST ensemble to calculate the AA 
    AA_MCMC_2{n} = P' ./ (1 + exp(S' .* (-GMST{n} + 0.5))); 

    AA_MCMC_2{n}(:,1) = repelem(0, 600)'; 
end

            
GMST_random_arrangment2 = GMST;







%% Recalculate the AA with pdfs: 3: sample from p and s using mean of whole ensemble of GMST 


unifrm_GMST = randperm(600);                % Randomly sample an ensemble from the MAGICC global mean temperature 600-member ensemble
unifrm_extrction_num_P = randperm(600);     % Randomly sample a value from distribution of 'p'
unifrm_extrction_num_S = randperm(600);     % Randomly sample a value from distribution of 's'

AA_MCMC_3 = [];
save_GMST = [];
for n = 1:3
    % Generate random GMST ensembles array
    GMST = GT_ensembles{n}(unifrm_GMST,:);

    % Randomly extract P
    P = new_x_vector_P_new(unifrm_extrction_num_P);

    % Randomly extract S
    S = new_x_vector_S_new(unifrm_extrction_num_S);

    % Use the randomly generate AA parameterisation values and randomly extracted GMST ensemble to calculate the AA 
    AA_MCMC_3{n} = P' ./ (1 + exp(S' .* (-GMST + 0.5))); 

    AA_MCMC_3{n}(:,1) = repelem(0, 600)'; 

    save_GMST{n} = GMST;
end
            
GMST_random_arrangment3 = GMST;









%% PUBLICATION: ALL 4 PLOTS TOGETHER: WITH FINAL AA



close all; clc
figure(42)
set(gcf, 'Units', 'Inches', 'Position', [.5 .5 20 11])
x_factor = 0.07;

% PLOT P: Merged prior
h = subplot(2,2,3); 
h.Position(4) = h.Position(4) - x_factor;
h.Position(3) = h.Position(3) - (x_factor+0.05);
h.Position(2) = h.Position(2) -- x_factor;
h.Position(1) = h.Position(1) - x_factor;

plot(xi_P, merged_prior_pdf_P, 'k', 'linewidth', 2); hold on         %  Plot prior informed pdf
title('a) Final Distribution: P')
ylabel('Probability')
xlabel('AA factor')
grid

ax = gca; ax.FontSize = 18;
ylim([0 1])



% PLOT S: Merged prior
h = subplot(2,2,4);
h.Position(4) = h.Position(4) - x_factor;
h.Position(3) = h.Position(3) - (x_factor+0.05);
h.Position(2) = h.Position(2) -- x_factor;
h.Position(1) = h.Position(1) - (x_factor*3.5);

plot(xi_S, merged_prior_pdf_S, 'k', 'linewidth', 2); hold on      
title('b) Final Distribution: S')
ylabel('Probability')
xlabel('AA factor change/ yr')
grid
ax = gca; ax.FontSize = 18; 
ylim([0 1])
% xlim([1 7])
xlim([1 8])



% AA PLOT
axes('Position',[0.6 0.29 0.3 0.48])


y1_obs = median(cell2mat(AA_obs_final), 'omitnan');
y1_obs = y1_obs(1:end-8);
y1_obs = movmean(y1_obs, 2);


% X
x_options = 1850:2300;
x_CMIP6_label = 1850:2100;
x = 1850:2012;


threshold = [0 100];

% Percentile: Option 1
stdY_option_1 = prctile(AA_MCMC{1}, threshold);


% Percentile: CMIP6
stdY_CMIP6 = prctile(cell2mat(CMIP6_AA), threshold);


% Percentile: Obs
stdY_obs = prctile(cell2mat(AA_obs_final), threshold);
stdY_obs = stdY_obs(:, 1:end-8);


% Create curves: Option 1
curve1 = stdY_option_1(2,:);
curve2 = stdY_option_1(1,:);


% Create curves: Obs 
curve3 = stdY_obs(2,:);
curve4 = stdY_obs(1,:);


% Create curves: CMIP6
curve7 = stdY_CMIP6(2,:);
curve8 = stdY_CMIP6(1,:);



% Shading std area: Option 1
x_O1 = [x_options, fliplr(x_options)];
inBetween_O1 = [curve1, fliplr(curve2)];


% Shading std area: CMIP6
x_CMIP6 = [x_CMIP6_label, fliplr(x_CMIP6_label)];
inBetween_CMIP6 = [curve7, fliplr(curve8)];


% Shading std area: Obs 
x_OBS = [x, fliplr(x)];
inBetween2 = [curve3, fliplr(curve4)];


% Fill: Option 1
c = [0 0 1];
% c = [0.5 0.5 1];
f = fill(x_O1, inBetween_O1, c, 'edgecolor','none');
set(f,'facealpha',.3)
hold on


% Fill: CMIP6
c = [1 0 0];
f = fill(x_CMIP6, inBetween_CMIP6, c, 'edgecolor','none');
set(f,'facealpha',.3)
hold on


% Fill: Obs
d = [0 0 0]+0.05;
f2 = fill(x_OBS, inBetween2, d, 'edgecolor','none');
set(f2,'facealpha',.3)
hold on


% Plot emul and cmip6 means
plt_CMIP6 = plot(x_CMIP6_label, y1_CMIP6, 'color', [1 0 0], 'LineWidth', 2);
plt_OBS = plot(x, y1_obs, 'color', [0 0 0], 'LineWidth', 2);
% plt_MCMC = plot(1850:2300, median(save_AA_MCMC{1}), 'color', [0 0 1], 'LineWidth', 2);
plt_MCMC = plot(1850:2300, median(AA_MCMC{1}), 'color', [0 0 1], 'LineWidth', 2);



% Add CMIP6 AA
for n = 1:3
    for j = 1:12
        rw = tas_global{n}{j};
        raa = tas_arctic_annual{n}{j};

        AA = raa ./ rw;
        AA = movmean(AA, 20);

        x_cmip6 = 1850:2100;
        plt_cmip6_aa = plot(x_cmip6, AA, 'color', [0 0 0]+0.7, 'LineWidth', 1); hold on

        uistack(plt_cmip6_aa,'bottom'); 
    end
end

text(1960, 6.8, 'c)', 'Fontweight', 'bold', 'Fontsize', 20)




%   Y Label
ylabel({'Arctic Amplification'})
xlabel('Year')

xlim([1970 2100])
ylim([0 6])
grid


set(gca,'FontSize', 20)


axes('Position',[.902 .29 .1 .48])
box on

colorss_new = [1 0.8 0.8; 0.5 0.5 1; 0 0.6906 0.5];
histogram(AA_MCMC{1}(:,end), 'facecolor', colorss_new(2,:),  'Orientation','horizontal'); ylim([0 6])
set(gca, 'Color', 'None');
box off
ax = gca;
ax.YAxis.Visible = 'off';
axis off



x = [0.56 0.56];
y = [0.02 0.98];
annotation('line', x, y, 'LineStyle', '--', 'Linewidth', 2)


% Save
set(gcf, 'PaperOrientation', 'landscape')
set(gcf,'PaperSize',[53 29.7000]);

% % Save Normalised
% temp=['AA_with_pdfs_2', '.pdf']; 
% saveas(gca,temp);  













