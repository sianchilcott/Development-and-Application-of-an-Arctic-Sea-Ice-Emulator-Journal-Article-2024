% ________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
%% Stage iii: SIA Parameterisation
%% Calibration: SIA_max 1.0


params_sia_max = [];
counter = 0;
x = linspace(0,2*pi,13);    % Parameterisation is a cosine curve so we define the appropriate x value to go with this
mov_mean_ind = 30;
SIA_max_param_1850 = [];
index_run = [1:12];      % Models to calibrate over
plot_run = 1;
SIA_max_param_1850_test = [];
params_save = [];

clc
close all
for n = 1    % We only calibrate for 1 SSP as all SSPs have the same historical emissions and therefore the 1850 annual SIA cycle is the same for all emission scenarios
    for j = index_run


    %   Counter for plotting
        counter = counter + 1;


    %   CMIP6 1850 SIA
        SIA_max_cmip6_running_mean = updated_hist_sia_annual_curve_all_models{j,n};
        SIA_max_cmip6 = SIA_max_cmip6_running_mean(1,:);
        jan_of_next_year = SIA_max_cmip6_running_mean(2,1);
        SIA_max_cmip6 = [SIA_max_cmip6, jan_of_next_year];      % Extract the SIA annual cycle from Jan-Jan of the following year



    %   Optimisation function
        ff = @(ppp)SIA_max_Calibration_Publication(ppp, x, SIA_max_cmip6);      % Calibration function defined in another script (SIA_max_Calibration_Publication.m)

        % Starting test parameters
        % f = 5;
        % a = 1;
        % e = 2.5;
        % d = 3;
        % g = 0.5;

        % f = 15;
        % a = 1;
        % e = 2;
        % d = 4;
        % g = 0.7;

                
        f = 5;
        a = 1;
        e = 5;
        d = 14;
        g = 0.7;



        % Testable matrix
        p0 = [f,a,e,d,g];
        params_hold = p0; 


        % Fminsearch
        opt_iter_options = optimset('Display','iter','TolX',1e-3,'MaxIter', 10e9);  % run calibration function
        [ppp, fval] = fminsearch(ff, p0, opt_iter_options);


        % Assign calibrated parameters
        f = ppp(1);
        a = ppp(2);
        e = ppp(3);
        d = ppp(4);
        g = ppp(5);



        % f = 4.3827;
        % a = 0.9242;
        % e = 2.0450;
        % d = 3.7568;
        % g = 0.6628;

        % g = 1;
        % a = 1;
        % e = 1;

    %   Run function with opimtisation parameters
        % e = 1;
        % y_SIA_max = f .* (-exp(sin(x .* a - e))) + sin(x.^g) + (f*d); 
        y_SIA_max = f .* (-exp(sin((x.^g) .* a - e))) + d; 

        % e = (a .* 2.9403) + -0.7302;
        e = (a .* 0.8003) + 1.5016;
        % e=0;
        % a = 0.5;
        y_SIA_max_test = f .* (-exp(sin((x.^g) .* a - e))) + d; 

        % SIA_max_param_1850_test = cat(1, SIA_max_param_1850_test, y_SIA_max);  
        SIA_max_param_1850 = cat(1, SIA_max_param_1850, y_SIA_max_test);  


        % Save calibrated parameters
        params_sia_max = cat(1, params_sia_max, [f,a,e,d,g]);



        % Test calibration for each model
        if plot_run == 1

            figure(2 + j)
            set(gcf, 'Units', 'Inches', 'Position', [.5 .5 20 11])

        %   Plot optimised values
            ss = plot(x, y_SIA_max, 'Color', [1 0 0], 'linewidth', 3);
            hold on
            plot(x, SIA_max_cmip6, '--', 'Color', 'k', 'linewidth', 3);
            plot(x, y_SIA_max_test, '-.', 'Color', 'b', 'linewidth', 3);


        %   Label X axis
            xlabel('Month')

        %   Label Y axis
            ylabel('SIA (km^2)')      

        %   Add title with calibrated parameters
            title(['SIA_max Emulation: ', '#', num2str(j), ', ', num2str(final_models_monthly_r1i1p1f1(j)), ', '], 'fontsize', 28, 'interpreter', 'none')

        %   Set the size of the plot font
            set(gca, 'fontsize', 22)

        %   X limts
            xlim([x(1) x(13)])       

        %   Add month names to the X label
            set(gca, 'xtick', x);
            set(gca,'xticklabel', {[ blanks(1) 'Jan'], [ blanks(1) 'Feb'], [ blanks(1) 'Mar'],[ blanks(1) 'Apr'],[ blanks(1) 'May'], ... 
            [ blanks(1) 'Jun'],[ blanks(1) 'Jul'], [ blanks(1) 'Aug'], [ blanks(1) 'Sep'], [ blanks(1) 'Oct'],[ blanks(1) 'Nov'], ... 
            [ blanks(1) 'Dec'], [ blanks(1) 'Jan'], [ blanks(1) 'Feb'], ''});


            legend('Emulation', 'CMIP6 Starting SIA', 'location', 'southwest')
        end

%         % Save
%         temp=[char(final_models_monthly_r1i1p1f1{j}), ' SIA_max 1850 ', '.png']; 
%         saveas(gca,temp);
    end
end


params_sia_max























%% Calibration: SIA_max 2.0


% params_sia_max = [];
% counter = 0;
% x = linspace(0,2*pi,13);    % Parameterisation is a cosine curve so we define the appropriate x value to go with this
% mov_mean_ind = 30;
% SIA_max_param_1850 = [];
% index_run = [1:4];      % Models to calibrate over
% plot_run = 1;
% SIA_max_param_1850_test = [];
% params_save = [];
% 
% clc
% close all
% for n = 1    % We only calibrate for 1 SSP as all SSPs have the same historical emissions and therefore the 1850 annual SIA cycle is the same for all emission scenarios
%     for j = index_run
% 
% 
%     %   Counter for plotting
%         counter = counter + 1;
% 
% 
%     %   CMIP6 1850 SIA
%         SIA_max_cmip6_running_mean = updated_hist_sia_annual_curve_all_models{j,n};
%         SIA_max_cmip6 = SIA_max_cmip6_running_mean(1,:);
%         jan_of_next_year = SIA_max_cmip6_running_mean(2,1);
%         SIA_max_cmip6 = [SIA_max_cmip6, jan_of_next_year];      % Extract the SIA annual cycle from Jan-Jan of the following year
% 
% 
% 
%     %   Optimisation function
%         ff = @(ppp)SIA_max_Calibration_Publication(ppp, x, SIA_max_cmip6);      % Calibration function defined in another script (SIA_max_Calibration_Publication.m)
% 
%         % Starting test parameters
%         f = 5;
%         a = 1;
%         d = 3;
%         g = 0.5;
% 
% 
%         % Testable matrix
%         p0 = [f,a,d,g];
%         params_hold = p0; 
% 
% 
%         % Fminsearch
%         opt_iter_options = optimset('Display','iter','TolX',1e-3,'MaxIter', 10e9);  % run calibration function
%         [ppp, fval] = fminsearch(ff, p0, opt_iter_options);
% 
% 
%         % Assign calibrated parameters
%         f = ppp(1);
%         a = ppp(2);
%         d = ppp(3);
%         g = ppp(4);
% 
% 
% 
%         % f = 4.3827;
%         % a = 0.9242;
%         % e = 2.0450;
%         % d = 3.7568;
%         % g = 0.6628;
% 
%         % g = 1;
%         % a = 1;
%         % e = 1;
% 
%     %   Run function with opimtisation parameters
%         y_SIA_max = f .* (-exp(sin(x .* a))) + sin(x.^g) + (f*d);
% 
%         % e = (a .* 2.9403) + -0.7302;
%         % e=1;
%         y_SIA_max_test = f .* (-exp(sin(x .* a))) + sin(x.^g) + (f*d);
% 
%         SIA_max_param_1850_test = cat(1, SIA_max_param_1850_test, y_SIA_max);
% 
% 
% 
%         % Save calibrated parameters
%         params_sia_max = cat(1, params_sia_max, [f,a,d,g]);
% 
% 
% 
%         % Test calibration for each model
%         if plot_run == 1
% 
%             figure(2 + j)
%             set(gcf, 'Units', 'Inches', 'Position', [.5 .5 20 11])
% 
%         %   Plot optimised values
%             ss = plot(x, y_SIA_max, 'Color', [1 0 0], 'linewidth', 3);
%             hold on
%             plot(x, SIA_max_cmip6, '--', 'Color', 'k', 'linewidth', 3);
%             plot(x, y_SIA_max_test, '-.', 'Color', 'b', 'linewidth', 3);
% 
% 
%         %   Label X axis
%             xlabel('Month')
% 
%         %   Label Y axis
%             ylabel('SIA (km^2)')      
% 
%         %   Add title with calibrated parameters
%             title(['SIA_max Emulation: ', '#', num2str(j), ', ', num2str(final_models_monthly_r1i1p1f1(j)), ', '], 'fontsize', 28, 'interpreter', 'none')
% 
%         %   Set the size of the plot font
%             set(gca, 'fontsize', 22)
% 
%         %   X limts
%             xlim([x(1) x(13)])       
% 
%         %   Add month names to the X label
%             set(gca, 'xtick', x);
%             set(gca,'xticklabel', {[ blanks(1) 'Jan'], [ blanks(1) 'Feb'], [ blanks(1) 'Mar'],[ blanks(1) 'Apr'],[ blanks(1) 'May'], ... 
%             [ blanks(1) 'Jun'],[ blanks(1) 'Jul'], [ blanks(1) 'Aug'], [ blanks(1) 'Sep'], [ blanks(1) 'Oct'],[ blanks(1) 'Nov'], ... 
%             [ blanks(1) 'Dec'], [ blanks(1) 'Jan'], [ blanks(1) 'Feb'], ''});
% 
% 
%             legend('Emulation', 'CMIP6 Starting SIA', 'location', 'southwest')
%         end
% 
% %         % Save
% %         temp=[char(final_models_monthly_r1i1p1f1{j}), ' SIA_max 1850 ', '.png']; 
% %         saveas(gca,temp);
%     end
% end
% 
% 
% params_sia_max

