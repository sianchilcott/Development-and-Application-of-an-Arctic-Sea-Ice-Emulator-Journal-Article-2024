%% Stage iii: SIA Parameterisation
% SIA_max Calibration (Chapter 2 Section 2.3.4)


SIA_max_param_1850 = [];
SIA_max_param_1850_test = [];
counter = 0;
x = linspace(0,2*pi,13);    % Parameterisation is a cosine curve so we define the appropriate x value to go with this
index_run = 1:12;      % Models to calibrate over
plot_run = 1;

clc; close all
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
        f = 5;
        a = 1;
        e = 5;
        d = 14;
        g = 0.7;



        % Testable matrix
        p0 = [f,a,e,d,g];
        params_hold = p0; 


        % The model calibration routine for SIA, Arctic seasonal and SIA_max utilises the Nelder and Mead simplex optimisation method (Nelder and Meadf, 1965; Lagarias et al., 1998)
        opt_iter_options = optimset('Display','iter','TolX',1e-3,'MaxIter', 10e9); 
        [ppp, fval] = fminsearch(ff, p0, opt_iter_options);


        % Assign calibrated parameters
        f = ppp(1);
        a = ppp(2);
        e = ppp(3);
        d = ppp(4);
        g = ppp(5);



        % Input the calibrated parameters into the function and run
        e = (a .* 0.8003) + 1.5016;
        y_SIA_max_test = f .* (-exp(sin((x.^g) .* a - e))) + d; 



        % Save optimised parameters to put into SIA parameterisation
        SIA_max_param_1850 = cat(1, SIA_max_param_1850, y_SIA_max_test); 



        % Plot seasonal SIA in 1850 using calibration parameters to check their match the CMIP6 data 
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


            %   Add legend
            legend('Emulation', 'CMIP6 Starting SIA', 'location', 'southwest')
        end

%         % Save plots if needed
%         temp=[char(final_models_monthly_r1i1p1f1{j}), ' SIA_max 1850 ', '.png']; 
%         saveas(gca,temp);
    end
end




















