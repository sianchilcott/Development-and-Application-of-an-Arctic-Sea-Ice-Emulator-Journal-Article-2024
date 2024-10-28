%% Appendix

%% Fig A.4


parameters_store = [0.4941   -8.4259   -0.0063    0.9309   -0.0780   -0.5773    0.1017
                        0.5613   -7.1602   -0.0071    0.9319   -0.0333    0.1384    0.1348
                        0.6322   -6.0756   -0.0031    0.9697   -0.0893   -0.6312    0.1086
                        0.5826   -6.2933   -0.0042    0.9557   -0.0964   -0.6432    0.1123
                        0.4154   -7.7789   -0.0007    0.9688   -0.0551   -0.0804    0.1267
                        0.5391   -7.4667   -0.0076    0.9505   -0.0875   -0.3211    0.0084
                        0.3465  -11.3899   -0.0009    1.0264   -0.0147    0.6405    0.3556
                        0.4858   -8.3263   -0.0062    0.9415   -0.0924   -0.5403    0.0691
                        0.4479   -9.0050   -0.0040    0.9766   -0.0973   -0.4534    0.0820
                        0.4378   -8.0920   -0.0057    0.9489   -0.1060   -0.3702    0.0921
                        0.5550   -6.4767   -0.0016    0.9877   -0.0218    0.9385    0.5679
                        0.5196   -7.3501   -0.0036    0.9572   -0.0943   -0.5025    0.0998];


f1 = parameters_store(j,1);
f2 = parameters_store(j,2);
g1 = parameters_store(j,3);
g2 = parameters_store(j,4);
a1 = parameters_store(j,5);
a2 = parameters_store(j,6);
a3 = parameters_store(j,7);
AAMST_hold = CMIP6_AMT_2300{5};

j=12;
counter = 0;
close all

f_hold = [];
g_hold = [];
a_hold = [];
h_hold = [];

for i = 1:length(AAMST_hold)
    counter = counter + 1;

    AAMST = AAMST_hold(i);

    f = (f1.*AAMST) + f2;
    g = (g1.*AAMST) + g2;
    e = 0.3;
    
    a = cos((AAMST.* a1) +a2) +a3; 

    ind_temp = 6;
    indf = AAMST > ind_temp;
    if indf == 1
        f(indf) = (f1.*ind_temp) + f2;
        g(indf) = (g1.*ind_temp) + g2;
        a(indf) = cos((ind_temp.* a1)+a2)+a3; 
    end

    starter_ind_temp = -15;
    indf = AAMST <= starter_ind_temp;
    f(indf) = (f1.*starter_ind_temp) + f2;

    
    x = linspace(0,2*pi,13);
    if ismember(counter, [9:11])
        c = [1 0 0];
    else
        c = [0 0 0];
    end
    calibrated_monthlytemp = f.*(cos(x .* g - e .* exp(cos(x.^a)))+((AAMST - mean(f.*(cos(x .* g - e .* exp(cos(x.^a)))), 2, 'omitnan')) ./ f));


    f_hold{i} = f;
    g_hold{i} = g;
    a_hold{i} = a;
    h_hold{i} = AAMST - mean(f.*(cos(x .* g - e .* exp(cos(x.^a)))));

end


% PLOT TO 2300
% make 1850:2100 black and 2100:2300 red
clc
close all
figure(3)
set(gcf, 'Units', 'Inches', 'Position', [.4 .4 19 10])
for i = 5
    subplot(1,2,1)
    for ii = 1:451
        if ismember(ii, 1:251)
            col = [0 0 0];
        else
            col = [1 0 0];
        end
        plot(CMIP6_AMT_2300{i}(ii,:), 'color', col, 'linewidth', 1); hold on
    end
    grid
    title(AMT_2300_models{i}, 'interpreter', 'none')
end
set(gca, 'xtick', 1:12);
tick_length = 1;
set(gca,'xticklabel', {[ blanks(tick_length) 'J'], [ blanks(tick_length) 'F'], [ blanks(tick_length) 'M'],...
    [ blanks(tick_length) 'A'],[ blanks(tick_length) 'M'],[ blanks(tick_length) 'J'],[ blanks(tick_length) 'J'],...
    [ blanks(tick_length) 'A'], [ blanks(tick_length) 'S'], [ blanks(tick_length) 'O'],[ blanks(tick_length) 'N'], ... 
[ blanks(tick_length) 'D'], ''});
ylabel('Arctic Temperature (\circC)')
xlabel('Month')
set(gca, 'Fontsize', 18)
xlim([1 12])



% get f_hold from Calobration_cmip6_test script
f_hold = cell2mat(f_hold);
g_hold = cell2mat(g_hold);
a_hold = cell2mat(a_hold);
h_hold = cell2mat(h_hold);


%  BOX 1 (f)
axes('Position',[.55 .59 .16 .33])
box on
scatter(mean(CMIP6_AMT_2300{5}, 2), f_hold, 'markerfacecolor', [0 0 0]+0.3, 'markeredgecolor', [0 0 0]+0.3); grid
xline(-8, 'linewidth', 1.5)
xline(8, 'linewidth', 1.5)
title('Calibration Parameter f')
xlabel('Arctic Annual Temperature (\circC)')
ylabel('Calibration Parameter f')
set(gca, 'Fontsize', 14)



%  BOX 2 (g)
axes('Position',[.78 .59 .16 .33])
box on
scatter(mean(CMIP6_AMT_2300{5}, 2), g_hold, 'markerfacecolor', [0 0 0]+0.3, 'markeredgecolor', [0 0 0]+0.3); grid
xline(8, 'linewidth', 1.5)
title('Calibration Parameter g')
xlabel('Arctic Annual Temperature (\circC)')
ylabel('Calibration Parameter g')
set(gca, 'Fontsize', 14)



%  BOX 3 (a)
axes('Position',[.66 .13 .16 .33])
box on
scatter(mean(CMIP6_AMT_2300{5}, 2), a_hold, 'markerfacecolor', [0 0 0]+0.3, 'markeredgecolor', [0 0 0]+0.3); grid
title('Calibration Parameter a')
xline(8, 'linewidth', 1.5)
xlabel('Arctic Annual Temperature (\circC)')
ylabel('Calibration Parameter a')
set(gca, 'Fontsize', 14)





% set(gcf, 'PaperOrientation', 'landscape')
% set(gcf,'PaperSize',[53 29.7000]);

% % Save
% temp=[' FigA.1 ', '.pdf']; 
% saveas(gca,temp);  









%% Fig A.10: Vary calibration parameters SIA

Save_sia_cal_params = new_vals;
cal_factors = {Save_sia_cal_params(:,1), Save_sia_cal_params(:,4), Save_sia_cal_params(:,3), Save_sia_cal_params(:,2)};
title_cfs = ["a", "b", "w1", "shift"];
 
month_index = 11;

close
for n = 2:2
    
    SIA_max = SIA_max_param_1850(1,1:12);
    
%   Vary each cf
    for j = 1:4
        %   Calibration Parameters to vary
        a = 0.2;
        b = 0;
        sm_shift = 0;
        w1 = 1;
    
        cf_to_vary = cal_factors{j}(n);
        
        figure(42)
        if ismember(j, 1:4)
            subplot(2, 2, j);
        else
            j_hold = j + 0.5;
            subplot(2, 2, (j+0.5));
        end
        set(gcf, 'Units', 'Inches', 'Position', [.5 .5 20 11])
        
        if j == 1
            vary_cf = [0:0.025:0.2];
            a = vary_cf;
        elseif j == 2
            vary_cf = [-3:0.25:0.5];
            b = vary_cf;
        elseif j == 3
            vary_cf = [1:-0.1:0];
            w1 = vary_cf;
        elseif j == 4
            vary_cf = [0:0.2:2];
            sm_shift = vary_cf;
        end
        
%       Length of test vectors
        for i = 1:length(vary_cf)

            x_tas = tas_store{n}{n};
            x_tas_orig = tas_store{n}{n};

            
           
            
            %   Calculate Weightings  
            if j == 3
                weightings_store = repmat([w1(i), 1-w1(i)], 12, 1);
           

                x_tas_hold = [];
                x_tas_diff = [];

                tt_store = 1:12;
                x_tas_shift = circshift(x_tas,1);
                x_tas_shift(1,:) = repelem(x_tas(1), 12);

                for tt = 1:12

                    weightings = weightings_store(tt,:); 
                    tt_hold = circshift(tt_store,2);

                    if tt == 1
                        hold_x_tas = ((x_tas(:,tt) .* weightings(1)) + (x_tas_shift(:,tt_hold(tt+1)) .* weightings(2))) ./ sum(weightings);
                    else
                        hold_x_tas = ((x_tas(:,tt) .* weightings(1)) + (x_tas(:,tt-1) .* weightings(2))) ./ sum(weightings);
                    end

                    x_tas_hold(:,tt) = hold_x_tas;
                end

                x_tas = x_tas_hold(:,month_index); 
            else
                x_tas = x_tas_orig(:,month_index); 
            end
            

            updated_fun2 = (((SIA_max(month_index)+sm_shift) .* (1 + exp(x_tas(1,:)-b))) - (a.*(x_tas - x_tas(1,:)))) ./ (1 + exp(x_tas-b));
            
            plot(x_tas_orig(:,month_index), updated_fun2, 'Color', [0 0 1], 'linewidth', 2);
            
            xx = 1;
            xx1 = x_tas_orig(1)+1;
            end_x1 = xx1;
            if j == 1
                end_x = -0.8;
                end_y = 9.8;
                
                end_x1 = -19.6;
                end_y1 = 6.0;               
                
                text(end_x1, end_y1, 'Default Parameters:', 'fontsize', 14)
                text_param = ["b = 0", "w1 = 1", "shift = 0"];
                text(end_x1, end_y1-1, text_param, 'fontsize', 14)
            elseif j == 2
                end_x = -0.6;
                end_y = 5.7;
                
                end_x1 = -19.6;
                end_y1 = 2.3;
                
                text(end_x1, end_y1, 'Default Parameters:', 'fontsize', 14)
                text_param = ["a = 0.2", "w1 = 1", "shift = 0"];
                text(end_x1, end_y1-1.2, text_param, 'fontsize', 14)
            elseif j == 3
                end_x = -0.7;
                end_y = 8.9;
                
                end_x1 = -19.6;
                end_y1 = 4.3;
                
                text(end_x1, end_y1, 'Default Parameters:', 'fontsize', 14)
                text_param = ["a = 0.2", "b = 0", "shift = 0"];
                text(end_x1, end_y1-1.2, text_param, 'fontsize', 14)
            elseif j == 4
                end_x = -0.5;
                end_y = 10.8;
                
                end_x1 = -19.6;
                end_y1 = 6.5;
                
                text(end_x1, end_y1, 'Default Parameters:', 'fontsize', 14)
                text_param = ["a = 0.2", "b = 0", "w1 = 1"];
                text(end_x1, end_y1-1.2, text_param, 'fontsize', 14)
            end
                
            text(end_x, end_y, string(vary_cf), 'fontsize', 12)
            hold on
        end
        sgtitle(['SIA Parameterisation: Varied Calibration Parameters'], 'fontsize', 30)
        title([title_cfs(j)])
        set(gca, 'fontsize', 22)
        xlim([-20, 1])

        
%       X and Y Label
        if ismember(j,[3,4])
            xlabel('November Arctic Temperature (\circC)')
        end
        if ismember(j,[1,3])            
            ylabel('SIA (million km^{2})')
        end
    end
end


% set(gcf, 'PaperOrientation', 'landscape')
% set(gcf,'PaperSize',[53 29.7000]);

% % Save
% temp=['FigA.4', '.pdf']; 
% saveas(gca,temp); 






%% Fig A.3: Vary calibration parameters AMT

xx2 = linspace(0,2*pi,13);
x = xx2';


AAT_tas_store = {emul_absolute_AAT_585; emul_absolute_AAT_245; emul_absolute_AAT_126};
cf = {'a', 'f', 'g', 'e'};

col_mat = hsv(5);
year_index = 12;

close
figure(42)
set(gcf, 'Units', 'Inches', 'Position', [.4 .4 20 11])

for n = 1:4

	hh = subplot(2, 2, n);

            
    for j = 12


    %   Counter for plotting
        counter = counter + 1;


        AAMST = AAT_emulation_absolute{j,1}(year_index);

        a = cos((AAMST.*-0.0728)-0.5531)+0.0741;


        kholda = parameters_store(j,:);
        f1 = kholda(1);
        f2 = kholda(2);
        g1 = kholda(3);
        g2 = kholda(4);
        e1 = kholda(5);
        e2 = kholda(6);

        f = (f1.*AAMST) + f2;
        g = (g1.*AAMST) + g2;
        e = (e1.*AAMST) + e2;



        if n == 1
            vary_cf = [0:0.05:1];
            a = vary_cf;

            f = -10;
            g = 1;
            e = 1;

        elseif n == 2
            vary_cf = [-15:2:5];
            f = vary_cf;

            a = 0.9;
            g=1;
            e=0;

        elseif n == 3
            vary_cf = [0.7:0.05:1];
            g = vary_cf;

            f = -10;
            a = 0.9;
            e = 0;

        elseif n == 4
            vary_cf = [0:0.1:1.5];
            e = vary_cf;

            f = -10;
            a = 0.9;
            g = 1;

        end


        monthlytemp2 = f.*(cos(x .* g - e .* exp(cos(x.^(a))))+((AAMST - mean(f.*(cos(x .* g - e .* exp(cos(x.^(a))))))) ./ f));

        ss = plot(x, monthlytemp2, 'color', [0 0 1], 'linewidth', 1);


        if x(13) == 13
            xlim([1 13])
            set(gca, 'xtick', 1:13);
        elseif x(13) == 12
            xlim([0 12])
            set(gca, 'xtick', 0:12);
        elseif x(13) == 1
            set(gca, 'xtick', [linspace(0,1-1/12,12),1]  );
        elseif x(13) == (2*pi)
            set(gca, 'xtick', [linspace(0,2*pi,13)]  );
            xlim([0 2*pi])
        end
        tick_length = 1;
        set(gca,'xticklabel', {[ blanks(tick_length) 'J'], [ blanks(tick_length) 'F'], [ blanks(tick_length) 'M'],...
            [ blanks(tick_length) 'A'],[ blanks(tick_length) 'M'],[ blanks(tick_length) 'J'],[ blanks(tick_length) 'J'],...
            [ blanks(tick_length) 'A'], [ blanks(tick_length) 'S'], [ blanks(tick_length) 'O'],[ blanks(tick_length) 'N'], ... 
        [ blanks(tick_length) 'D'], [ blanks(tick_length) 'J'], ''});
        xtickangle(0)

       


%           Length of test vectors
        font_siza = 15;
        if n == 1
            end_x1 = 6.4;
            end_y1 = -14;

            text_param = {'Default Parameters:', ['f=', num2str(f)], ['g=', num2str(g)], ['e=', num2str(e)]};
            text(end_x1, end_y1, text_param, 'fontsize', font_siza)

        elseif n == 2
            end_x1 = 6.4;
            end_y1 = -10;

            text_param = {'Default Parameters:', ['a=', num2str(a)], ['g=', num2str(g)], ['e=', num2str(e)]};
            text(end_x1, end_y1, text_param, 'fontsize', font_siza)

        elseif n == 3               
            end_x1 = 6.4;
            end_y1 = -12;

            text_param = {'Default Parameters:', ['a=', num2str(a)], ['f=', num2str(f)], ['e=', num2str(e)]};
            text(end_x1, end_y1, text_param, 'fontsize', font_siza)

        elseif n == 4
            end_x1 = 6.4;
            end_y1 = -15;

            text_param = {'Default Parameters:', ['a=', num2str(a)], ['f=', num2str(f)], ['g=', num2str(g)]};
            text(end_x1, end_y1, text_param, 'fontsize', font_siza)

        elseif n == 5
            end_x1 = x(end)+0.1;
            end_y1 = monthlytemp2(end)+0.15;

            text_param = {'Default Parameters:', ['a=', num2str(a)], ['f=', num2str(f)], ['g=', num2str(g)], ['e=', num2str(e)]};
            text(end_x1, end_y1, text_param, 'fontsize', font_siza)
        end

        sgtitle(['Arctic Monthly Temperature Parameterisation: Varied Calibration Parameters'], 'fontsize', 28)
        title([cf(n)])
        
        if ismember(n,[3,4]) 
        	xlabel('Month')
        end
        if ismember(n,[1,3])            
            ylabel({'December'; 'Arctic Temperature (\circC)'})
        end
        
        set(gca, 'fontsize', 22)


        if ismember(n, [1,3])
            pos = get(hh, 'Position');
            posnew = pos; 
            posnew(1) = posnew(1) - 0.04; 
            posnew(2) = posnew(2) - 0.03; 
            set(hh, 'Position', posnew);
        else
            pos = get(hh, 'Position');
            posnew = pos; 
            posnew(1) = posnew(1) - 0.02; 
            posnew(2) = posnew(2) - 0.03; 
            set(hh, 'Position', posnew);
        end
        
    end
end



% set(gcf, 'PaperOrientation', 'landscape')
% set(gcf,'PaperSize',[53 29.7000]);

% % Save
% temp=[' Fig A.2 ', '.pdf']; 
% saveas(gca,temp);








%% Calculate the gof: AMT


gof_store_final_AMT = [];
gof_store = [];
gof_store_hold_final = [];
gof_store_hold_final2 = [];

index_new = [1:12];
for n = 1
    gof_store_hold = [];
    for i = index_new

        AMT = amt_mon_rearrange_final_MCMC{n}(:, i);


        AMT_obs = monthly_obs_temp{i}';           
        AMT_obs = movmean(AMT_obs, 20);



        A = cellfun(@(v) (AMT_obs - v(:,1:171)).^2, AMT, 'un', 0);
        A2 = cellfun(@(v) sum(v, 2, 'omitnan'), A,  'un', 0);
        gof = cellfun(@(v) v ./ length(AMT_obs), A2, 'un', 0);
        gof = cellfun(@(v) mean(v), gof, 'UniformOutput', true);

        gof_store{i} = gof;

      
    end
    gof_store_final_AMT = gof_store;
end

gof_store_final_AMT_final = cellfun(@(v) mean(v), gof_store_final_AMT, 'UniformOutput', true);



% Get vars ready for table
A_LB_march = [ gof_store_final_AMT{3}; mean(gof_store_final_AMT{3})];
A_LB_sept = [ gof_store_final_AMT{9}; mean(gof_store_final_AMT{9})];




% Create first set of tables
A_LB = table(A_LB_march, A_LB_sept);



A_LB.Properties.VariableNames{'A_LB_march'} = 'March';
A_LB.Properties.VariableNames{'A_LB_sept'} = 'Sept';



model_names = [ MAGICC_ensemble_name_string; "Mean" ];

T1 = table(model_names, A_LB);
T1.Properties.VariableNames{'model_names'} = 'Model Names';

T1.Properties.VariableNames{'A_LB'} = 'AA LB: AMT GOF';


T1.(1) = categorical(T1.(1))







%% Calculate the gof: SIA


gof_store_final_SIA = [];
% gof_store = zeros(12,12);
gof_store = [];
gof_store_hold_final = [];
gof_store_hold_final2 = [];

index_new = [1:12];
for n = 1
    gof_store_hold = [];
    for i = index_new
        for m = 1:6

            SIA = sia_mon_rearrange_final_MCMC{n}(:, i);


            SIA_OBS = monthly_obs_SIA{m}{i};              
            SIA_OBS = movmean(SIA_OBS, 20);
            


            A = cellfun(@(v) (SIA_OBS - v(:,1:170)).^2, SIA, 'un', 0);
            A2 = cellfun(@(v) sum(v, 2, 'omitnan'), A,  'un', 0);
            gof = cellfun(@(v) v ./ length(SIA_OBS), A2, 'un', 0);
            gof = cellfun(@(v) mean(v), gof, 'UniformOutput', true);

            gof_store{m,i} = gof;

        end       
    end
    gof_store_final_SIA = gof_store;
end




% Get vars ready for table
A_LB_march = [ gof_store_final_SIA{3}(:,1); mean(gof_store_final_SIA{3}(:,1))];
A_LB_sept = [ gof_store_final_SIA{9}(:,1); mean(gof_store_final_SIA{9}(:,1))];




% Create first set of tables
A_LB = table(A_LB_march, A_LB_sept);



A_LB.Properties.VariableNames{'A_LB_march'} = 'March';
A_LB.Properties.VariableNames{'A_LB_sept'} = 'Sept';




model_names = [ MAGICC_ensemble_name_string; "Mean" ];

% T1 = table(model_names, A_LB, A_MB, A_UB);
T1 = table(model_names, A_LB);
T1.Properties.VariableNames{'model_names'} = 'Model Names';

T1.Properties.VariableNames{'A_LB'} = 'AA LB: SIA GOF';

T1.(1) = categorical(T1.(1))






%% Only Means in Table

mean_AMST_months = [];
mean_SIA_months = [];
for i = 1:12
    mean_SIA_months{i} = mean(cell2mat(gof_store_final_SIA(:,i)));
    mean_AMST_months{i} = mean(cell2mat(gof_store_final_AMT(:,i)));
end
mean_SIA_months = cell2mat(mean_SIA_months)';
mean_AMST_months = cell2mat(mean_AMST_months)';


% Create first set of tables
model_names = [string(month_label)' ];
T1 = table(model_names, mean_SIA_months, mean_AMST_months);
T1.Properties.VariableNames{'mean_SIA_months'} = 'RSS: SIA';
T1.Properties.VariableNames{'mean_AMST_months'} = 'RSS: AMST';
T1.Properties.VariableNames{'model_names'} = 'Month'; 
T1.(1) = categorical(T1.(1))






%% CMIP6: Calculate the gof: AMT

gof_store_final_AMT = [];
gof_store = zeros(12,12);

index_new = [1:12];
for n = 1:3
    for j = 1:12
        for i = index_new
            
            AMT_CMIP6 = tas_cmip6_rearrange_2{n}{j,i};
            AMT_emulated = tas_emululation_rearrange_NO_BC{n}{j,i};
            
            
            AMT_CMIP6 = movmean(AMT_CMIP6, 20);
            AMT_emulated = movmean(AMT_emulated, 20);
            
            
            A = (AMT_CMIP6 - AMT_emulated).^2;
            A = sum(A);
            gof = A ./ length(AMT_CMIP6);

            mdl = fitlm(AMT_CMIP6, AMT_emulated);
            r_squared = mdl.Rsquared.Ordinary;
            
            gof_store(j,i) = gof;
%             gof_store(j,i) = r_squared;

        end
    end
    gof_store_final_AMT{n} = gof_store;
end




% Get vars ready for table
A_585_march = [ gof_store_final_AMT{1}(:,3); mean(gof_store_final_AMT{1}(:,3))];
A_585_sept = [ gof_store_final_AMT{1}(:,9); mean(gof_store_final_AMT{1}(:,9))];


A_245_march = [ gof_store_final_AMT{2}(:,3); mean(gof_store_final_AMT{2}(:,3))];
A_245_sept = [ gof_store_final_AMT{2}(:,9); mean(gof_store_final_AMT{2}(:,9))];


A_126_march = [ gof_store_final_AMT{3}(:,3); mean(gof_store_final_AMT{3}(:,3))];
A_126_sept = [ gof_store_final_AMT{3}(:,9); mean(gof_store_final_AMT{3}(:,9))];


% Create first set of tables
A_585 = table(A_585_march, A_585_sept);
A_245 = table(A_245_march, A_245_sept);
A_126 = table(A_126_march, A_126_sept);

A_585.Properties.VariableNames{'A_585_march'} = 'March';
A_585.Properties.VariableNames{'A_585_sept'} = 'Sept';

A_245.Properties.VariableNames{'A_245_march'} = 'March';
A_245.Properties.VariableNames{'A_245_sept'} = 'Sept';

A_126.Properties.VariableNames{'A_126_march'} = 'March';
A_126.Properties.VariableNames{'A_126_sept'} = 'Sept';

% model_names = [ final_models_monthly_r1i1p1f1; "Mean" ];
model_names = [ models_alone; "Mean" ];

T1 = table(model_names, A_585, A_245, A_126);
T1.Properties.VariableNames{'model_names'} = 'Model Names';

T1.Properties.VariableNames{'A_585'} = 'SSP585: AMT GOF';
T1.Properties.VariableNames{'A_245'} = 'SSP245: AMT GOF';
T1.Properties.VariableNames{'A_126'} = 'SSP126: AMT GOF';

% T1.Properties.VariableNames{'A_585'} = 'SSP585: AMT R-Squared';
% T1.Properties.VariableNames{'A_245'} = 'SSP245: AMT R-Squared';
% T1.Properties.VariableNames{'A_126'} = 'SSP126: AMT R-Squared';


T1.(1) = categorical(T1.(1))




%% Table of parameters: SIA




model_names = [ models_alone; "Mean" ];


new_vals = [0.3795    0.3961    0.5408   -1.3240
            0.4006   -0.0411    0.5451   -0.9881
            0.2952    0.0687    0.2140   -4.4366
            0.2958   -0.1261    0.1251   -4.4320
            0.4096   -0.0108    1.0000   -0.7661
            0.4010    0.5436    0.5505    0.6585
            0.3473   -0.1438    0.4485    0.6471
            0.3242    0.4528    0.2728   -0.9961
            0.3240   -0.4721    0.4894   -1.1122
            0.2555    0.2036         0   -1.3647
            0.3889   -0.4409         0   -1.7137
            0.3848    0.0896    0.1904   -1.8352];


a = [new_vals(:,1); mean(new_vals(:,1))];
b = [new_vals(:,4); mean(new_vals(:,4))];
shift = [new_vals(:,2); mean(new_vals(:,2))];
w1 = [new_vals(:,3); mean(new_vals(:,3))];

T1 = table(model_names, a,b,shift,w1);

T1.Properties.VariableNames{'model_names'} = 'Model Names'




%% Table of parameters: AMST




model_names = [ models_alone; "Mean" ];

parameters_store_new = [0.4941   -8.4259   -0.0063    0.9309   -0.0780   -0.5773    0.1017
                        0.5613   -7.1602   -0.0071    0.9319   -0.0333    0.1384    0.1348
                        0.6322   -6.0756   -0.0031    0.9697   -0.0893   -0.6312    0.1086
                        0.5826   -6.2933   -0.0042    0.9557   -0.0964   -0.6432    0.1123
                        0.4154   -7.7789   -0.0007    0.9688   -0.0551   -0.0804    0.1267
                        0.5391   -7.4667   -0.0076    0.9505   -0.0875   -0.3211    0.0084
                        0.3465  -11.3899   -0.0009    1.0264   -0.0147    0.6405    0.3556
                        0.4858   -8.3263   -0.0062    0.9415   -0.0924   -0.5403    0.0691
                        0.4479   -9.0050   -0.0040    0.9766   -0.0973   -0.4534    0.0820
                        0.4378   -8.0920   -0.0057    0.9489   -0.1060   -0.3702    0.0921
                        0.5550   -6.4767   -0.0016    0.9877   -0.0218    0.9385    0.5679
                        0.5196   -7.3501   -0.0036    0.9572   -0.0943   -0.5025    0.0998];




f1 = [parameters_store_new(:,1); mean(parameters_store_new(:,1))];
f2 = [parameters_store_new(:,2); mean(parameters_store_new(:,2))];

g1 = [parameters_store_new(:,3); mean(parameters_store_new(:,3))];
g2 = [parameters_store_new(:,4); mean(parameters_store_new(:,4))];

a1 = [parameters_store_new(:,5); mean(parameters_store_new(:,5))];
a2 = [parameters_store_new(:,6); mean(parameters_store_new(:,6))];
a3 = [parameters_store_new(:,7); mean(parameters_store_new(:,7))];

T1 = table(model_names, f1,f2, g1,g2, a1,a2,a3);

T1.Properties.VariableNames{'model_names'} = 'Model Names'






%% 2300 MODELS TABLE


T = table(SIA_2300_models');
T.Properties.VariableNames{'Var1'} = 'Models with Data to 2300'






%% TABLE OF ALL EMULATORS 


emulator_names = ["CE"; "OCCAA"; "OC"]
emulator_description = ["CMIP6 calibrated emulator"; "Bias correction ONLY of global warming"; "Bias correction of both global warming and AA"];

T = table( char(emulator_names), char(emulator_description) )
T.Properties.VariableNames{'Var1'} = 'Emulator Names'
T.Properties.VariableNames{'Var2'} = 'Emulator Description'




%% Make a table of the means and std of each month


close all; clc

% AMST
p_save = zeros(1,12);
k_save = zeros(1,12);
for i = 1:12
    subplot(4,3,i)

    % Emulation
    pd = fitdist(AMT_emulated_2080_2100_MEAN(:,i), 'Kernel', 'Kernel', 'normal');
    x_values_emul = linspace(min(AMT_emulated_2080_2100_MEAN(:,i)), max(AMT_emulated_2080_2100_MEAN(:,i)), 100);
    pd_emul = pdf(pd, x_values_emul);
    scatter(x_values_emul, pd_emul); hold on
    
    
    
    % CMIP6
    pd = fitdist(AMT_CMIP6_2080_2100_MEAN(:,i), 'Kernel', 'Kernel', 'normal');
    x_values_cmip6 = linspace(min(AMT_CMIP6_2080_2100_MEAN(:,i)), max(AMT_CMIP6_2080_2100_MEAN(:,i)), 100);
    pd_cmip6 = pdf(pd, x_values_cmip6);
    scatter(x_values_cmip6, pd_cmip6)

    % disp(i)
    mean_emul = mean(AMT_emulated_2080_2100_MEAN(:,i));
    std_emul = std(AMT_emulated_2080_2100_MEAN(:,i));
    
    mean_cmip6 = mean(AMT_CMIP6_2080_2100_MEAN(:,i));
    std_cmip6 = std(AMT_CMIP6_2080_2100_MEAN(:,i));

    mean_save_cmip6(i) = mean_cmip6;
    mean_save_emul(i) = mean_emul;
    std_save_cmip6(i) = std_cmip6;
    std_save_emul(i) = std_emul;
    
end

% SIA
close all; clc

p_save = zeros(1,12);
k_save = zeros(1,12);
for i = 1:12
    subplot(4,3,i)

    % Emulation
    pd = fitdist(SIA_emulated_2080_2100_MEAN(:,i), 'Kernel', 'Kernel', 'normal');
    x_values_emul = linspace(min(SIA_emulated_2080_2100_MEAN(:,i)), max(SIA_emulated_2080_2100_MEAN(:,i)), 100);
    pd_emul = pdf(pd, x_values_emul);
    % scatter(x_values_emul, pd_emul); hold on
    
    
    
    % CMIP6
    pd = fitdist(SIA_CMIP6_2080_2100_MEAN(:,i), 'Kernel', 'Kernel', 'normal');
    x_values_cmip6 = linspace(min(SIA_CMIP6_2080_2100_MEAN(:,i)), max(SIA_CMIP6_2080_2100_MEAN(:,i)), 100);
    pd_cmip6 = pdf(pd, x_values_cmip6);
    % scatter(x_values_cmip6, pd_cmip6)

    % disp(i)
    mean_emul = mean(SIA_emulated_2080_2100_MEAN(:,i));
    std_emul = std(SIA_emulated_2080_2100_MEAN(:,i));
    
    mean_cmip6 = mean(SIA_CMIP6_2080_2100_MEAN(:,i));
    std_cmip6 = std(SIA_CMIP6_2080_2100_MEAN(:,i));

    mean_save_cmip6_sia(i) = mean_cmip6;
    mean_save_emul_sia(i) = mean_emul;
    std_save_cmip6_sia(i) = std_cmip6;
    std_save_emul_sia(i) = std_emul;
    
end



clc
% SIA
std_emulation = [ std_save_emul_sia'; mean(std_save_emul_sia)];
std_cmip6_models = [ std_save_cmip6_sia'; mean(std_save_cmip6_sia)];

mean_emulation = [ mean_save_emul_sia'; mean(mean_save_emul_sia)];
mean_cmip6_models = [ mean_save_cmip6_sia'; mean(mean_save_cmip6_sia)];


% Create first set of tables: std
A_LB = table(std_emulation, std_cmip6_models);
A_LB.Properties.VariableNames{'std_emulation'} = 'Emulation';
A_LB.Properties.VariableNames{'std_cmip6_models'} = 'CMIP6';


% Create first set of tables: mean
A_LB_2 = table(mean_emulation, mean_cmip6_models);
A_LB_2.Properties.VariableNames{'mean_emulation'} = 'Emulation';
A_LB_2.Properties.VariableNames{'mean_cmip6_models'} = 'CMIP6';


month_names = [ month_label' ; "Mean" ];

T1 = table(month_names, A_LB, A_LB_2);
T1.Properties.VariableNames{'month_names'} = 'Months';

T1.Properties.VariableNames{'A_LB'} = 'SIA: Standard Deviation';
T1.Properties.VariableNames{'A_LB_2'} = 'SIA: Mean';


T1.(1) = categorical(T1.(1))







% AMST
std_emulation = [ std_save_emul'; mean(std_save_emul)];
std_cmip6_models = [ std_save_cmip6'; mean(std_save_cmip6)];

mean_emulation = [ mean_save_emul'; mean(mean_save_emul)];
mean_cmip6_models = [ mean_save_cmip6'; mean(mean_save_cmip6)];


% Create first set of tables: std
A_LB = table(std_emulation, std_cmip6_models);
A_LB.Properties.VariableNames{'std_emulation'} = 'Emulation';
A_LB.Properties.VariableNames{'std_cmip6_models'} = 'CMIP6';


% Create first set of tables: mean
A_LB_2 = table(mean_emulation, mean_cmip6_models);
A_LB_2.Properties.VariableNames{'mean_emulation'} = 'Emulation';
A_LB_2.Properties.VariableNames{'mean_cmip6_models'} = 'CMIP6';


month_names = [ month_label' ; "Mean" ];

T1 = table(month_names, A_LB, A_LB_2);
T1.Properties.VariableNames{'month_names'} = 'Months';

T1.Properties.VariableNames{'A_LB'} = 'AMST: Standard Deviation';
T1.Properties.VariableNames{'A_LB_2'} = 'AMST: Mean';


T1.(1) = categorical(T1.(1))






%% SIA_max parameterisation comparison to CMIP6: FigA.8


close; clc
figure(12)
h = subplot(1,1,1);
set(gcf, 'Units', 'Inches', 'Position', [.5 .5 20 11])

plota_emulation = zeros(1,12);
for j = 1:12

    % Plot our emulation of SIA in 1850
    plota_emulation(j) = plot(SIA_max_param_1850(j, 1:12), 'color', colorblind(j,:), 'linewidth', 2.5); hold on

     % Plot CMIP6 SIA in 1850
    plot_cmip6 = plot(updated_hist_sia_annual_curve_all_models{j,1}(1,:), 'color', colorblind(j,:), 'linewidth', 2.5); hold on
    plot_cmip6.LineStyle = '-.';
    hhplota = plot_cmip6;
end

set(gca, 'xtick', 1:12);
tick_length = 1;
set(gca,'xticklabel', {[ blanks(tick_length) 'Jan'], [ blanks(tick_length) 'Feb'], [ blanks(tick_length) 'Mar'],...
    [ blanks(tick_length) 'Apr'],[ blanks(tick_length) 'May'],[ blanks(tick_length) 'Jun'],[ blanks(tick_length) 'Jul'],...
    [ blanks(tick_length) 'Aug'], [ blanks(tick_length) 'Sep'], [ blanks(tick_length) 'Oct'],[ blanks(tick_length) 'Noc'], ... 
[ blanks(tick_length) 'Dec'], ''});
xlabel('Month')
set(gca, 'Fontsize', 18)
xlim([1 12])
% ylim([4 25])

title('SIA_max', 'interpreter', 'none')
ylabel('Initial SIA parameteriation (year=1850) (million km^2)')
set(gca, 'Fontsize', 20)
grid 
tt2 = legend(plota_emulation, models_alone, 'location', 'southwest', 'NumColumns', 2);
tt2.Box = 'off';
set(gcf, 'PaperOrientation', 'landscape')
set(gcf,'PaperSize',[53 29.7000]);


% % Save
% temp=[' Fig2.6  ', '.pdf']; 
% saveas(gca,temp);  




%% SENSITIVITY for completion seminar: march and sept
% ONLY MARCH AND SEPT


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

close
figure(2)
set(gcf, 'Units', 'Inches', 'Position', [.4 .4 19 10])
hold on
h = subplot(1,1,1);
for month_ind = [3,9]

    
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




    % Calculate the Sensitivity: Observationally constrained SIA using CMIP6 calibrated AA
    clear coefficients
    magicc_sens_hold = [];
    magicc_sens_CMIP6_AA = [];
    for n = 1:1
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
        magicc_sens_CMIP6_AA{n} = cell2mat(magicc_sens_hold);   
    end



    % Calculate the Sensitivity: OC emulator: Observationally constrained GMST & AA
    clear coefficients
    magicc_sens_hold = [];
    magicc_sens_O1 = [];
    magicc_sens_hold_1 = [];
    for n = 1:1
        SIA = sia_mon_rearrange_final_MCMC{n}(:, month_ind);
        for j = 1:12
            for ii = 1:600
                SIA_emul = movmean(SIA{j}(ii,:), movmean_ind);
                SIA_emul = SIA_emul(index_1974_2014);

                % GT = GT_ensembles{n};
                GT = GMST_random_arrangment{n};
                GT_emul = movmean(GT(ii,:), movmean_ind); 
                GT_emul = GT_emul(index_1974_2014);

                coefficients = polyfit(GT_emul, SIA_emul, 1);
                magicc_sens_hold_1{ii} = coefficients(1);
            end  
            magicc_sens_hold{j} = (cell2mat(magicc_sens_hold_1))';
        end
        magicc_sens_O1{n} = cell2mat(magicc_sens_hold(:));   
    end


    % Calculate the CMIP6 Sensitivity
    clear coefficients
    emul_cmip6_sens = [];
    emul_cmip6_sens_hold = [];
    for n = 1
        for j = 1:12
            SIA_cmip6 = sia_emul_final_NO_BC{n}{j,month_ind};
            SIA_cmip6 = movmean(SIA_cmip6, movmean_ind);
            SIA_cmip6 = SIA_cmip6(index_1974_2014);

            GT_cmip6 = tas_global{n}{j};
            GT_cmip6 = movmean(GT_cmip6, movmean_ind);
            GT_cmip6 = GT_cmip6(index_1974_2014);

            coefficients = polyfit(GT_cmip6, SIA_cmip6, 1);
            emul_cmip6_sens_hold{j,n} = coefficients(1); 
        end
    end
    emul_cmip6_sens = cell2mat(emul_cmip6_sens_hold);




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



    mods_legend = [];
    
    if month_ind == 3
        ind = 0.25;
    elseif month_ind == 9
        ind = 0.55;
    end
    
    % Plot sensitivity as boxplots    
    obs_sens_bp = boxplot(obs_sens(:), 'positions', [ind-0.06], 'color', [0 0 0]+0.7, 'plotstyle', 'compact'); hold on
    hAx = gca;
    lines = hAx.Children;
    lines = lines(1);
    bp_cmip6 = findobj(lines, 'tag', 'Box');
    bp_cmip6.YData(1,1) = quantile(obs_sens(:), 1.0);
    bp_cmip6.YData(1,2) = quantile(obs_sens(:), 0.0);

    
    
    % Dirk's Sens Data
    if month_ind == 3
        boxplot([-1.35, -1.85], 'positions', [ind-0.03], 'color', colorss1(1,:), 'Labels', [' '], 'plotstyle','compact'); hold on 
    else
        boxplot([-3.81, -4.39], 'positions', [ind-0.03], 'color', colorss1(1,:), 'Labels', [' '], 'plotstyle','compact'); hold on 
    end

        
    boxplot(cmip6_sens, 'positions', [ind], 'color', [1 0 0], 'plotstyle','compact'); hold on    % CMIP6
    hAx = gca;
    lines = hAx.Children;
    lines = lines(1);
    bp_cmip6 = findobj(lines, 'tag', 'Box');
    bp_cmip6.YData(1,1) = quantile( cmip6_sens, 0.25);
    bp_cmip6.YData(1,2) = quantile( cmip6_sens, 0.75);


    boxplot(emul_cmip6_sens, 'positions', [ind+0.03], 'color', colorss1(2,:), 'plotstyle','compact'); hold on    % CMIP6
    hAx = gca;
    lines = hAx.Children;
    lines = lines(1);
    bp_cmip6 = findobj(lines, 'tag', 'Box');
    bp_cmip6.YData(1,1) = quantile( emul_cmip6_sens, 0.25);
    bp_cmip6.YData(1,2) = quantile( emul_cmip6_sens, 0.75);


        
    % cmip6_aa_sens = boxplot(magicc_sens_CMIP6_AA{1}(:), 'positions', [ind], 'color', [0    0.5000    0.5000], 'Labels', [' '], 'plotstyle','compact'); hold on % MAGICC using the UPPER BOUNDARY of our observationally constrained AA
    % set(cmip6_aa_sens(length(cmip6_aa_sens),:),'Visible','off')     % Remove the outliers
    % hAx = gca;
    % lines = hAx.Children;
    % lines = lines(1);
    % bp_cmip6 = findobj(lines, 'tag', 'Box');
    % bp_cmip6.YData(1,1) = quantile( magicc_sens_CMIP6_AA{1}(:), 0.25);
    % bp_cmip6.YData(1,2) = quantile( magicc_sens_CMIP6_AA{1}(:), 0.75);
    % 
    % 
    % SIE_constrained_bp = boxplot(magicc_sens_O1{1}, 'positions', [ind+0.03], 'color', colorss1(2,:), 'Labels', [' '], 'plotstyle','compact'); hold on % MAGICC using the UPPER BOUNDARY of our observationally constrained AA
    % set(SIE_constrained_bp(length(SIE_constrained_bp),:),'Visible','off')     % Remove the outliers
    % hAx = gca;
    % lines = hAx.Children;
    % lines = lines(1);
    % bp_cmip6 = findobj(lines, 'tag', 'Box');
    % bp_cmip6.YData(1,1) = quantile( magicc_sens_O1{1}, 0.25);
    % bp_cmip6.YData(1,2) = quantile( magicc_sens_O1{1}, 0.75);



    

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
    sgtitle([num2str(years(index_1974_2014(1))), '-', ...
    num2str(years(index_1974_2014(end))), ' MAGICC Sensitivity (dSIA/dGTA)'], 'fontsize', 28, 'interpreter', 'none')

    if ismember(month_ind, [1,5,9])
        ylabel('million km^2/ \circC')
    end
%     title([num2str(month_label{month_ind})])

    set(gca, 'fontsize', 22)
        
    
    length_curve_2 = 0.7;
    xlim([0 length_curve_2])    
    ylim([-6 0])    
    xticks([0.25,0.55])
    xticklabels(string({month_label{3}, month_label{9}}));
    

 
    grid
    hold on
    
end

h.Position = [0.3 0.075 0.4 0.85];
    
grid 
hLegend = legend(findall(gca,'Tag','Box'), {'CMIP6 Emulator', 'CMIP6', 'Niederdrenk and Notz, (2018)', 'Plausible Range (from Observations)'}, 'location', 'southwest', 'fontsize', 18);
hLegend.Box = 'off';

% % Save
% temp=[' Sensitivity: March and September 4 ', '.png']; 
% saveas(gca,temp);



% Save Normalised
set(gcf, 'PaperOrientation', 'landscape')
set(gcf,'PaperSize',[51 30]);

% % Save Normalised
% temp=['Fig14', '.pdf']; 
% saveas(gca,temp);















%% FigA.9: SIA VS AMST to show need for bias corrections

close all
figure(41);
set(gcf, 'Units', 'Inches', 'Position', [.5 .5 20 11])

% BOX 3
% axes('Position',[.61 .3 .38 .4])
axes('Position',[.4 .3 .5 .5])
box on

legend_h = [];
legend_obs = [];
for i = 1:12
    h = plot(movmean(mean_CMIP6_AMT_ssp{1}(:,i), 10), movmean(mean_CMIP6_SIA_ssp{1}(:,i), 10), 'linewidth', 1.5); hold on
    h.Color = [colorblind(i,:), 0.7];
    legend_h = cat(1, legend_h, h);

    % plot obs
    for m = 3:6
        plot_obs = plot(movmean(monthly_obs_temp{i}(1:170), 15), movmean(monthly_obs_SIA{m}{i}, 15), 'color', colorblind(i,:), 'linewidth', 2); hold on
        plot_obs.LineStyle = '-.';
    end
    legend_obs = cat(1, legend_obs, plot_obs);
end
ylabel('Arctic SIA (million km^{2})')
xlabel('Arctic Temperature (\circC)')
set(gca, 'Fontsize', 16)

grid
leg_labels = [month_label, 'Observations'];
tt2 = legend([legend_h; legend_obs(5)], leg_labels, 'fontsize', 12)
tt2.Box = 'off';




% set(gcf, 'PaperOrientation', 'landscape')
% set(gcf,'PaperSize',[53 29.7000]);

% % Save
% temp=[' Fig2.12 ', '.pdf']; 
% saveas(gca,temp);  





%% Fig B4| Average 2080-2100 Calibrated AMST and SIA Annual Cycle

% Plot initialisation
index = [231:251];  % 2080-2100 index
emulated_percentile_2080_2100_L = [];
cmip6_percentile_2080_2100_L = [];

cmip6_percentile_2080_2100_VL = [];
emulated_percentile_2080_2100_VL = [];

median_emul_2080_2100 = [];
median_cmip6_2080_2100 = [];
AMT_emulated_2080_2100_MEAN_hold = [];
AMT_CMIP6_2080_2100_MEAN_hold = [];

tas_store = {ssp585_tas_absolute_pi; ssp245_tas_absolute_pi; ssp126_tas_absolute_pi};

for n = 1:3
   
    AMT_CMIP6 = tas_store{n};
    if size(AMT_CMIP6,2) ~= 1
        AMT_CMIP6 = AMT_CMIP6';
    end
    AMT_CMIP6_2080_2100 = cellfun(@(v) v(231:251,:), AMT_CMIP6, 'un', 0);

    % MEAN OF 2080-2100 EMULATION
    AMT_CMIP6_2080_2100_MEAN = cellfun(@(v) mean(v), AMT_CMIP6_2080_2100, 'un', 0);
    AMT_CMIP6_2080_2100_MEAN = cell2mat(AMT_CMIP6_2080_2100_MEAN);

    AMT_CMIP6_2080_2100_MEAN_hold{n} = AMT_CMIP6_2080_2100_MEAN;

    % AMT_EMULATED = tas_cal_new_bias_corrected(:,n);
    % AMT_EMULATED = AMST_calibration(:,n);
    AMT_EMULATED = AMST_calibration_bias_corrected(:,n);
    AMT_emulated_2080_2100 = cellfun(@(v) v(231:251,:), AMT_EMULATED, 'un', 0);

    % MEAN OF 2080-2100 EMULATION
    AMT_emulated_2080_2100_MEAN = cellfun(@(v) mean(v), AMT_emulated_2080_2100, 'un', 0);
    AMT_emulated_2080_2100_MEAN = cell2mat(AMT_emulated_2080_2100_MEAN);

    AMT_emulated_2080_2100_MEAN_hold{n} = AMT_emulated_2080_2100_MEAN;
    
    % Likely range
    threshold = [17 83];
     
    % 2080-2100 mean
    emulated_percentile_2080_2100_likely_hold = prctile(cell2mat(AMT_emulated_2080_2100), threshold);
    cmip6_percentile_2080_2100_likely_hold = prctile(cell2mat(AMT_CMIP6_2080_2100), threshold);

    % Save arrays
    emulated_percentile_2080_2100_L{n} = emulated_percentile_2080_2100_likely_hold;
    cmip6_percentile_2080_2100_L{n} = cmip6_percentile_2080_2100_likely_hold;


    
    
    
    % Very likely range
    threshold = [5 95];
    
    % 2080-2100 mean 
    emulated_percentile_2080_2100_hold_VL = prctile(cell2mat(AMT_emulated_2080_2100), threshold);
    cmip6_percentile_2080_2100_hold_VL = prctile(cell2mat(AMT_CMIP6_2080_2100), threshold);
    
    % Save arrays
    emulated_percentile_2080_2100_VL{n} = emulated_percentile_2080_2100_hold_VL;
    cmip6_percentile_2080_2100_VL{n} = cmip6_percentile_2080_2100_hold_VL;
    
    
    
    % Median
    median_emul_2080_2100_hold = prctile(cell2mat(AMT_emulated_2080_2100), 50);
    median_cmip6_2080_2100_hold = prctile(cell2mat(AMT_CMIP6_2080_2100), 50);
    

    median_emul_2081{n} = median_emul_2080_2100_hold;
    median_cmip6_2080_2100{n} = median_cmip6_2080_2100_hold;
end



% Plot Initialisation
hnew = [];
hnew_L = [];
hnew_VL = [];
hnew_cmip6 = [];
hnew_cmip6_VL = [];
hnew_cmip6_L = [];
ssp_lists = {'SSP5-8.5', 'SSP2-4.5', 'SSP1-2.6'};
counter = 0;
x = 1:12;
close
figure(40)
set(gcf, 'Units', 'Inches', 'Position', [.5 .5 20 11])
for n = 1:3
    counter = counter + 1;
    h = subplot(2,3,n);


    % Very Likely Range
    percentile_emul = emulated_percentile_2080_2100_VL{n};
    percentile_cmip6 = cmip6_percentile_2080_2100_VL{n};


    % Likely Range
    percentile_L_emul = emulated_percentile_2080_2100_L{n};
    percentile_L_cmip6 = cmip6_percentile_2080_2100_L{n};


    % Very Likely Range Percentile Curves
    curve1 = percentile_emul(2,:);
    curve2 = percentile_emul(1,:);
    
    curve3 = percentile_cmip6(2,:);
    curve4 = percentile_cmip6(1,:);

    % Likely Range Percentile Curves
    curve5 = percentile_L_emul(2,:);
    curve6 = percentile_L_emul(1,:);
    
    curve7 = percentile_L_cmip6(2,:);
    curve8 = percentile_L_cmip6(1,:);


    % Create shaded Likely percentile area
    x4 = [x, fliplr(x)];
    inBetween3 = [curve5, fliplr(curve6)];
    inBetween4 = [curve7, fliplr(curve8)];

    
    
    % Fill: Calibration
    c = colmat(n,:);
    f3 = fill(x4, inBetween3, c, 'edgecolor','none');
    set(f3,'facealpha',.3)
    hold on
    
    % Fill: CMIP6
    c = [0 0 0]+0.3;
    f4 = fill(x4, inBetween4, c, 'edgecolor','none');
    set(f4,'facealpha',.3)
    hold on


    % Plot median
    median_emulation = median_emul_2081{n};
    median_cmip6 = median_cmip6_2080_2100{n};

    m_emulation = plot(x, median_emulation, 'color', colmat(n,:), 'LineWidth', 2);
    hold on   
    m_cmip6 = plot(x, median_cmip6, 'color', 'k', 'LineWidth', 2);
    
    % Plot the very likely range as line plots
    VL_range = plot(x, curve1, '--', 'color', colmat(n,:), 'LineWidth', 2);     % Very likely calibrated range
    VL_range2 = plot(x, curve2, '--', 'color', colmat(n,:), 'LineWidth', 2);

    VL_range_cmip6 = plot(x, curve3, '--', 'color', 'k', 'LineWidth', 2);       % Very likely CMIP6 range
    VL_range2_cmip6 = plot(x, curve4, '--', 'color', 'k', 'LineWidth', 2);
    
    
    hnew_VL = cat(1, hnew_VL, VL_range);
    hnew_L = cat(1, hnew_L, f3);
    hnew = cat(1, hnew, m_emulation);
    
    if n == 3
        hnew_cmip6_medium = cat(1, hnew_cmip6, m_cmip6);
        hnew_cmip6_VL = cat(1, hnew_cmip6, VL_range_cmip6);
        hnew_cmip6_L = cat(1, hnew_cmip6, f4);
    end

    % Plot labels
    grid
    xlim([1 12])
    ylim([-38 13])
    if n == 1
        ylabel(['Arctic Temperature (\circC)']) 
    end
    title(['AMT: ', ssp_lists{n}])
    set(gca, 'fontsize', 18)


    set(gca, 'xtick', 1:12);
    tick_length = 1;
    set(gca,'xticklabel', {[ blanks(tick_length) 'J'], [ blanks(tick_length) 'F'], [ blanks(tick_length) 'M'],...
        [ blanks(tick_length) 'A'],[ blanks(tick_length) 'M'],[ blanks(tick_length) 'J'],[ blanks(tick_length) 'J'],...
        [ blanks(tick_length) 'A'], [ blanks(tick_length) 'S'], [ blanks(tick_length) 'O'],[ blanks(tick_length) 'N'], ... 
    [ blanks(tick_length) 'D'], ''});
    xtickangle(0)


    posnew = get(h, 'Position');
    posnew(1) = posnew(1) - 0.08;
    posnew(2) = posnew(2) - 0.07;
    posnew(4) = posnew(4) + 0.02;
    posnew(3) = posnew(3) + 0.02;
    set(h, 'Position', posnew);
end








%__________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
% Add SIA to plot  
% Calculate the percentiles for the 2080-2100 mean: CMIP6 vs Emulation:

index = [231:251];
emulated_percentile_2080_2100_L = [];
cmip6_percentile_2080_2100_L = [];

cmip6_percentile_2080_2100_VL = [];
emulated_percentile_2080_2100_VL = [];

median_emul_2080_2100 = [];
median_cmip6_2080_2100 = [];

for n = 1:3     % Each SSP emission scenario
   
    SIA_CMIP6 = updated_hist_sia_annual_curve_all_models(:,n);
    SIA_CMIP6_2080_2100 = cellfun(@(v) v(231:251,:), SIA_CMIP6, 'un', 0);
    
    % SIA_EMULATED = sia_emulation_no_bc(:,n);
    SIA_EMULATED = calibrated_SIA_store(:,n);
    SIA_emulated_2080_2100 = cellfun(@(v) v(231:251,:), SIA_EMULATED, 'un', 0);

       
    % MEAN OF 2080-2100 EMULATION
    SIA_emulated_2080_2100_MEAN = cellfun(@(v) mean(v), SIA_emulated_2080_2100, 'un', 0);
    SIA_emulated_2080_2100_MEAN = cell2mat(SIA_emulated_2080_2100_MEAN);
    

    % MEAN OF 2080-2100 CMIP6
    SIA_CMIP6_2080_2100_MEAN = cellfun(@(v) mean(v), SIA_CMIP6_2080_2100, 'un', 0);
    SIA_CMIP6_2080_2100_MEAN = cell2mat(SIA_CMIP6_2080_2100_MEAN);

    % Likely range
    threshold = [17 83];
   
    % 2080-2100 mean
    emulated_percentile_2080_2100_likely_hold = prctile(cell2mat(SIA_emulated_2080_2100), threshold);
    cmip6_percentile_2080_2100_likely_hold = prctile(cell2mat(SIA_CMIP6_2080_2100), threshold);   
    
    % Save arrays
    emulated_percentile_2080_2100_L{n} = emulated_percentile_2080_2100_likely_hold;
    cmip6_percentile_2080_2100_L{n} = cmip6_percentile_2080_2100_likely_hold;


    
    
    
    % Very likely range
    threshold = [5 95];
    
    % 2081
    emulated_percentile_2080_2100_hold_VL = prctile(cell2mat(SIA_emulated_2080_2100), threshold);
    cmip6_percentile_2080_2100_hold_VL = prctile(cell2mat(SIA_CMIP6_2080_2100), threshold);
   
    % Save arrays
    emulated_percentile_2080_2100_VL{n} = emulated_percentile_2080_2100_hold_VL;
    cmip6_percentile_2080_2100_VL{n} = cmip6_percentile_2080_2100_hold_VL;
    
    
   
    % Median
    median_emul_2080_2100_hold = prctile(cell2mat(SIA_emulated_2080_2100), 50);
    median_cmip6_2080_2100_hold = prctile(cell2mat(SIA_CMIP6_2080_2100), 50);


    median_emul_2081{n} = median_emul_2080_2100_hold;
    median_cmip6_2080_2100{n} = median_cmip6_2080_2100_hold;
end




% Plot Initialisation
clear ff
hnew = [];
hnew_L = [];
hnew_VL = [];
hnew_cmip6 = [];
hnew_cmip6_VL = [];
hnew_cmip6_L = [];
x = 1:12;
ssp_lists = {'SSP5-8.5', 'SSP2-4.5', 'SSP1-2.6'};
for n = 1:3
    counter = counter + 1;
    h = subplot(2,3,counter);


    % Very Likely Range
    stdY_emul = emulated_percentile_2080_2100_VL{n};
    stdY_cmip6 = cmip6_percentile_2080_2100_VL{n};


   % Likely Range
    stdY3_emul = emulated_percentile_2080_2100_L{n};
    stdY3_cmip6 = cmip6_percentile_2080_2100_L{n};


    % Very Likely Range: Percentile Curves
    curve1 = stdY_emul(2,:);
    curve2 = stdY_emul(1,:);
    
    curve3 = stdY_cmip6(2,:);
    curve4 = stdY_cmip6(1,:);

    % Likely Range: Percentile Curves
    curve5 = stdY3_emul(2,:);
    curve6 = stdY3_emul(1,:);
    
    curve7 = stdY3_cmip6(2,:);
    curve8 = stdY3_cmip6(1,:);


    % Create shaded Likely percentile area
    x4 = [x, fliplr(x)];
    inBetween3 = [curve5, fliplr(curve6)];
    inBetween4 = [curve7, fliplr(curve8)];

    
    
    % Fill: Likely Calibrated Range
    c = colmat(n,:);
    f3 = fill(x4, inBetween3, c, 'edgecolor','none');
    set(f3,'facealpha',.3)
    hold on
    
    % Fill: Likely CMIP6 Range
    c = [0 0 0]+0.3;
    f4 = fill(x4, inBetween4, c, 'edgecolor','none');
    set(f4,'facealpha',.3)
    hold on


    % Plot median
    median_cal = median_emul_2081{n};
    median_cmip6 = median_cmip6_2080_2100{n};

    mm = plot(x, median_cal, 'color', colmat(n,:), 'LineWidth', 2);
    hold on
    
    mm2 = plot(x, median_cmip6, 'color', 'k', 'LineWidth', 2);
    hold on
    
     % Plot the very likely range as line plots
    VL_range = plot(x, curve1, '--', 'color', colmat(n,:), 'LineWidth', 2);         % Very likely calibrated range
    VL_range2 = plot(x, curve2, '--', 'color', colmat(n,:), 'LineWidth', 2);

    VL_range_cmip6 = plot(x, curve3, '--', 'color', 'k', 'LineWidth', 2);           % Very likely cmip6 range
    VL_range2_cmip6 = plot(x, curve4, '--', 'color', 'k', 'LineWidth', 2);
    
    
    hnew_VL = cat(1, hnew_VL, VL_range);
    hnew_L = cat(1, hnew_L, f3);
    hnew = cat(1, hnew, mm);
    
    if n == 3
        hnew_cmip6_medium = cat(1, hnew_cmip6, mm2);
        hnew_cmip6_VL = cat(1, hnew_cmip6, VL_range_cmip6);
        hnew_cmip6_L = cat(1, hnew_cmip6, f4);
    end
    
    % Plot labels
    grid
    xlim([1 12])
    ylim([0 18])  
    if n == 1
        ylabel(['SIA (million km^2)']) 
    end
    sgtitle(['Mean 2080-2100 Calibrated AMT & SIA',], 'fontsize', 22, 'interpreter', 'none')
    title(['SIA: ', ssp_lists{n}])
    set(gca, 'fontsize', 18)


    set(gca, 'xtick', 1:12);
    tick_length = 1;
    set(gca,'xticklabel', {[ blanks(tick_length) 'J'], [ blanks(tick_length) 'F'], [ blanks(tick_length) 'M'],...
        [ blanks(tick_length) 'A'],[ blanks(tick_length) 'M'],[ blanks(tick_length) 'J'],[ blanks(tick_length) 'J'],...
        [ blanks(tick_length) 'A'], [ blanks(tick_length) 'S'], [ blanks(tick_length) 'O'],[ blanks(tick_length) 'N'], ... 
    [ blanks(tick_length) 'D'], ''});
    xtickangle(0)

    posnew = get(h, 'Position');
    posnew(1) = posnew(1) - 0.08;
    posnew(2) = posnew(2) - 0.05;
    posnew(4) = posnew(4) + 0.02;
    posnew(3) = posnew(3) + 0.02;
    set(h, 'Position', posnew);
        
        
end


hnew_final = [hnew_VL(1); hnew_L(1); hnew(1); hnew_cmip6_medium(1); hnew_cmip6_VL(1); hnew_cmip6_L(1)]; 
legend_deeets = ["Emulated VLR", "Emulated LR", "Emulated Median", "CMIP6 Median", "CMIP6 VLR", "CMIP6 LR"];
[tt2] = legend(hnew_final, legend_deeets, 'Location', 'eastoutside', 'fontsize', 18, 'NumColumns', 1);
tt2.Box = 'off';

% Posistion legend in the correct place
newPosition = [0.92 0.5 0.01 0.01];
newUnits = 'normalized';
set(tt2,'Position', newPosition,'Units', newUnits);



% set(gcf, 'PaperOrientation', 'landscape')
% set(gcf,'PaperSize',[53 29.7000]);

% % Save
% temp=[' Fig2.13', '.pdf']; 
% saveas(gca,temp);  






%% FigB6: SHOW UNCONSTRAINED, CONSTRAINED, OBS:

% Constrained
AMT_C_1974_2014 = cellfun(@(v) v(125:165,:), final_final_save_MAGICC_MCMC{1}(:), 'un', 0);
AMT_C_1974_2014 = cellfun(@(v) median(v), AMT_C_1974_2014, 'un', 0);
AMT_C_1974_2014 = median(cell2mat(AMT_C_1974_2014));



% Un-constrained
AMT_CMIP6_1974_2014 = cellfun(@(v) v(125:165,:), AMT_CMIP6, 'un', 0);
AMT_CMIP6_1974_2014 = cellfun(@(v) median(v), AMT_CMIP6_1974_2014, 'un', 0);
AMT_CMIP6_1974_2014 = median(cell2mat(AMT_CMIP6_1974_2014));


% Obs
amt_obs_1974_2014 = cellfun(@(v) v(125:165,:), monthly_obs_temp, 'un', 0);
amt_obs_1974_2014 = cellfun(@(v) median(v), amt_obs_1974_2014, 'un', 0);
amt_obs_1974_2014 = cell2mat(amt_obs_1974_2014);



clc
close
figure(42)
set(gcf, 'Units', 'Inches', 'Position', [.4 .4 10 7])

plot_obs = plot(amt_obs_1974_2014, 'k', 'linewidth', 3); hold on
plot_obs.LineStyle = '--';

plot_unconstrained = plot(AMT_CMIP6_1974_2014, 'r', 'linewidth', 3);

plot_constrained = plot(AMT_C_1974_2014, 'b', 'linewidth', 3);

tt2 = legend('Observations', 'Unconstrained', 'Constrained');
tt2.Box = 'off';

set(gca, 'xtick', 1:12);
tick_length = 1;
set(gca,'xticklabel', {[ blanks(tick_length) 'Jan'], [ blanks(tick_length) 'Feb'], [ blanks(tick_length) 'March'],...
    [ blanks(tick_length) 'April'],[ blanks(tick_length) 'May'],[ blanks(tick_length) 'June'],[ blanks(tick_length) 'July'],...
    [ blanks(tick_length) 'Aug'], [ blanks(tick_length) 'Sept'], [ blanks(tick_length) 'Oct'],[ blanks(tick_length) 'Nov'], ... 
[ blanks(tick_length) 'Dec'], ''});


title(' Median 1979-2014 Seasonal Arctic Temperature Cycle ')
ylabel('Arctic Temperature (\circC)') 

grid
set(gca, 'Fontsize', 16)
xlim([1 12])



% set(gcf, 'PaperOrientation', 'landscape')
% set(gcf,'PaperSize',[53 29.7000]);

% % Save
% temp=[' Fig2.15 ', '.pdf']; 
% saveas(gca,temp);   










%% FigB5: MAGICC SIA 1979-2014 Annual Curve: OUR EMULATION COMPARISON WITH OBSERVED


SIA_MAGICC_ssps = [];
SIA_CMIP6_ssps = [];
SIA_Calibration_ssps = [];


SIA_emulated_1979_2014_L = [];
AMT_emulated_1979_2014_L = [];
SIA_emulated_1979_2014_VL = [];
AMT_emulated_1979_2014_VL = [];
SIA_median_emululated_1979_2014 = [];
AMT_median_emululated_1979_2014 = [];
SIA_emulated_CMIP6_AA_1979_2014_L = [];
AMT_emulated_CMIP6_AA_1979_2014_L = [];
SIA_emulated_CMIP6_AA_1979_2014_VL = [];
AMT_emulated_CMIP6_AA_1979_2014_VL = [];
SIA_median_emululated_CMIP6_AA_1979_2014 = [];
AMT_median_emululated_CMIP6_AA_1979_2014 = [];
        

years_index = [130:165];

for n = 1:1
    for m = 1:3

        %   MAGICC: SIA
        % SIA_MAGICC = emulated_sia_1{n}{m}(:);
        SIA_MAGICC = emulated_sia_MCMC{n}(:);
        SIA_emulated_1979_2014 = cellfun(@(v) v(years_index,:), SIA_MAGICC, 'un', 0);
        SIA_emulated_1979_2014 = cell2mat(SIA_emulated_1979_2014);



        %   MAGICC: AMT
        % AMT_MAGICC = final_final_save_MAGICC_1{n}{m}(:);
        AMT_MAGICC = final_final_save_MAGICC_MCMC{n}(:);
        AMT_emulated_1979_2014 = cellfun(@(v) v(years_index,:), AMT_MAGICC, 'un', 0);
        AMT_emulated_1979_2014 = cell2mat(AMT_emulated_1979_2014);





        % Likely range
        threshold = [17 83];

        % SIA
        SIA_emulated_1979_2014_L_hold = prctile(SIA_emulated_1979_2014, threshold);
        SIA_emulated_1979_2014_L{m,n} = SIA_emulated_1979_2014_L_hold;

        % AMT
        AMT_emulated_1979_2014_L_hold = prctile(AMT_emulated_1979_2014, threshold);
        AMT_emulated_1979_2014_L{m,n} = AMT_emulated_1979_2014_L_hold;




        % Very likely range
        threshold = [5 95];

        % SIA
        SIA_emulated_1979_2014_VL_hold = prctile(SIA_emulated_1979_2014, threshold);
        SIA_emulated_1979_2014_VL{m,n} = SIA_emulated_1979_2014_VL_hold;

        % AMT
        AMT_emulated_1979_2014_VL_hold = prctile(AMT_emulated_1979_2014, threshold);
        AMT_emulated_1979_2014_VL{m,n} = AMT_emulated_1979_2014_VL_hold;



        % Median
        % SIA
        SIA_median_emululated_1979_2014_hold = prctile(SIA_emulated_1979_2014, 50);
        SIA_median_emululated_1979_2014{m,n} = SIA_median_emululated_1979_2014_hold;

        % AMT
        AMT_median_emululated_1979_2014_hold = prctile(AMT_emulated_1979_2014, 50);
        AMT_median_emululated_1979_2014{m,n} = AMT_median_emululated_1979_2014_hold;


        
    end
    
    % MAGICC: CMIP AA: SIA
    SIA_MAGICC_CMIP6_AA = emulated_sia_CMIP6_AA{n}(:);
    SIA_emulated_CMIP6_AA_1979_2014 = cellfun(@(v) v(years_index,:), SIA_MAGICC_CMIP6_AA, 'un', 0);
    SIA_emulated_CMIP6_AA_1979_2014 = cell2mat(SIA_emulated_CMIP6_AA_1979_2014);
    
    
    % MAGICC: CMIP AA: AMT
    AMT_MAGICC_CMIP6_AA = final_final_save_MAGICC_CMIP6_AA{n}(:);
    AMT_emulated_CMIP6_AA_1979_2014 = cellfun(@(v) v(years_index,:), AMT_MAGICC_CMIP6_AA, 'un', 0);
    AMT_emulated_CMIP6_AA_1979_2014 = cell2mat(AMT_emulated_CMIP6_AA_1979_2014);
    
    
    
    % Likely range
    threshold = [17 83];
    % SIA
    SIA_emulated_CMIP6_AA_1979_2014_L_hold = prctile(SIA_emulated_CMIP6_AA_1979_2014, threshold);
    SIA_emulated_CMIP6_AA_1979_2014_L{n} = SIA_emulated_CMIP6_AA_1979_2014_L_hold;
    % AMT
    AMT_emulated_CMIP6_AA_1979_2014_L_hold = prctile(AMT_emulated_CMIP6_AA_1979_2014, threshold);
    AMT_emulated_CMIP6_AA_1979_2014_L{n} = AMT_emulated_CMIP6_AA_1979_2014_L_hold;
    
           
    % Very likely range
    threshold = [5 95];
    % SIA
    SIA_emulated_CMIP6_AA_1979_2014_VL_hold = prctile(SIA_emulated_CMIP6_AA_1979_2014, threshold);
    SIA_emulated_CMIP6_AA_1979_2014_VL{n} = SIA_emulated_CMIP6_AA_1979_2014_VL_hold;
    % AMT
    AMT_emulated_CMIP6_AA_1979_2014_VL_hold = prctile(AMT_emulated_CMIP6_AA_1979_2014, threshold);
    AMT_emulated_CMIP6_AA_1979_2014_VL{n} = AMT_emulated_CMIP6_AA_1979_2014_VL_hold;
    
    
    % Median
    % SIA
    SIA_median_emululated_CMIP6_AA_1979_2014_hold = prctile(SIA_emulated_CMIP6_AA_1979_2014, 50);
    SIA_median_emululated_CMIP6_AA_1979_2014{n} = SIA_median_emululated_CMIP6_AA_1979_2014_hold;
    % AMT
    AMT_median_emululated_CMIP6_AA_1979_2014_hold = prctile(AMT_emulated_CMIP6_AA_1979_2014, 50);
    AMT_median_emululated_CMIP6_AA_1979_2014{n} = AMT_median_emululated_CMIP6_AA_1979_2014_hold;

        
end


hnew = [];
hnew_L = [];
hnew_obs = [];
hnewmm1 = [];
hnew_obs2 = [];

 
col_mat = colorblind;
col_mat = [1 0.8 0.8; 0.8 0.8 1; 0 0.6906 0.5];
colorss2 = [1 0 0; 0 0 1; 0 0.3906 0];

SIA_observations = monthly_obs_SIA;

close

% Plot 12 panel figure
figure(40)
set(gcf, 'Units', 'Inches', 'Position', [.5 .5 20 11])
h = subplot(2,1,1);
x = 1:12;
    
for n = 1
    for m = [1]
    
        % L: AMT
        AMT_stdY_L = AMT_emulated_1979_2014_L{m,n};
        
        % L: AMT: CMIP6 AA
        AMT_stdY_L_cmip_aa = AMT_emulated_CMIP6_AA_1979_2014_L{n};

        % VL: AMT
        AMT_stdY_VL = AMT_emulated_1979_2014_VL{m,n};


        % L
        curve1 = AMT_stdY_L(2,:);
        curve2 = AMT_stdY_L(1,:);

        % VL
        curve3 = AMT_stdY_VL(2,:);
        curve4 = AMT_stdY_VL(1,:);
        
        % VL
        curve5 = AMT_stdY_L_cmip_aa(2,:);
        curve6 = AMT_stdY_L_cmip_aa(1,:);


        % L
        x_L = [x, fliplr(x)];
        inBetween_L = [curve1, fliplr(curve2)];
        inBetween_L_cmip6_aa = [curve5, fliplr(curve6)];



        % L
        c = colorss2(m,:);
        hh(1) = fill(x_L, inBetween_L, c, 'edgecolor','none');
        set(hh(1),'facealpha',.3)
        hold on


        median_y1 = AMT_median_emululated_1979_2014{m,n};
        c = col_mat(m,:);
        hh(2) = plot(x, median_y1, 'color', c, 'LineWidth', 3);
        hh(3) = plot(x, curve3, '--', 'color', c, 'LineWidth', 2);
    	plot(x, curve4, '--', 'color', c, 'LineWidth', 2);


        grid
        xlim([1 12])
        ylabel(['Arctic Monthly Temperature (\circC)']) 
        title(['Mean 1979-2014 Annual MAGICC AMT Cycle with Observations',], 'fontsize', 30, 'interpreter', 'none')
        set(gca, 'fontsize', 18)


        set(gca, 'xtick', 1:12);
        tick_length = 1;
        set(gca,'xticklabel', {[ blanks(tick_length) 'Jan'], [ blanks(tick_length) 'Feb'], [ blanks(tick_length) 'Mar'],...
            [ blanks(tick_length) 'Apr'],[ blanks(tick_length) 'May'],[ blanks(tick_length) 'Jun'],[ blanks(tick_length) 'Jul'],...
            [ blanks(tick_length) 'Aug'], [ blanks(tick_length) 'Sep'], [ blanks(tick_length) 'Oct'],[ blanks(tick_length) 'Nov'], ... 
        [ blanks(tick_length) 'Dec'], ''});



        posnew = get(h, 'Position');
        posnew(1) = posnew(1) - 0.02;
        posnew(2) = posnew(2) - 0.01;
        posnew(4) = posnew(4) + 0.02;
        posnew(3) = posnew(3) - 0.02;
        set(h, 'Position', posnew);
        
    
    end
    text(11.5, 3, '(a)', 'FontSize' , 22)
end


hold on
ylim([-33 5]); 
hold on
hh(4) = plot(1:12, median(monthly_obs_temp_2(years_index,:)), '-.', 'color', [0 0 0], 'LineWidth', 2);
AMT_OBS_PERCENTILE = prctile(monthly_obs_temp_2(years_index,:), threshold);

% L
curve1 = AMT_OBS_PERCENTILE(2,:);
curve2 = AMT_OBS_PERCENTILE(1,:);

x_L = [x, fliplr(x)];
inBetween_L = [curve1, fliplr(curve2)];
      
% L
c = [0 0 0]+ 0.6;
hh(5) = fill(x_L, inBetween_L, c, 'edgecolor','none');
set(hh(5),'facealpha',.3)
hold on
        
        
clear allChildren
allChildren = get(gca, 'Children');


legend_labels = ["Observations", "Observational Median", "SSP5-8.5 Emulated LR", "SSP5-8.5 Emulated VLR", "SSP5-8.5 Emulated Median"];
tt2 = legend(allChildren([1,2,7,4,6]), legend_labels, 'location', 'eastoutside', 'NumColumns', 1, 'fontsize', 17);
tt2.Box = 'off';      

h.Position(1) = h.Position(1)-0.03;
h.Position(3) = h.Position(3)-0.05;

% % Posistion legend in the correct place
newPosition = [0.8 0.67 0.1837 0.1591];
newUnits = 'normalized';
set(tt2,'Position', newPosition,'Units', newUnits);





% SIA
hnew = [];
hnew_L = [];
hnew_obs = [];
hnewmm1 = [];
hnew_obs2 = [];

 
col_mat = [1 0.8 0.8; 0.8 0.8 1; 0 0.6906 0.5];
colorss2 = [1 0 0; 0 0 1; 0 0.3906 0];



% Plot 12 panel figure
h = subplot(2,1,2);
x = 1:12;

    
for n = 1
    for m = [1]
    
        % L: SIA
        SIA_stdY_L = SIA_emulated_1979_2014_L{m,n};
        
        % L: SIA: CMIP6 AA
        SIA_stdY_L_cmip_aa = SIA_emulated_CMIP6_AA_1979_2014_L{n};

        % VL: SIA
        SIA_stdY_VL = SIA_emulated_1979_2014_VL{m,n};


        % L
        curve1 = SIA_stdY_L(2,:);
        curve2 = SIA_stdY_L(1,:);

        % VL
        curve3 = SIA_stdY_VL(2,:);
        curve4 = SIA_stdY_VL(1,:);
        
        % VL
        curve5 = SIA_stdY_L_cmip_aa(2,:);
        curve6 = SIA_stdY_L_cmip_aa(1,:);


        % L
        x_L = [x, fliplr(x)];
        inBetween_L = [curve1, fliplr(curve2)];
        inBetween_L_cmip6_aa = [curve5, fliplr(curve6)];



        % L
        % c = col_mat(m,:);
        c = colorss2(m,:);
        hh(1) = fill(x_L, inBetween_L, c, 'edgecolor','none');
        set(hh(1),'facealpha',.3)
        hold on
        
      


        median_y1 = SIA_median_emululated_1979_2014{m,n};
        c = col_mat(m,:);
        hh(2) = plot(x, median_y1, 'color', c, 'LineWidth', 3);
        hh(3) = plot(x, curve3, '--', 'color', c, 'LineWidth', 1);
    	plot(x, curve4, '--', 'color', c, 'LineWidth', 1);


        grid
        xlim([1 12])
        ylabel(['SIA (million km^2)']) 
        title(['Mean 1979-2014 Annual MAGICC SIA Cycle with Observations',], 'fontsize', 30, 'interpreter', 'none')
        set(gca, 'fontsize', 18)


        set(gca, 'xtick', 1:12);
        tick_length = 1;
        set(gca,'xticklabel', {[ blanks(tick_length) 'Jan'], [ blanks(tick_length) 'Feb'], [ blanks(tick_length) 'Mar'],...
            [ blanks(tick_length) 'Apr'],[ blanks(tick_length) 'May'],[ blanks(tick_length) 'Jun'],[ blanks(tick_length) 'Jul'],...
            [ blanks(tick_length) 'Aug'], [ blanks(tick_length) 'Sep'], [ blanks(tick_length) 'Oct'],[ blanks(tick_length) 'Nov'], ... 
        [ blanks(tick_length) 'Dec'], ''});
    
        posnew = get(h, 'Position');
        posnew(1) = posnew(1) - 0.02;
        posnew(2) = posnew(2) - 0.02;
        posnew(4) = posnew(4) + 0.02;
        posnew(3) = posnew(3) - 0.02;
        set(h, 'Position', posnew);



    end
    text(11.5, 16.6, '(b)', 'FontSize' , 22)
end



hold on

index = [-2:0.1:12];
index2 = linspace(-2,12, length(index));
for i = 1:length(index)
    x = linspace(index(i),index2(i)+0.5,2);
    y = [0, 1];
    plot(x, y, 'color', [0 0 0])
    ylim([0 18]); 
    hold on
end

text(1.1, 0.5, 'Ice-free', 'FontSize', 18)
 

hold on
hnew_obs = [];
for m = 1:6  
    hh(4) = plot(mean(monthly_obs_SIA_new_final{m}(years_index,:), 'omitnan'), '-.', 'color', [0 0 0], 'LineWidth', 1);
end
OBS_SIA = monthly_obs_SIA_new_final;
OBS_SIA = cellfun(@(v) v(years_index,:), OBS_SIA, 'un', 0);
OBS_SIA = cell2mat(OBS_SIA);
       
hh(5) = plot(median(OBS_SIA, 'omitnan'), 'color', [0 0 0], 'LineWidth', 2);

clear allChildren
allChildren = get(gca, 'Children');


legend_labels = ["Observational Median", "Observations", "SSP5-8.5 Emulated LR", "SSP5-8.5 Emulated Median",  "SSP5-8.5 Emulated VLR"];
tt2 = legend(allChildren([1, 2, 154, 153, 152]), legend_labels, 'location', 'eastoutside', 'NumColumns', 1, 'fontsize', 17);
tt2.Box = 'off';      

h.Position(1) = h.Position(1)-0.03;
h.Position(3) = h.Position(3)-0.05;

% % Posistion legend in the correct place
newPosition = [0.8 0.19 0.1837 0.1591];
newUnits = 'normalized';
set(tt2,'Position', newPosition,'Units', newUnits);
    



% set(gcf, 'PaperOrientation', 'landscape')
% set(gcf,'PaperSize',[53 29.7000]);
% 
% % Save Normalised
% temp=[' Fig2.14 ', '.pdf']; 
% saveas(gca,temp);  









%% 95% confidence band

for i = 1:12

    mu = mean(emulated_percentile_2080_2100_L{1});
    sigma = std(emulated_percentile_2080_2100_L{1}); 
    alpha = 0.05; % Significance level
    n = 1;
    z = norminv(1 - alpha/2); % Z-score for the normal distribution
    lower_bound = mu - z * sigma / sqrt(n);
    upper_bound = mu + z * sigma / sqrt(n);
    
    range_of_interest = [cmip6_percentile_2080_2100_L{1}(1,i), cmip6_percentile_2080_2100_L{1}(2,i)];
   
    
    if range_of_interest(1) >= lower_bound(i) && range_of_interest(2) <= upper_bound(i)
        disp('The range falls within the 95% confidence band.');
    else
        disp('The range does not fall within the 95% confidence band.');
    end
    
end







%%  Fig A6: Plot SIA and AMST bias correct and summer underestimation timeseries them together with observations


% Find the mean SIA
mean_CMIP6_SIA_ssp = [];
for n = 1:3
    mean_CMIP6_SIA_month = [];
    for i = 1:12
        mean_CMIP6_SIA = [];
        for j = 1:12
            mean_CMIP6_SIA = cat(1, mean_CMIP6_SIA, updated_hist_sia_annual_curve_all_models{j,n}(:,i)');
        end
        mean_CMIP6_SIA_month = cat(2, mean_CMIP6_SIA_month, mean(mean_CMIP6_SIA)');
    end
    mean_CMIP6_SIA_ssp{n} = mean_CMIP6_SIA_month;
end


% Mean AMT
mean_CMIP6_AMT_ssp = [];
for n = 1:3
    mean_CMIP6_AMT_month = [];
    for i = 1:12
        mean_CMIP6_AMT = [];
        for j = 1:12
            mean_CMIP6_AMT = cat(1, mean_CMIP6_AMT, tas_store{n}{j}(:,i)');
        end
        mean_CMIP6_AMT_month = cat(2, mean_CMIP6_AMT_month, mean(mean_CMIP6_AMT)');
    end
    mean_CMIP6_AMT_ssp{n} = mean_CMIP6_AMT_month;
end


close all
figure43 = figure(43);
set(gcf, 'Units', 'Inches', 'Position', [.5 .5 20 11])


% RE-ARRANGE
tas_cmip6_rearrange_hold = [];
sia_cmip6_rearrange_hold = [];
tas_cmip6_rearrange = [];
sia_cmip6_rearrange = [];
for n = 1:3  
    for j = 1:12  
        for i = 1:12
            tas_cmip6_rearrange_hold{j,i} = tas_store{n}{j}(:,i)';
            sia_cmip6_rearrange_hold{j,i} = updated_hist_sia_annual_curve_all_models{j,n}(:,i)';
        end
    end
    tas_cmip6_rearrange{n} = tas_cmip6_rearrange_hold;
    sia_cmip6_rearrange{n} = sia_cmip6_rearrange_hold;
end


% TAS
for i = 1:12
    h = subplot(4,3,i, 'Parent', figure43);


    % CMIP6 Median 
    y1_cmip6 = mean(cell2mat(tas_cmip6_rearrange{1}(:,i)), 'omitnan');

    x = 1850:2100;
    x_cmip6 = 1850:2100;


    % Percentile: CMIP6
    stdY_cmip6 = prctile(cell2mat(tas_cmip6_rearrange{1}(:,i)), [17 83]);


    % Create curves: CMIP6 
    curve3 = stdY_cmip6(2,:);
    curve4 = stdY_cmip6(1,:);


    % Shading std area: CMIP6 
    x3 = [x_cmip6, fliplr(x_cmip6)];
    inBetween2 = [curve3, fliplr(curve4)];


    % Fill: CMIP6
    d = [0 0 0]+0.1;
    f2 = fill(x3, inBetween2, d, 'edgecolor','none');
    set(f2,'facealpha',.3)
    hold on


%   Plot emul and cmip6 means
    plt_cmip6 = plot(x_cmip6, y1_cmip6, 'color', [0 0 0], 'LineWidth', 1.5);
    x_obs = 1850:2020;
    x_to_x_obs = [130: 170];
    plt_obs = plot(x_obs(x_to_x_obs), movmean(monthly_obs_temp{i}(x_to_x_obs), 1), 'color', [1 0 0], 'LineWidth', 1);


    title([month_label(i)]);
    set(gca,'FontSize', 14)

    
    
    pos = get(h, 'Position');
    posnew = pos;
    posnew(2) = posnew(2) - 0.02;
    posnew(3) = posnew(3) - 0.09;

    if ismember(i, [1,4,7,10])
        posnew(1) = posnew(1) - 0.08;
    elseif ismember(i, [2,5,8,11])
        posnew(1) = posnew(1) - 0.21;
    elseif ismember(i, [3,6,9,12])
        posnew(1) = posnew(1) - 0.345;
    end
    set(h, 'Position', posnew);

%   Y Label
    if ismember(i, [4])
        ylab = ylabel('Arctic Absolute Temperature (\circC)');
    end

%   X Label
    if ismember(i, [11])
        xlab = xlabel('Year');
    end    

    xlim([1850 2100])
    grid

end
ylab.Position(2) = ylab.Position(2)-15;




% SIA
counter = [3,3,3,6,6,6,9,9,9,12,12,12];
for i = 1:12
    h = subplot(4,6,i+counter(i), 'Parent', figure43);


    % CMIP6 Median 
    y1_cmip6 = mean(cell2mat(sia_cmip6_rearrange{1}(:,i)), 'omitnan');

    x = 1850:2100;
    x_cmip6 = 1850:2100;


    % Percentile: CMIP6
    stdY_cmip6 = prctile(cell2mat(sia_cmip6_rearrange{1}(:,i)), [17 83]);


    % Create curves: CMIP6 
    curve3 = stdY_cmip6(2,:);
    curve4 = stdY_cmip6(1,:);


    % Shading std area: CMIP6 
    x3 = [x_cmip6, fliplr(x_cmip6)];
    inBetween2 = [curve3, fliplr(curve4)];


    % Fill: CMIP6
    d = [0 0 0]+0.1;
    f2 = fill(x3, inBetween2, d, 'edgecolor','none');
    set(f2,'facealpha',.3)
    hold on


%   Plot emul and cmip6 means
    plt_cmip6 = plot(x_cmip6, y1_cmip6, 'color', [0 0 0], 'LineWidth', 1.5);
    x_obs = 1850:2019;
    x_to_x_obs = [130: 170];
    for m = 3:6
        plt_obs = plot(x_obs(x_to_x_obs), movmean(monthly_obs_SIA{m}{i}(x_to_x_obs), 1), 'color', [0 0 1], 'LineWidth', 1); hold on
    end

   title([month_label(i)]);
    set(gca,'FontSize', 14)
    
    
    pos = get(h, 'Position');
    posnew = pos;
    posnew(2) = posnew(2) - 0.02;
    posnew(3) = posnew(3) + 0.02;

    if ismember(i, [1,4,7,10])
        posnew(1) = posnew(1) - 0.01;
    elseif ismember(i, [2,5,8,11])
    elseif ismember(i, [3,6,9,12])
        posnew(1) = posnew(1) + 0.018;
    end
    set(h, 'Position', posnew);

%   Y Label
    if ismember(i, [4])
        ylab = ylabel('Arctic SIA (million km^{2})');
    end

%   X Label
    if ismember(i, [11])
        xlab = xlabel('Year');
    end    

    xlim([1850 2100])
    grid

end
ylab.Position(2) = ylab.Position(2)-10;
set(gcf, 'PaperOrientation', 'landscape')
set(gcf,'PaperSize',[53 29.7000]);

% % Save
% temp=[' Fig2.11 ', '.pdf']; 
% saveas(gca,temp); 





