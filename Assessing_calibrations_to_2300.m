%% Assessing Calibrations to 2300: Section 3.1


%% Part 1: Extract temperature using ARO IPCC WG1 AR6 shape file
% Open shape file to extract ARO Temp Data

shapefilePath = '/Users/smchilcott/Documents/PhD UniMelb/THESIS/Chapte 2_Winter sea ice to 2300/DATA_2300/IPCC-WGI-reference-regions-v4_shapefile/IPCC-WGI-reference-regions-v4.shp';

% Read the shapefile
shapefileData = shaperead(shapefilePath)

Continent_hold = [];
for index = 1:size(shapefileData,1)
    % disp(shapefileData(index).Acronym);
    Continent_hold{index} = shapefileData(index).Acronym;     
end

ARO_ind = find(Continent_hold == "ARO")

% Double check it is the correct ARO
disp(Continent_hold(ARO_ind))


aroIndex = find(strcmpi({shapefileData.Acronym}, 'ARO'));
aroCoordinates_lon = shapefileData(aroIndex).X;
aroCoordinates_lat = shapefileData(aroIndex).Y;

figure(2); scatter(aroCoordinates_lon, aroCoordinates_lat)



%% ________ Extract the CMIP6 variable: SIAREAN ________ %

% Initialisation
rootdir = cd;
filelist = dir(fullfile(rootdir, '**/*.nc'));

full_file_path = [];
name_filelist = [];
for k = 1:length(filelist)
    name_filelist{k} = filelist(k).name;
    siarean_folder{k} = filelist(k).folder;
    full_file_path{k} = cat(2, siarean_folder{k}, '/', name_filelist{k});
end
siarean = full_file_path'
name_filelist = name_filelist';





% Get all the model names together
model = [];
for i = 1 : length(name_filelist)  

    info = ncinfo(siarean{i});
    variables = info.Variables;

    filename = strsplit(char( name_filelist{i}), '_');
    model{i} = filename{3};

end
model = string(model);
[~, ind_models] = unique(model)
siarean_2300_model_names = model(ind_models);       % save the names of all the models with data available to 2300




%% Open SIAREAN netcdf files

siarean_store = [];
for i = 1 : length(siarean)  

    info = ncinfo(siarean{i});
    variables = info.Variables;

    filename = strsplit(char( name_filelist{i}), '_');  
    variable = filename{1};
    ensemble = filename{5};
    model_new = filename{3};
    dates = filename{end};


    names = [];
    for index = 1:length(variables)                     
        names{index} = variables(index).Name;     
    end
    names = string(names);


    ind_1 = find(names == "siarean");   
    var_name = "siarean";

    % Open files
    file_open = netcdf.open(siarean{i}); 
    varID = netcdf.inqVarID(file_open, var_name);    % Open SIA data
    siarean_store{i} = double(netcdf.getVar(file_open,varID));
end





% _____________ Splice time periods (historical and SSP data etc) together _____________ %
ind_models_new = [ind_models; length(siarean)+1];
siarean_store_combined = [];
for i = 1:length(ind_models_new)-1

    siarean_store_combined_hold = siarean_store(ind_models_new(i):ind_models_new(i+1)-1);

    disp(model(ind_models_new(i):ind_models_new(i+1)-1))

    siarean_store_combined_hold_2 = [];
    for j = 1:length(siarean_store_combined_hold)
        siarean_store_combined_hold_2 = cat(2, siarean_store_combined_hold_2, siarean_store_combined_hold{j}(:)');
    end

    siarean_store_combined{i} = siarean_store_combined_hold_2;
end
cellfun(@size, siarean_store_combined, 'UniformOutput', false)      % Check the sizes of all models are correct




%% RESHAPE SIAREAN INTO EACH MONTH


% Create Labels
month_label = {'January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December'};
siarean_final = [];
for i = 1:length(siarean_store_combined)

    sia_monthly = siarean_store_combined{i};
   
    month_idx_year = 1:12;
    len_sia_vector = length(sia_monthly);
    num_years_simulation = len_sia_vector/12;
    month_idx = repmat(month_idx_year, 1, num_years_simulation);
    
    monthly_sia_store = [];
    for nn = 1:12
        idx = find(month_idx == nn)';
        monthly_sia = sia_monthly(idx)';

        monthly_sia_store = cat(2, monthly_sia_store, monthly_sia);
    
    end
    siarean_final{i} = monthly_sia_store;
end
cellfun(@size, siarean_final, 'UniformOutput', false)




%% PLOT SIAREAN IN MARCH TO 2300 TO TEST 

close all
for j = 1:length(siarean_final)
    for i = 3
        plot(siarean_final{j}(:,3)); hold on
    end
end










%% _____________________ Get tas (TEMPERATURE) to 2300 _____________________%

GMST_path = [];
full_file_path = [];
mag_folder = [];
name_filelist = [];

rootdir = cd;
filelist = dir(fullfile(rootdir, '**/*.nc'));

for k = 1:length(filelist)
    name_filelist{k} = filelist(k).name;
    mag_folder{k} = filelist(k).folder;
    full_file_path{k} = cat(2, mag_folder{k}, '/', name_filelist{k});
end
tas_path = full_file_path'





% Get all model name data together for temperature data 
model = [];
for i = 1 : length(name_filelist)  

    info = ncinfo(tas_path{i});
    variables = info.Variables;

    filename = strsplit(char( name_filelist{i}), '_');
    model{i} = filename{3};

end

model = string(model);
[~, ind_models] = unique(model);
tas_2300_model_names = model(ind_models); 








%%  Open tas netcdf files

tas_store_2300 = [];
lat_hold = [];
lon_hold = [];
for i = 1 : length(tas_path)  

    info = ncinfo(tas_path{i});
    variables = info.Variables;

    filename = strsplit(char( name_filelist{i}), '_');  
    variable = filename{1};
    ensemble = filename{5};
    model_new = filename{3};
    dates = filename{end};


    names = [];
    for index = 1:length(variables)    
        % disp(variables(index).Name);
        names{index} = variables(index).Name;     
    end
    names = string(names);


    ind_1 = find(names == "tas");   
    var_name = "tas";

    % Open files
    file_open = netcdf.open(tas_path{i}); 
    varID = netcdf.inqVarID(file_open, var_name);    % Open tas data
    tas_store_2300{i} = double(netcdf.getVar(file_open,varID));


    % Extract lat and lon data
    file_open = netcdf.open(tas_path{i}); 
    varID_lat = netcdf.inqVarID(file_open, "lat");    % Open tas data
    varID_lon = netcdf.inqVarID(file_open, "lon");    % Open tas data
    lat_hold{i} = double(netcdf.getVar(file_open,varID_lat));
    lon_hold{i} = double(netcdf.getVar(file_open,varID_lon));

end









%% ________________ EXTRACT GLOBAL-MEAN ________________ %%

clc
Kelvin = -273.15;
tas_global_2300 = [];
tas_global_2300_final = [];
test_tas_weighted = [];
for i = 1:length(tas_store_2300)

    siza = size(tas_store_2300{i});

    lat = lat_hold{i};
    lon = lon_hold{i};


    % add lat weighting
    weights = cos(deg2rad(lat_hold{i})); 
    weights = weights';
    weights = repmat(weights, [length(lon_hold{i}), 1]);



    tas_global_2300_hold = [];
    regriddedTempData_store = [];
    for j = 1:siza(3)

        tempData = tas_store_2300{i}(:,:,j);
        tempData = sum(sum(tempData.*weights,1),2) ./ sum(sum(weights,1),2);
        regriddedTempData_store{j} = tempData + Kelvin;

    end
   
    tas_global_2300_final{i} = cell2mat(regriddedTempData_store);
end

cellfun(@size, tas_global_2300_final, 'UniformOutput', false)







%% Splice GMST timeperiods together (historical/ SSPs etc, as some of the SSP data is separated into multiple cells)

clc
ind_models_new = [ind_models; length(tas_global_2300_final)+1];
tas_store_combined_2300 = [];
for i = 1:length(ind_models_new)-1

    tas_store_combined_2300_hold = tas_global_2300_final(ind_models_new(i):ind_models_new(i+1)-1);

    disp(model(ind_models_new(i):ind_models_new(i+1)-1))

    tas_store_combined_2300_hold_2 = [];
    for j = 1:length(tas_store_combined_2300_hold)
        tas_store_combined_2300_hold_2 = cat(2, tas_store_combined_2300_hold_2, tas_store_combined_2300_hold{j}(:)');
    end

    tas_store_combined_2300{i} = tas_store_combined_2300_hold_2; 
end

cellfun(@size, tas_store_combined_2300, 'UniformOutput', false)

tas_2300_model_names = model(ind_models);           % save the temperature model names






%% RESHAPE THE GMST INTO EACH MONTH


% Create Labels
month_label = {'January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December'};
tas_global_2300 = [];

for i = 1:length(tas_store_combined_2300)

    tas_monthly = tas_store_combined_2300{i};
   
    month_idx_year = 1:12;
    len_sia_vector = length(tas_monthly);
    num_years_simulation = len_sia_vector/12;
    month_idx = repmat(month_idx_year, 1, num_years_simulation);
    
    monthly_tas_store = [];
    for nn = 1:12
        idx = find(month_idx == nn)';
        monthly_tas = tas_monthly(idx)';

        monthly_tas_store = cat(2, monthly_tas_store, monthly_tas);
    
    end
    tas_global_2300{i} = monthly_tas_store;
end


cellfun(@size, tas_global_2300, 'UniformOutput', false)


tas_global_2300 = cellfun(@(v) mean(v,2), tas_global_2300, 'UniformOutput', false);





%% Convert tas global to tas global anomaly

tas_global_2300_anom = [];
for j = 1:length(tas_global_2300)
    pi = mean(tas_global_2300{j}(1:50));
    tas_global_2300_anom{j} = tas_global_2300{j} - pi;
end




%% SCATTER TAS GLOBAL ANOMALY WITH TAS.MAG FILES TO CHECK THEY ARE THE SAME

tas_store_extraced = [ tas_global{1}(1), tas_global{1}(2), tas_global{1}(3), tas_global{1}(6), {NaN}, tas_global{1}(12) ];      % Extract the models from the 2100 calibration that match the 2300 models (IPSL was not used in the 2100 calibration, hence there is a NaN here)
close all
for j = 1
    plot(tas_global_2300_anom{j}, 'r'); hold on
    plot(tas_store_extraced{j}, 'k')
end



%% PLOT GMST TO 2300 TO TEST 

close all
for j = 1:length(tas_global_2300)
    for i = 3
        plot(tas_global_2300{j}); hold on
    end
end




%% Scatter 2300 GMST vs SIAREAN (March)

tas_and_siarean_2300_model_names = tas_2300_model_names([1,2,3,5,6]);
tas_global_2300_same_models_siarean = tas_global_2300([1,2,3,5,6]);
tas_global_2300_anom_same_models_siarean = tas_global_2300_anom([1,2,3,5,6]);

clc
close all
for j = 1:length(tas_global_2300_same_models_siarean)
    for i = 3
        scatter( tas_global_2300_same_models_siarean{j}, siarean_final{j}(:,3) ); hold on
    end
end









%% ________________ EXTRACT ARCTIC TAS AND WEIGHT ________________ %%


% 65-90 degrees north = 25 degrees
% assume lat is in the second index

clc
Kelvin = -273.15;
tas_amst_2300_final_just_weighted = [];
for i = 1:length(tas_store_2300)
    siza = size(tas_store_2300{i});

    % Define the area of the Arctic (as used from the ARO data)
    arctic_lat = find(lat_hold{i} >= 67.5); 


    % Calculate weights for lat weighting
    weights = cos(deg2rad(lat_hold{i})); 
    weights = weights';
    weights = repmat(weights, [length(lon_hold{i}), 1]);  

    tas_amst_2300_hold = []; 
    for j = 1:siza(3)
        

        % Weight the temperature based on its latitude
        tempData = tas_store_2300{i}(:,arctic_lat,j);
        tempData = sum(sum(tempData.*weights,1),2) ./ sum(sum(weights,1),2);
        tas_amst_2300_hold{j} = tempData + Kelvin;

    end
   
    tas_amst_2300_final_just_weighted{i} = cell2mat(tas_amst_2300_hold);
end

cellfun(@size, tas_amst_2300_final_just_weighted, 'UniformOutput', false)





%% EXTRACT ARCTIC TAS, WEIGHT AND EXTRACT ARO COORDINATES


ind = find(isnan(aroCoordinates_lat) == 1);
aroCoordinates_lat(ind) = [];
aroCoordinates_lon(ind) = [];
tas_amst_2300_final = [];

for i = 1:length(tas_amst_2300_final_just_weighted)

    siza = size(tas_amst_2300_final_just_weighted{i});
    arctic_lat = find(lat_hold{i} >= 67.5);

    % Get lat and lon data from each model
    lon = lon_hold{i}; 
    lon = lon - 180;
    lat = lat_hold{i};


    % Find the lat and lon indices closest to the ARO on my grid
    [~, latIndex] = min(abs(lat - aroCoordinates_lat));
    [~, lonIndex] = min(abs(lon - aroCoordinates_lon));

    % Extract the closest indices from my temp data grid
    lat_hold_new = lat(latIndex);
    lon_hold_new = lon(lonIndex);

    lat = flip(lat_hold{i});

    % Create a grid of the original lat and lon data
    [grid_lon_mesh, grid_lat_mesh] = meshgrid(lon, lat);

    % Extract the data from the shape of the ARO grid
    in_shape = inpolygon(grid_lat_mesh(:), grid_lon_mesh(:), lat_hold_new, lon_hold_new);
    in_shape_data_extraction = in_shape;


    for j = 1:siza(3)
        
        % Extract ARO region
        data = tas_amst_2300_final_just_weighted{i}(:,:,j);
        data = rot90(data);

        data_within_shape = data(in_shape_data_extraction); 
    end
   
    tas_amst_2300_final{i} = cell2mat(data_within_shape);
end

cellfun(@size, tas_amst_2300_final, 'UniformOutput', false)







%% Splice AMST timeperiods together

clc
ind_models_new = [ind_models; length(tas_amst_2300_final)+1];
tas_amst_store_combined_2300 = [];
for i = 1:length(ind_models_new)-1

    tas_store_combined_2300_hold = tas_amst_2300_final(ind_models_new(i):ind_models_new(i+1)-1);

    disp(model(ind_models_new(i):ind_models_new(i+1)-1))

    tas_store_combined_2300_hold_2 = [];
    for j = 1:length(tas_store_combined_2300_hold)
        tas_store_combined_2300_hold_2 = cat(2, tas_store_combined_2300_hold_2, tas_store_combined_2300_hold{j}(:)');
    end

    tas_amst_store_combined_2300{i} = tas_store_combined_2300_hold_2;

end

cellfun(@size, tas_amst_store_combined_2300, 'UniformOutput', false)







%% RESHAPE THE AMST INTO EACH MONTH

% Create Labels
month_label = {'January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December'};
tas_amst_2300 = [];

for i = 1:length(tas_amst_store_combined_2300)

    tas_monthly = tas_amst_store_combined_2300{i};
   
    month_idx_year = 1:12;
    len_sia_vector = length(tas_monthly);
    num_years_simulation = len_sia_vector/12;
    month_idx = repmat(month_idx_year, 1, num_years_simulation);
    
    monthly_tas_store = [];
    for nn = 1:12
        idx = find(month_idx == nn)';
        monthly_tas = tas_monthly(idx)';

        monthly_tas_store = cat(2, monthly_tas_store, monthly_tas);
    
    end
    tas_amst_2300{i} = monthly_tas_store;
end

cellfun(@size, tas_amst_2300, 'UniformOutput', false)





%% Arctic Annual temperature: 2300


aat_2300 = [];
aat_2300_anom = [];
for j = 1:length(tas_amst_2300)
    aat_2300{j} = mean(tas_amst_2300{j},2);

    pi = mean(aat_2300{j}(1:50));
    aat_2300_anom{j} = aat_2300{j} - pi;
end

cellfun(@size, aat_2300, 'UniformOutput', false)




%% PLOT AMST MARCH TO 2300 TO TEST

tas_store_extracted = [tas_store{1}(1:3), tas_store{1}(6), tas_store{1}(12)];
tas_amst_2300_extracted2 = [tas_amst_2300(1:4), tas_amst_2300(6)];
amst_2300_models_extraced = [amst_2300_models(1:4), amst_2300_models(6)];

close all; clc
figure(40)
set(gcf, 'Units', 'Inches', 'Position', [.4 .4 19 10])

for j = 1:length(tas_amst_2300)
    for i = 1:12
        subplot(3,4,i)
        plot(tas_amst_2300_extracted2{j}(:,i), 'r', 'linewidth', 1); hold on
        plot(tas_store_extracted{j}(:,i), 'k', 'linewidth', 2)
    end
    sgtitle(amst_2300_models_extraced(j))
end








%% INITIALISATION: TESTING 2100 CALIBRATIONS ON 2300 DATA

colmat = [1 0 0; 0 0 1; 0.9290 0.6940 0.1250];

SIA_2300_models = ["ACCESS-CM2"    "ACCESS-ESM1-5"    "CESM2-WACCM"    "IPSL-CM6A-LR"    "MRI-ESM2-0"];
amst_2300_models = [ "ACCESS-CM2"    "ACCESS-ESM1-5"    "CESM2-WACCM"    "CanESM5"    "IPSL-CM6A-LR"    "MRI-ESM2-0"];

tas_amst_combined_2300_ruby = tas_amst_2300;
aat_2300_anom_ruby = aat_2300_anom;
aat_2300_absolute_ruby = aat_2300;
tas_global_2300_anom_ruby = tas_global_2300_anom;

aat_2300_anom_ruby = cellfun(@(v) v', aat_2300_anom_ruby, 'UniformOutput', false);
aat_2300_absolute_ruby  = cellfun(@(v) v', aat_2300_absolute_ruby, 'UniformOutput', false);







%% Assign 

% CMIP6 2100
amst_2100 = [tas_store{1}(1:3),    tas_store{1}(12)];
updated_hist_sia_annual_curve_all_models_extracted =  [updated_hist_sia_annual_curve_all_models(1:3, 1);    updated_hist_sia_annual_curve_all_models(12, 1)];

% 2300 data
tas_amst_2300_nocanesm = [tas_amst_combined_2300_ruby(1:3), tas_amst_combined_2300_ruby(6)];
if size(tas_amst_2300_nocanesm{3})<451
    tas_amst_2300_nocanesm{3} = [tas_amst_2300_nocanesm{3}; repelem(NaN, 12)];
end
if size(siarean_final{3})<451
    siarean_final{3} = [siarean_final{3}; repelem(NaN, 12)];
end
siarean_final_extraced = [siarean_final(1:3), siarean_final(5)];



%% Calibrate IPSL to 2100 (as we didn't use this model in the original calibration)

model_ind = [1, 2, 3, 6, 12]; % Need to calibrate to IPSL and need to remove CanESM5 (4) when applying to SIA

AAT_IPSL_anomaly_2300 = [];
AAT_IPSL_absolute_2300 = [];
if length(tas_global_2300_anom_ruby{3}) < 451
    tas_global_2300_anom_ruby{3} = [tas_global_2300_anom_ruby{3}; NaN];
else
    tas_global_2300_anom_ruby{3} = tas_global_2300_anom_ruby{3}(1:451); 
end
if length(aat_2300_anom_ruby{3}) < 451
    aat_2300_anom_ruby{3} = [aat_2300_anom_ruby{3}, NaN];
else
    aat_2300_anom_ruby{3} = aat_2300_anom_ruby{3}(1:451); 
end

tas_global_extracted = [tas_global{1}(1:3); tas_global{1}(6); {NaN}; tas_global{1}(12)];
tas_arctic_extracted = [tas_arctic_annual{1}(1:3); tas_arctic_annual{1}(6); {NaN}; tas_arctic_annual{1}(12)];


close all; clc
for j = 1:6                                      % j = CMIP6 models used in this analysis                                        


    rw = tas_global_2300_anom_ruby{j}(1:251)';                          % CMIP6 global temperature anomaly (1850:2100) - GMST
    raa = aat_2300_anom_ruby{j}(1:251);                  % CMIP6 Arctic annual surface air temperature anomaly (1850:2100) - AAST


    RC = polyfit(rw, raa, 1); 

    rc_save_2300{j} = RC(1);


    scatter( tas_global_2300_anom_ruby{j}, aat_2300_anom_ruby{j}, 15, 'b', 'filled' ); hold on

    if ismember(j, [1:4,6])
        scatter( tas_global_extracted{j}, tas_arctic_extracted{j}, 15, 'k', 'filled' ); hold on
    end

end






%% _______________ EMULATOR FORCED WITH 2300 GMST _______________ %% 

%% Step i: GMST - AAST
% Calculate the Arctic Amplification to pass through emulator using 2300 GMST with 2100 calibrated values

rc_save_2300_extraced = [rc_save(model_ind(1:4),1); rc_save_2300(5); rc_save(12,1)];
rc_save_2300_extraced = rc_save_2300;

AAT_emulation_anomaly_2300 = [];
AAT_emulation_absolute_2300 = [];

for j = 1:length(aat_2300_anom_ruby)                                    % j = CMIP6 models used in this analysis                                        
    
    
    rw = tas_global_2300_anom_ruby{j};                          % CMIP6 global temperature anomaly (1850:2100) - GMST
    raa = aat_2300_anom_ruby{j};                  % CMIP6 Arctic annual surface air temperature anomaly (1850:2100) - AAST
    
    
    % Calculate the Arctic annual temperature anomaly         
    AAT_emulation_anomaly_2300{j} = rw .* rc_save_2300_extraced{j}; 
    
    
    
    % Calculate the absolute Arctic annual temperature              
    PI_cmip6 = mean(aat_2300_absolute_ruby{j}(1:50));               % Extract the CMIP6 absolute AAST pre-indsutrial mean
    AAT_emulation_absolute_2300{j} = (rw .* rc_save_2300_extraced{j}) + PI_cmip6;         % Use the CMIP6 absolute AAST pre-indsutrial mean to emulate the absolute AAST

end
% ABSOLUTE 
if size(AAT_emulation_absolute_2300{1}, 1) ~= 1
    AAT_emulation_absolute_2300 = cellfun(@transpose, AAT_emulation_absolute_2300, 'UniformOutput', false);
end
if size(AAT_emulation_absolute_2300, 2) ~= 1
    AAT_emulation_absolute_2300 = AAT_emulation_absolute_2300';
end
if length(AAT_emulation_absolute_2300{3}) ~= 451
    AAT_emulation_absolute_2300{3} = [AAT_emulation_absolute_2300{3}, NaN];
end


% ANOMALY
if size(AAT_emulation_anomaly_2300{1}, 1) ~= 1
    AAT_emulation_anomaly_2300 = cellfun(@transpose, AAT_emulation_anomaly_2300, 'UniformOutput', false);
end
if size(AAT_emulation_anomaly_2300, 2) ~= 1
    AAT_emulation_anomaly_2300 = AAT_emulation_anomaly_2300';
end
if length(AAT_emulation_anomaly_2300{3}) ~= 451
    AAT_emulation_anomaly_2300{3} = [AAT_emulation_anomaly_2300{3}, NaN];
end



% CMIP6 ABSOLUTE 
if size(aat_2300_absolute_ruby{1}, 1) ~= 1
    aat_2300_absolute_ruby = cellfun(@transpose, aat_2300_absolute_ruby, 'UniformOutput', false);
end
if size(aat_2300_absolute_ruby, 2) ~= 1
    aat_2300_absolute_ruby = aat_2300_absolute_ruby';
end
if length(aat_2300_absolute_ruby{3}) ~= 451
    aat_2300_absolute_ruby{3} = [aat_2300_absolute_ruby{3}, NaN];
end


% CMIP6 ANOM
if size(aat_2300_anom_ruby{1}, 1) ~= 1
    aat_2300_anom_ruby = cellfun(@transpose, aat_2300_anom_ruby, 'UniformOutput', false);
end
if size(aat_2300_anom_ruby, 2) ~= 1
    aat_2300_anom_ruby = aat_2300_anom_ruby';
end
if length(aat_2300_anom_ruby{3}) ~= 451
    aat_2300_anom_ruby{3} = [aat_2300_anom_ruby{3}, NaN];
end
clc

cellfun(@size, aat_2300_absolute_ruby, 'UniformOutput', false)
cellfun(@size, aat_2300_anom_ruby, 'UniformOutput', false)
cellfun(@size, AAT_emulation_anomaly_2300, 'UniformOutput', false)
cellfun(@size, AAT_emulation_absolute_2300, 'UniformOutput', false)










%% Stage ii: AAT - AMST: CALIBRATION


% Initialise
year_ind = [1:251];
clear pi
x = linspace(0,2*pi,13);
year = 1850:2100;      
parameters_store = []; 
calibrated_monthlytemp_save = []; 
counter = 0;   


% Get AMST data all ssps as an array
tas_store_hold = cellfun(@(v) v(1:251,:), tas_amst_combined_2300_ruby, 'UniformOutput', false);


tic
% Calibration routine for IPSL as this model was not in the original 2100 calibration routine so calibrate this model to 2100
for j = [1:6]


%   Counter for plotting
    counter = counter + 1;


    tas = tas_store_hold{j};     % extract the CMIP6 AMST of SSP5-8.5
    add_col = tas(2:end,1);             % make the annual temperature curve Jan-Jan of following year
    add_col_final = tas_amst_combined_2300_ruby{j}(252,1);
    add_col = [add_col; add_col_final];
    tas2 = [tas, add_col];
    y = tas2;


    AAMST = AAT_emulation_absolute_2300{j}(1:251)';                  % Get the emulated AAST from all SSP scenarios used to force the AMST parameterisation
    AAMST = repmat(AAMST, [1 13]);                          % Make it the same size as the CMIP6 annual temperature curve we are calibrating to

    e = 0.3;

    ffunc = @(ppp)AMST_parameterisation_publication2(ppp,x,y,AAMST,e);    % Calibration function (in another script (AMST_parameterisation_publication.m))


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


    % Cap the linear regressions of the calibration parameters so they don't change after the AAST rises above 6 degrees celsius
    ind_temp = 6;
    indf = AAMST > ind_temp;
    indf = find(indf==1);
    f(indf) = (f1.*ind_temp) + f2;
    g(indf) = (g1.*ind_temp) + g2;
    a(indf) = cos((ind_temp.* a1)+a2)+a3;



    % Put the calibrated parameters into the Arctic annual temperature cycle parameterisation
    calibrated_monthlytemp = f.*(cos(x .* g - e .* exp(cos(x.^(a)))) + ((AAMST - mean(f.*(cos(x .* g - e .* exp(cos(x.^(a))))), 2, 'omitnan')) ./ f));


    calibrated_monthlytemp_save{j} = calibrated_monthlytemp;


    % Save the calibration parameters    
    parameters_store{j} = [f1,f2,g1,g2,a1,a2,a3];

end
parameters_store = cell2mat(parameters_store');
parameters_store_2300 = parameters_store;
toc






%% Stage ii: AAT - AMST: TEST CALIBRATION PARAMETERS 


% Get calibration parameters from 2100 calibration
parameter_store_2100 = [0.4941   -8.4259   -0.0063    0.9309   -0.0780   -0.5773    0.1017
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


parameter_store_2100_extraced = [parameter_store_2100(1:3, :); parameter_store_2100(6, :); parameters_store_2300(5,:); parameter_store_2100(12, :)];


% Rename 2300 variables so as not to write over previous vals
AAT_emulation_absolute_2300_extracted = AAT_emulation_absolute_2300;
tas_amst_combined_2300_ruby_extracted = tas_amst_combined_2300_ruby;


% Re-calibrated data to extract IPSL
% parameter_store_2100_extraced = parameters_store_2300;



%% ________ Apply AMST calibration parameters to all years ________ %% 


clc
x = linspace(0,2*pi,13);
AMST_calibration_2300 = [];
years = [1850:2100];

for j = [1:6]
    hnew = [];

    % Emulated AAT
    AAMST = AAT_emulation_absolute_2300_extracted{j}';
    AAMST = mean(tas_amst_combined_2300_ruby_extracted{j},2);
    

    e = 0.3;

    % Calibrated Regression Parameters
    f1 = parameter_store_2100_extraced(j,1);
    f2 = parameter_store_2100_extraced(j,2);
    g1 = parameter_store_2100_extraced(j,3);
    g2 = parameter_store_2100_extraced(j,4);
    a1 = parameter_store_2100_extraced(j,5);
    a2 = parameter_store_2100_extraced(j,6);
    a3 = parameter_store_2100_extraced(j,7);

   

    % Calculate the evolution of the parameters with AAST 
    f = (f1.*AAMST) + f2;
    g = (g1.*AAMST) + g2;
    a = cos((AAMST.*a1)+a2)+a3;



    % Cap the linear regressions of the calibration parameters so they don't change after the AAST rises above 6 degrees celsius
    ind_temp = 6;
    indf = AAMST > ind_temp;
    indf = find(indf==1);
    f(indf) = (f1.*ind_temp) + f2;
    g(indf) = (g1.*ind_temp) + g2;
    a(indf) = cos((ind_temp.* a1)+a2)+a3;


    % Put the calibrated parameters into the Arctic annual temperature cycle parameterisation
    AMST_calibration_hold = f.*((cos((x .* g) - e .* exp(cos(x.^(a)))))+((AAMST - mean(f.*((cos((x .* g) - e .* exp(cos(x.^(a)))))), 2, 'omitnan')) ./ f));


    AMST_calibration_2300{j} = AMST_calibration_hold(:,1:12);


end

% Apply bias correction to correct for the pre-industrial offset
calls_for_mean_calibrated = [];
calls_for_mean_cmip6 = [];
AMST_calibration_bias_corrected_2300 = [];

mean_correction_year_index = 50;
for j = 1:6
    BC_cmip6 = tas_amst_combined_2300_ruby_extracted{j}(1:mean_correction_year_index,:);     % Extract the 1850:1900 (pre-industrial) temperature 
    calls_for_mean_cmip6 = cat(1, calls_for_mean_cmip6, BC_cmip6);  % concatenate pre-industrial tas from all models

    calls = AMST_calibration_2300{j}(1:mean_correction_year_index,:);  % Extract our emulated 1850:1900 (pre-industrial) temperature
    calls_for_mean_calibrated = cat(1, calls_for_mean_calibrated, calls);                 % concatenate temperatures for comparison with CMIP6
end

calls_for_mean_calibrated = mean(calls_for_mean_calibrated);      % Calculate the pre-industrial mean
calls_for_mean_cmip6 = mean(calls_for_mean_cmip6);
residuals = calls_for_mean_cmip6 - calls_for_mean_calibrated;     % Subtract the mean CMIP6 pre-industrial from our emulation to generate residuals

for j = 1:6
    AMST_calibration_bias_corrected_2300{j} = AMST_calibration_2300{j} + residuals;       % Add residuals to our calibrated AMST
end




% Re-arrange the shape of our emulated AMST for plotting
tic
tas_cal_rearrange_final_2300 = [];
tas_cmip6_rearrange_final_2300 = [];
for j = 1:6
    disp('NEW MODEL')

    for i = 1:12
        tas_cal_rearrange_final_2300{j,i} = AMST_calibration_bias_corrected_2300{j}(:,i)';
        tas_cmip6_rearrange_final_2300{j,i} = tas_amst_combined_2300_ruby_extracted{j}(:,i)';
    end

end
toc

if size(tas_cmip6_rearrange_final_2300{3,1}, 2) < 451
    tas_cmip6_rearrange_final_2300_hold = cellfun(@(v) [v, NaN],  tas_cmip6_rearrange_final_2300(3,:), 'UniformOutput', false);
    tas_cmip6_rearrange_final_2300(3,:) = tas_cmip6_rearrange_final_2300_hold;
else
    tas_cmip6_rearrange_final_2300_hold = cellfun(@(v) v(1:451),  tas_cmip6_rearrange_final_2300(3,:), 'UniformOutput', false);
    tas_cmip6_rearrange_final_2300(3,:) = tas_cmip6_rearrange_final_2300_hold;
end


cellfun(@size, tas_cmip6_rearrange_final_2300, 'UniformOutput', false)
cellfun(@size, tas_cal_rearrange_final_2300, 'UniformOutput', false)












%% Step iii: Calibrate SIA_max Parameterisation
 
params_sia_max = [];
counter = 0;
x = linspace(0,2*pi,13);    % Parameterisation is a cosine curve so we define the appropriate x value to go with this
mov_mean_ind = 30;
SIA_max_param_1850_extraced_2300 = [];
index_run = [1:5];     % Models to calibrate over
plot_run = 1;

close all
for j = index_run


%   Counter for plotting
    counter = counter + 1;
    

%   CMIP6 1850 SIA
    SIA_max_cmip6_running_mean = siarean_final{j};
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


    % Fminsearch
    opt_iter_options = optimset('Display','iter','TolX',1e-3,'MaxIter', 10e9);  % run calibration function
    [ppp, fval] = fminsearch(ff, p0, opt_iter_options);


    % Assign calibrated parameters
    f = ppp(1);
    a = ppp(2);
    e = ppp(3);
    d = ppp(4);
    g = ppp(5);
   


%   Run function with opimtisation parameters
    e = (a .* 0.8003) + 1.5016;
    y_SIA_max = f .* (-exp(sin((x.^g) .* a - e))) + d; 
    SIA_max_param_1850_extraced_2300 = cat(1, SIA_max_param_1850_extraced_2300, y_SIA_max);



    % Save calibrated parameters
    params_sia_max = cat(1, params_sia_max, [f,a,e,d]);


    
    % Test calibration for each model
    if plot_run == 1

        figure(2 + j)
        set(gcf, 'Units', 'Inches', 'Position', [.5 .5 20 11])
    
    %   Plot optimised values
        ss = plot(x, y_SIA_max, 'Color', [1 0 0], 'linewidth', 3);
        hold on
        plot(x, SIA_max_cmip6, '--', 'Color', 'k', 'linewidth', 3);


    %   Label X axis
        xlabel('Month')

    %   Label Y axis
        ylabel('SIA (km^2)')      

    %   Add title with calibrated parameters
        title(['SIA_max Emulation: ', '#', num2str(j), ', ', num2str(SIA_2300_models(j)), ', '], 'fontsize', 28, 'interpreter', 'none')

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
end








%% Step iii: Calibrate SIA Parameterisation 
 
close all
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


new_vals_extracted = [new_vals(1:4,:); new_vals(12,:)];
AMST_calibration_bias_corrected_2300_extraced = [AMST_calibration_bias_corrected_2300(1:3), AMST_calibration_bias_corrected_2300(5:6)];
tas_amst_2300_extraced = [tas_amst_combined_2300_ruby(1:3), tas_amst_combined_2300_ruby(5:6)];

% Opimisation Initiialisation
Save_sia_cal_params = [];
x_tas_hold = [];
tt_store = 1:12;
index_run = [1:12];

close
hnew3 = [];
counter = 0;
clear lines; cols = lines(12);

for j = 1:5

%   Counter for plotting
    counter = counter + 1;

%   CMIP6 SIA
    SIA_store = siarean_final{j}(1:251,:);


%   Emulated AMST
    x_tas_store = AMST_calibration_bias_corrected_2300_extraced{j}(1:251,:);


%   Rename everything for plotting later
    SIA = SIA_store;
    x_tas = x_tas_store;


%   Use parameterised SIA_max (SIA in 1850) - TAKEN STRAIGHT FROM THE IPSL SIAREAN DATA FOR NOW- NEED TO CALIBRATE PROPERLY SOON
    SIA_max = SIA_max_param_1850(j,1:12);


    % Starting test parameters
    a = new_vals_extracted(j,1);
    sm_shift = new_vals_extracted(j,2);
    w1_store = new_vals_extracted(j,3);
    b = new_vals_extracted(j,4);


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



    % TEST CALIBRATIONS
    cmip_sens = [];
    ind = [1,3,5,6,12];
    set(gcf, 'Units', 'Inches', 'Position', [.5 .5 20 11])
    for i = ind
        cmip_sens_hold = scatter(tas_amst_2300_extraced{j}(1:251,i), siarean_final{j}(1:251,i), 15, 'markerfacecolor', cols(i,:), 'markeredgecolor', cols(i,:)); hold on
        cmip_sens = cat(1, cmip_sens, cmip_sens_hold);

        plot(x_tas_store(:,i), calibrated_SIA(:,i),  'color', cols(i,:), 'linewidth', 1.5); hold on
    end
    legend(cmip_sens, month_label(ind))
end

Save_sia_cal_params = cell2mat(Save_sia_cal_params');
IPSL_SIA_params = Save_sia_cal_params;













%% Apply SIA Calibration Parameters to our SIA Parameterisation 

       
% Final calibrated parameters
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
            0.3848    0.0896    0.1904   -1.8352]; %    CMIP6 Calibrated


% EXTRACT ONLY NEEDED MODELS
model_ind = [1, 2, 3, 6, 12]; % Need to calibrate to IPSL and need to remove CanESM5 (4) when applying to SIA

% AMST Extract
% tas_amst_2300_extraced = [tas_amst_2300(1:3), tas_amst_2300(5:6)];
if size(tas_amst_combined_2300_ruby{3})<451
    tas_amst_combined_2300_ruby{3} = [tas_amst_combined_2300_ruby{3}; repelem(NaN, 12)];
end
tas_amst_2300_extraced = [tas_amst_combined_2300_ruby(1:3), tas_amst_combined_2300_ruby(5:6)];
AMST_calibration_bias_corrected_2300_extraced = [AMST_calibration_bias_corrected_2300(1:3), AMST_calibration_bias_corrected_2300(5:6)];

% SIA Extract
SIA_max_param_1850_extraced = SIA_max_param_1850; 
new_vals_extracted = [Save_sia_cal_params(1:2,:); new_vals(3,:); Save_sia_cal_params(4,:); new_vals(12,:)];

if size(siarean_final{3})<451
    siarean_final{3} = [siarean_final{3}; repelem(NaN, 12)];
end
if size(tas_amst_2300_extraced{3})<451
    tas_amst_2300_extraced{3} = [tas_amst_2300_extraced{3}; repelem(NaN, 12)];
end
tas_global_2300_anom_ruby_extracted = [tas_global_2300_anom_ruby(1:3), tas_global_2300_anom_ruby(5:6)];
        

% Opimisation Initiialisation
calibrated_SIA_store_2300 = [];
x_tas_hold = [];
Save_sia_cal_params_new_2300 = [];
save_x_tas_weighted = [];
clear lines; cols = lines(12);

clc
close all
for j = 1:5
    
    hnew3 = []; 
    tt_store = 1:12;

    
%   SIA CMIP6
    SIA_store = siarean_final{j};
    
    
%   Emulated tas
    x_tas_store = AMST_calibration_bias_corrected_2300_extraced{j};
    % x_tas_store = tas_amst_2300_extraced{j};          % Test: force 2100 calibrations with CMIP6 2300 AMST data to see the difference between the AMST overestimation and no overesimtation
    

%   Rename everything for plotting later
    SIA = SIA_store;
    x_tas = x_tas_store;

%   Use parameterised SIA_max
    SIA_max = SIA_max_param_1850_extraced_2300(j,1:12);

    
    % Calibrated parameters
    a = new_vals_extracted(j,1);
    sm_shift = new_vals_extracted(j,2);
    w1_store = new_vals_extracted(j,3);
    b = new_vals_extracted(j,4);
    

    Save_sia_cal_params_new_2300{j} = [a, sm_shift, w1_store, b];
            
    
%   Extraxt dec value before start date
    tas_hold = x_tas_store(end-length(x_tas_store)+1,:);
    dec_previous_year_tas = tas_hold(12);
    dec_tas_test = x_tas(:,12);
    dec_tas_test = circshift(dec_tas_test,1);
    x_tas_shift = [dec_tas_test, x_tas];
    x_tas_shift(1) = dec_previous_year_tas;


    % The Dependency of SIA on Temperature - section 2.2.3
    x_tas_hold = [];
    for tt = 1:12 

        if tt == 1
            hold_x_tas = ((x_tas(:,tt) .* w1_store) + (x_tas_shift(:,1) .* (1-w1_store))) ./ (w1_store+(1-w1_store));
        else
            hold_x_tas = ((x_tas(:,tt) .* w1_store) + (x_tas(:,tt-1) .* (1-w1_store))) ./ (w1_store+(1-w1_store));
        end

        x_tas_hold(:,tt) = hold_x_tas;
    end
    x_tas = x_tas_hold;     % Updated AMST
    save_x_tas_weighted{j} = x_tas;
            
    
    % Run the equation with optimised parameters
    calibrated_SIA = ( ((SIA_max+sm_shift) .* (1 + exp((x_tas(1,:)-b)))) - (a.*(x_tas - x_tas(1,:))) ) ./ ( 1 + exp(x_tas-b) );
    calibrated_SIA_store_2300{j} = calibrated_SIA;

        
    % TEST CALIBRATIONS
    % close  
    figure(42)
    set(gcf, 'Units', 'Inches', 'Position', [.5 .5 20 11])
    % ind = [1:12];
    ind = [12, 3, 5];
    cmip6_sens = [];
    for i = ind
        cmip6_sens_hold = scatter( movmean(tas_amst_2300_extraced{j}(:,i), 1), movmean(siarean_final{j}(:,i), 1), 20, 'markerfacecolor', cols(i,:), 'markeredgecolor', cols(i,:)); hold on
        cmip6_sens = cat(1, cmip6_sens, cmip6_sens_hold);

        plot(x_tas_store(:,i), calibrated_SIA(:,i), 'color', cols(i,:), 'linewidth', 3); hold on
  
    end
    legend(cmip6_sens, month_label(ind))
    
end


% Reshape SIA for plotting
tic
sia_cmip6_rearrange = [];
sia_cmip6_rearrange_final = [];
for j = 1:5
    disp('NEW MODEL')
    
    for i = 1:12
        sia_cal_rearrange{j,i} = calibrated_SIA_store_2300{j}(:,i)';
        sia_cmip6_rearrange{j,i} = siarean_final{j}(:,i)';
    end
end
sia_cal_rearrange_final_2300 = sia_cal_rearrange;
sia_cmip6_rearrange_final_2300 = sia_cmip6_rearrange;
toc


cellfun(@size, sia_cal_rearrange_final_2300, 'UniformOutput', false)
cellfun(@size, sia_cmip6_rearrange_final_2300, 'UniformOutput', false)






%% Plot SIA with CMIP6 for comparison Fig3.7 and FigB.3 

threshold = [5 95];        % Likely range
counter = 0;
SSP_runs = [1];

close  
figure(42)
set(gcf, 'Units', 'Inches', 'Position', [.5 .5 20 11])

sia_cmip6 = sia_cmip6_rearrange_final_2300;
emul = sia_cal_rearrange_final_2300;

for i = [1:12]
    counter = counter + 1;
    h = subplot(3,4,i);     


    % Mean Calibration 
    cal_mean = mean(cell2mat(emul(:,i)), 'omitnan');

    % CMIP6 DATA 
    cmip6_mean = mean(cell2mat(sia_cmip6(:,i)), 'omitnan');


    % Define years of timeseries
    x = 1850:2300;


    % Percentile: Calibration
    percentile_cal = prctile(cell2mat(emul(:,i)), threshold);

    % Percentile: CMIP6
    percentile_cmip6 = prctile(cell2mat(sia_cmip6(:,i)), threshold);


    % Create curves: Calibration
    curve1 = percentile_cal(2,:);
    curve2 = percentile_cal(1,:);


    % Create curves: CMIP6 
    curve3 = percentile_cmip6(2,:);
    curve4 = percentile_cmip6(1,:);


    % Shading std area: Calibration
    x2 = [x, fliplr(x)];
    inBetween_cal = [curve1, fliplr(curve2)];
    inBetween_cmip6 = [curve3, fliplr(curve4)];


    % Fill: Calibration
    c = colmat(1,:);
    f(1) = fill(x2, inBetween_cal, c, 'edgecolor','none');
    set(f(1),'facealpha',.2)
    hold on


    % Fill: CMIP6
    d = [0 0 0]+0.1;
    f2(1) = fill(x2, inBetween_cmip6, d, 'edgecolor','none');
    set(f2(1),'facealpha',.3)
    hold on


%   Plot Means SSP
    plt(1) = plot(x, cal_mean, 'color', colmat(1,:), 'LineWidth', 2);
    plt_cmip6(1) = plot(x, cmip6_mean, 'color', [0 0 0], 'LineWidth', 2);

    
    % Set size of axis labels
    set(gca,'FontSize', 19)

    % Title to each panel (each representing a month)
    title([month_label(i)], 'fontsize', 22);


    % Main plot title
    % sgtitle(['2100 CMIP6 Calibrations of AMST vs SIA on 2300 Data'], 'fontsize', 28);
    % sgtitle(['2100 CMIP6 Calibrations of AMST vs SIA on 2300 Data (using raw AMST data)'], 'fontsize', 28);


    h.Position(1) = h.Position(1) - 0.052;


%   Y Label
    if ismember(i, [1,5,9])
        ylabel({'Arctic SIA'; '(million km^{2})'})
    end

%   X Label
    if ismember(i, [9:12])
        xlabel({'Year'})
    end    

    % Axis limits
    xlim([1850 2300])
    ylim([0 18])
    yline(1)

    grid
    
    xline(2100, 'linewidth', 1.5)


end

% Legend data
hnew_final = [f, f2, plt, plt_cmip6]; 

% Legend labels
legend_deeets = ["Emulated Likely Range", "CMIP6 Likely Range","Emulated Mean", "CMIP6 Mean"];

 % Create legend
[tt2] = legend(hnew_final, legend_deeets, 'Location', 'eastoutside', 'fontsize', 16, 'NumColumns', 1);
tt2.Box = 'off';

% Move posistion of the legend
newPosition = [0.9 0.15 0.05 0.05];
newUnits = 'normalized';
set(tt2,'Position', newPosition,'Units', newUnits);


% % Save
% temp=[' 2100 SIA on 2300 Data ','.png']; 
% saveas(gca,temp);



set(gcf, 'PaperOrientation', 'landscape')
set(gcf,'PaperSize',[53 29.7000]);

% % Save 
% temp=['Fig3.7', '.pdf']; 
% saveas(gca,temp); 



% % Save 
% temp=['FigB.3', '.pdf']; 
% saveas(gca,temp); 









%% ____________________ USING 2300 CALIBRATIONS FOR ANALYSIS ____________________ %% 


%% Understanding the impact of the weighting scheme on SIA parameterisation: Fig3.9, Section 3.3.2.1

close; clc
figure(42)
set(gcf, 'Units', 'Inches', 'Position', [.5 .5 20 11])

ind_months = [10:12, 2:5];
ind_months = 5;

% WITH WEIGHTING
h = subplot(1,2,2);
for j = 3
    for i = ind_months
        x_tas_store = AMST_calibration_bias_corrected_2300_extraced{j};
        x_tas = x_tas_store;
        
        
        
        a = new_vals_extracted(j,1);
        sm_shift = new_vals_extracted(j,2);
        w1_store = new_vals_extracted(j,3);
        b = new_vals_extracted(j,4);
        
        
        
        %   Extraxt dec value before start date
        tas_hold = x_tas_store(end-length(x_tas_store)+1,:);
        dec_previous_year_tas = tas_hold(12);
        dec_tas_test = x_tas(:,12);
        dec_tas_test = circshift(dec_tas_test,1);
        x_tas_shift = [dec_tas_test, x_tas];
        x_tas_shift(1) = dec_previous_year_tas;
        
        
        % The Dependency of SIA on Temperature
        x_tas_hold = [];
        for tt = 1:12 
        
            if tt == 1
                hold_x_tas = ((x_tas(:,tt) .* w1_store) + (x_tas_shift(:,1) .* (1-w1_store))) ./ (w1_store+(1-w1_store));
            else
                hold_x_tas = ((x_tas(:,tt) .* w1_store) + (x_tas(:,tt-1) .* (1-w1_store))) ./ (w1_store+(1-w1_store));
            end
        
            x_tas_hold(:,tt) = hold_x_tas;
        end
        x_tas_new = x_tas_hold;     % Updated AMST
        
        
        
        
        
        
        % 1st month
        SIA_max = SIA_max_param_1850_extraced_2300(j,i);
        x_tas = x_tas_new(:,i);
        calibrated_SIA = ( ((SIA_max+sm_shift) .* (1 + exp((x_tas(1,:)-b)))) - (a.*(x_tas - x_tas(1,:))) )  ./ ( 1 + exp(x_tas-b) ); 
        scatter( x_tas_store(:,i), calibrated_SIA, 30, 'markerfacecolor', colorblind(i,:), 'markeredgecolor', colorblind(i,:) ); hold on
        
    end
end
grid
h.Position(3) = h.Position(3) + 0.05;
h.Position(4) = h.Position(4) - 0.25;

ylabel({ 'Arctic SIA (million km^{2})'} )
xlabel( 'Arctic Temperature (\circC)' )

title(' b) SIA Parameterisation with Weighting Scheme' )

legend(month_label(ind_months)); legend boxoff 

ax = gca;
ax.FontSize = 16;



% NO WEIGHTING
h = subplot(1,2,1);
for j = 3
    for i = ind_months

        x_tas_store = AMST_calibration_bias_corrected_2300_extraced{j}; 
        x_tas = x_tas_store(:,i);
        
        
        
        a = new_vals_extracted(j,1);
        sm_shift = new_vals_extracted(j,2);
        w1_store = new_vals_extracted(j,3);
        b = new_vals_extracted(j,4);
        
       

        % 1st month
        SIA_max = SIA_max_param_1850_extraced_2300(j,i);
        calibrated_SIA = ( ((SIA_max+sm_shift) .* (1 + exp((x_tas(1,:)-b)))) - (a.*(x_tas - x_tas(1,:))) )  ./ ( 1 + exp(x_tas-b) ); 
        scatter( x_tas_store(:,i), calibrated_SIA, 30, 'markerfacecolor', colorblind(i,:), 'markeredgecolor', colorblind(i,:) ); hold on

        
    end
end
grid
h.Position(1) = h.Position(1) - 0.03;
h.Position(3) = h.Position(3) + 0.05;
h.Position(4) = h.Position(4) - 0.25;

ylabel({'Arctic SIA (million km^{2})'})
xlabel('Arctic Temperature (\circC)')

title(' a) SIA Parameterisation without Weighting Scheme' )

legend(month_label(ind_months)); legend boxoff 

ax = gca;
ax.FontSize = 16;


set(gcf, 'PaperOrientation', 'landscape')
set(gcf,'PaperSize',[53 29.7000]);

% % Save Normalised
% temp=['Fig3.9', '.pdf']; 
% saveas(gca,temp);







%% Plot against AMST: to assess how well our calibrations match the CMIP6 AMST vs SIA data to 2300 (Fig6, Section 3.1)

tas_global_2300_anom_ruby_extracted = [tas_global_2300_anom_ruby(1:3), tas_global_2300_anom_ruby(5:6)];

if size(calibrated_SIA_store_2300{3}, 1) < 451
    calibrated_SIA_store_2300{3} = [calibrated_SIA_store_2300{3}; repelem(NaN, 12)];
else
    calibrated_SIA_store_2300{3} = calibrated_SIA_store_2300{3}(1:451,:);
end

if size(AMST_calibration_bias_corrected_2300_extraced{3}, 1) < 451
    AMST_calibration_bias_corrected_2300_extraced{3} = [AMST_calibration_bias_corrected_2300_extraced{3}; repelem(NaN, 12)];
else
    AMST_calibration_bias_corrected_2300_extraced{3} = AMST_calibration_bias_corrected_2300_extraced{3}(1:451,:);
end


cols = jet(12);
close all; clc
figure(3)
set(gcf, 'Units', 'Inches', 'Position', [.4 .4 19 10])
ind_months = [12, 3, 5];
counter_made = -1;
for j = 1:5
    counter_made = counter_made + 2;
    h = subplot(5,2,counter_made);
    cmip6_scat = [];
    emulation_scat = [];
    for i = ind_months
        % Emulation
        emulation_scat_hold = plot( AMST_calibration_bias_corrected_2300_extraced{j}(:,i), calibrated_SIA_store_2300{j}(:,i), 'color', cols(i,:), 'linewidth', 1.3 ); hold on
        emulation_scat = cat(1, emulation_scat, emulation_scat_hold);

        % CMIP6
        cmip6_scat_hold = scatter( tas_amst_2300_extraced{j}(:,i), siarean_final{j}(:,i),  5, 'markerfacecolor', cols(i,:), 'markeredgecolor', cols(i,:) ); hold on
        cmip6_scat = cat(1, cmip6_scat, cmip6_scat_hold);

    end
    grid
    h.FontSize = 16;
    title(SIA_2300_models(j), 'FontSize', 14)
    yline(1)
    xlim([-40 15])

    h.Position(3) = h.Position(3) - 0.15;
    h.Position(1) = h.Position(1) - 0.05;

    if j == 1
     h.Position(2) = h.Position(2) + 0.04;
    elseif j == 2
        h.Position(2) = h.Position(2) + 0.03;
    elseif j == 3
        h.Position(2) = h.Position(2) + 0.02;
    elseif j == 4
        h.Position(2) = h.Position(2) + 0.01;
    elseif j == 5
        h.Position(2) = h.Position(2) + 0.00;
    end
end
y_pos = ylabel({'Arctic SIA (million km^{2})'})
y_pos.Position(2) = y_pos.Position(2) + 55;
xlabel('Arctic Temperature (\circC)')





% GMST
ind_months = [12, 5, 3];
counter_made = 0;
for j = 1:5
    counter_made = counter_made + 2;
    h = subplot(5,2,counter_made);
    cmip6_scat = [];
    emulation_scat = [];

    for i = ind_months

        % Emulation to 2300 from 2300 GMST
        emulation_scat_hold = plot( tas_global_2300_anom_ruby_extracted{j}, calibrated_SIA_store_2300{j}(:,i), 'color', cols(i,:), 'linewidth', 1.3 ); hold on
        emulation_scat = cat(1, emulation_scat, emulation_scat_hold);

        % CMIP6
        cmip6_scat_hold = scatter( tas_global_2300_anom_ruby_extracted{j}, siarean_final{j}(:,i),  5, 'markerfacecolor', cols(i,:), 'markeredgecolor', cols(i,:) ); hold on
        cmip6_scat = cat(1, cmip6_scat, cmip6_scat_hold);

    end
    xlim([-2 12])
    grid
    h.FontSize = 16;
    title(SIA_2300_models(j), 'FontSize', 14)
    yline(1)
    % ylabel({'Arctic SIA'; '(million km^{2})'})
    % xlabel(' Global-Mean Temperature Anomaly ')

    h.Position(3) = h.Position(3) - 0.15;
    h.Position(1) = h.Position(1) - 0.25;

    if j == 1
     h.Position(2) = h.Position(2) + 0.04;
    elseif j == 2
        h.Position(2) = h.Position(2) + 0.03;
    elseif j == 3
        h.Position(2) = h.Position(2) + 0.02;
    elseif j == 4
        h.Position(2) = h.Position(2) + 0.01;
    elseif j == 5
        h.Position(2) = h.Position(2) + 0.00;
    end
end
xlabel(' Global-Mean Temperature Anomaly ')

legend_text = [month_label(ind_months); month_label(ind_months)];
legend_text(1,:) = cellfun(@(v) [v, ': CMIP6'], legend_text(1,:), 'UniformOutput', false);
legend_text(2,:) = cellfun(@(v) [v, ': Emulator'], legend_text(2,:), 'UniformOutput', false);
legend_text = legend_text';
legend_text = legend_text(:);

legend_data = [cmip6_scat, emulation_scat];
legend_data = legend_data(:);

tt2 = legend(legend_data, legend_text, 'location', 'southoutside', 'NumColumns', 1);
tt2.Position(1) = tt2.Position(1)+0.17;
tt2.Position(2) = tt2.Position(2)+0.1;
tt2.FontSize = 16;
legend boxoff


x = [0.29 0.29];
y = [0.02 0.98];
annotation('line', x, y, 'LineStyle', '--', 'Linewidth', 2)


set(gcf, 'PaperOrientation', 'portrait')
set(gcf,'PaperSize',[45 27.7000]); 














%% RATE OF SEASONAL ARCTIC SEA ICE LOSS WITH SEASONAL AA TO 2300: Fig3.6, Chapter 3, Section 3.3.1.3


tas_aat_extracted = [aat_2300_anom_ruby(1:3); aat_2300_anom_ruby(5:6)];
tas_gmst_extracted = [tas_global_2300_anom_ruby(1:3), tas_global_2300_anom_ruby(5:6)];
tas_amst_2300_extracted = [tas_amst_combined_2300_ruby(1:3), tas_amst_combined_2300_ruby(5:6)];

if size(tas_gmst_extracted{1},1) ~= 1
    tas_gmst_extracted = cellfun(@transpose, tas_gmst_extracted, 'UniformOutput', false)';
end
if size(tas_aat_extracted{1},1) ~= 1
    tas_aat_extracted = cellfun(@transpose, tas_aat_extracted, 'UniformOutput', false)';
end

% GET MEAN SIA (reshape)
reshaped_siarean_final = [];
reshaped_amst_2300 = [];
for j = 1:5
    for i = 1:12
        reshaped_siarean_final{j,i} = siarean_final{j}(:,i)';
        reshaped_amst_2300{j,i} = tas_amst_2300_extracted{j}(:,i)';
    end
end
% GET AMST ANOMALY
pi_amst = cellfun(@(v) mean(v(1:50)),  reshaped_amst_2300,  'UniformOutput', false);
reshaped_amst_2300 = cellfun(@(v1,v2) v1-v2,  reshaped_amst_2300,  pi_amst, 'UniformOutput', false);


tas_store_2100 = tas_store{1};
pi_amst_2100 = cellfun(@(v) mean(v(1:50,:)),  tas_store_2100,  'UniformOutput', false);
tas_store_2100 = cellfun(@(v1,v2) v1-v2,  tas_store_2100,  pi_amst_2100, 'UniformOutput', false);



left_color = [1 0 0];
right_color = [0 0 1];
colorss1 = [1 0.7 0.7; 0.8 0.8 1; 0 0.6906 0.5];
close all; clc
figure(2)
set(gcf, 'Units', 'Inches', 'Position', [.4 .4 19 10])
for j = 1:5
    for i = 1:12
        h = subplot(3,4,i);

        yyaxis left 
        sia = reshaped_siarean_final{j,i};
        plot(1850:2300-1, movmean(diff(sia), 40), '-', 'color', colorss1(1,:), 'linewidth', 1 ); hold on
        plot(1850:2300-1, movmean(diff(mean(cell2mat(reshaped_siarean_final(:,i)))), 40), '-', 'color', [1 0 0], 'linewidth', 2 ); hold on
        if ismember(i, [1,5,9])
            ylabel('Rate of SIA loss')
        end
        ylim([-0.3 0])



        yyaxis right
        amst = reshaped_amst_2300{j,i};
        gmst = tas_gmst_extracted{j};
        aat = tas_aat_extracted{j};
        plot(1850:2300, movmean( amst ./ gmst , 25), '-', 'color', colorss1(2,:), 'linewidth', 1 ); hold on
        plot(1850:2300, movmean( mean(cell2mat(reshaped_amst_2300(:,i))) ./ mean(cell2mat(tas_gmst_extracted)) , 25), '-', 'color', [0 0 1], 'linewidth', 2 ); hold on

        if ismember(i, [11])
            ylim([2 7])
        elseif ismember(i, [1, 2, 10, 12])
            ylim([2 6])
        elseif ismember(i, [3:4])
            ylim([2 4])
        elseif ismember(i, [6:8])
            ylim([0 2])
        elseif ismember(i, [9])
            ylim([1 3])
        elseif ismember(i, [5])
            ylim([1 2])
        end
        % ylim([2 7])
        if ismember(i, [4,8,12])
            ylabel(' Seasonal AA ')
        end



        xlim([1980 2300])
        if ismember(i, 9:12)
            xlabel( 'Year' )       
        end

        grid 

        title(month_label(i))

        ax = gca;
        ax.YAxis(1).Color = 'r';
        ax.YAxis(2).Color = 'b';

        ax = gca; ax.FontSize = 16; 
    end
end


set(gcf, 'PaperOrientation', 'landscape')
set(gcf,'PaperSize',[53 29.7000]);

% % Save
% temp=['Fig3.6', '.pdf']; 
% saveas(gca,temp); 





%% RATE OF ANNUAL ARCTIC SEA ICE LOSS WITH ANNUAL AA TO 2300: Fig3.5, Chapter 3, Section 3.3.1.3

tas_aat_extracted = [aat_2300_anom_ruby(1:3); aat_2300_anom_ruby(5:6)];
tas_gmst_extracted = [tas_global_2300_anom_ruby(1:3), tas_global_2300_anom_ruby(5:6)];
tas_amst_2300_extracted = [tas_amst_combined_2300_ruby(1:3), tas_amst_combined_2300_ruby(5:6)];

if size(tas_gmst_extracted{1},1) ~= 1
    tas_gmst_extracted = cellfun(@transpose, tas_gmst_extracted, 'UniformOutput', false)';
end
if size(tas_aat_extracted{1},1) ~= 1
    tas_aat_extracted = cellfun(@transpose, tas_aat_extracted, 'UniformOutput', false)';
end

siarean_annual = cellfun(@(v) mean(v,2),  siarean_final,  'UniformOutput', false);
if size(siarean_annual{1},1) ~= 1
    siarean_annual = cellfun(@transpose,  siarean_annual,  'UniformOutput', false)';
end


left_color = [1 0 0];
right_color = [0 0 1];
colorss1 = [1 0.7 0.7; 0.8 0.8 1; 0 0.6906 0.5];
close all; clc
figure(2)
set(gcf, 'Units', 'Inches', 'Position', [.4 .4 19 10])
for j = 1:5
    h = subplot(1,1,1);

    yyaxis left 
    sia = siarean_annual{j};
    plot(1850:2300-1, movmean(diff(sia), 40), '-', 'color', colorss1(1,:), 'linewidth', 1 ); hold on
    plot(1850:2300-1, movmean(diff(mean(cell2mat(siarean_annual))), 40), '-', 'color', [1 0 0], 'linewidth', 2 ); hold on
    ylabel('Rate of SIA Annual Loss')
    ylim([-0.16 0])



    yyaxis right
    gmst = tas_gmst_extracted{j};
    aat = tas_aat_extracted{j};
    plot(1850:2300, movmean( aat ./ gmst , 25), '-', 'color', colorss1(2,:), 'linewidth', 1 ); hold on
    plot(1850:2300, movmean( mean(cell2mat(tas_aat_extracted)) ./ mean(cell2mat(tas_gmst_extracted)) , 25), '-', 'color', [0 0 1], 'linewidth', 2 ); hold on

    ylim([1.5 4])
    ylabel(' Annual AA ')


end
xlim([1980 2300])
xlabel( 'Year' )   
grid 
ax = gca;
ax.YAxis(1).Color = 'r';
ax.YAxis(2).Color = 'b';
ax = gca; ax.FontSize = 18; 
set(gcf, 'PaperOrientation', 'landscape')
set(gcf,'PaperSize',[53 29.7000]);

% % Save
% temp=['Fig3.5', '.pdf']; 
% saveas(gca,temp); 










%% ____________________________________________________________________________________________________________________________________________________________ % 















%% UNDERSTANDING THE SIA PARAMETERISATION FOR CHAPTER 3: WITH WEIGHTING

j = 3; % Use as an example model
x_tas_store = AMST_calibration_bias_corrected_2300_extraced{j};
x_tas = x_tas_store;



a = new_vals_extracted(j,1);
sm_shift = new_vals_extracted(j,2);
w1_store = new_vals_extracted(j,3);
b = new_vals_extracted(j,4);



%   Extraxt dec value before start date
tas_hold = x_tas_store(end-length(x_tas_store)+1,:);
dec_previous_year_tas = tas_hold(12);
dec_tas_test = x_tas(:,12);
dec_tas_test = circshift(dec_tas_test,1);
x_tas_shift = [dec_tas_test, x_tas];
x_tas_shift(1) = dec_previous_year_tas;


% The Dependency of SIA on Temperature - section 2.2.3
x_tas_hold = [];
for tt = 1:12 

    if tt == 1
        hold_x_tas = ((x_tas(:,tt) .* w1_store) + (x_tas_shift(:,1) .* (1-w1_store))) ./ (w1_store+(1-w1_store));
    else
        hold_x_tas = ((x_tas(:,tt) .* w1_store) + (x_tas(:,tt-1) .* (1-w1_store))) ./ (w1_store+(1-w1_store));
    end

    x_tas_hold(:,tt) = hold_x_tas;
end
x_tas_new = x_tas_hold;     % Updated AMST






% 1st month
i = 5;
SIA_max = SIA_max_param_1850_extraced_2300(j,i);
x_tas = x_tas_new(:,i);
close; calibrated_SIA = ( ((SIA_max+sm_shift) .* (1 + exp((x_tas(1,:)-b)))) - (a.*(x_tas - x_tas(1,:))) )  ./ ( 1 + exp(x_tas-b) ); plot(x_tas_store(:,i), calibrated_SIA); hold on



% 2nd month
i = 3;
x_tas = x_tas_new(:,i); 
SIA_max = SIA_max_param_1850_extraced_2300(j,i);

calibrated_SIA = ( ((SIA_max+sm_shift) .* (1 + exp((x_tas(1,:)-b)))) - (a.*(x_tas - x_tas(1,:))) )  ./ ( 1 + exp(x_tas-b) ); plot(x_tas_store(:,i), calibrated_SIA); hold on


% 3rd month
i = 12;
x_tas = x_tas_new(:,i);
SIA_max = SIA_max_param_1850_extraced_2300(j,i);

calibrated_SIA = ( ((SIA_max+sm_shift) .* (1 + exp((x_tas(1,:)-b)))) - (a.*(x_tas - x_tas(1,:))) )  ./ ( 1 + exp(x_tas-b) ); plot(x_tas_store(:,i), calibrated_SIA); hold on


% 4th month
i = 4;
x_tas = x_tas_new(:,i);
SIA_max = SIA_max_param_1850_extraced_2300(j,i);

calibrated_SIA = ( ((SIA_max+sm_shift) .* (1 + exp((x_tas(1,:)-b)))) - (a.*(x_tas - x_tas(1,:))) )  ./ ( 1 + exp(x_tas-b) ); plot(x_tas_store(:,i), calibrated_SIA); hold on


% 5th month
i = 10;
x_tas = x_tas_new(:,i);
SIA_max = SIA_max_param_1850_extraced_2300(j,i);

calibrated_SIA = ( ((SIA_max+sm_shift) .* (1 + exp((x_tas(1,:)-b)))) - (a.*(x_tas - x_tas(1,:))) )  ./ ( 1 + exp(x_tas-b) ); plot(x_tas_store(:,i), calibrated_SIA); hold on


yline(8)
yline(10)



















%% UNDERSTANDING THE SIA PARAMETERISATION FOR CHAPTER 3: WITH OUT WEIGHTING

j = 3; i = 5;
x_tas = AMST_calibration_bias_corrected_2300_extraced{j}(:,i);
SIA_max = SIA_max_param_1850_extraced_2300(j,i);


a = new_vals_extracted(j,1);
sm_shift = new_vals_extracted(j,2);
w1_store = new_vals_extracted(j,3);
b = new_vals_extracted(j,4);

% [a, sm_shift, w1_store, b]

% 1st month
% close; calibrated_SIA = SIA_max ./ ( 1 + exp(x_tas-b) ); plot(x_tas-b, calibrated_SIA); hold on
close; calibrated_SIA = ( ((SIA_max+sm_shift) .* (1 + exp((x_tas(1,:)-b)))) - (a.*(x_tas - x_tas(1,:))) )  ./ ( 1 + exp(x_tas-b) ); plot(x_tas, calibrated_SIA); hold on



% 2nd month
i = 3;
x_tas = AMST_calibration_bias_corrected_2300_extraced{j}(:,i); 
SIA_max = SIA_max_param_1850_extraced_2300(j,i);

calibrated_SIA = ( ((SIA_max+sm_shift) .* (1 + exp((x_tas(1,:)-b)))) - (a.*(x_tas - x_tas(1,:))) )  ./ ( 1 + exp(x_tas-b) ); plot(x_tas, calibrated_SIA); hold on


% 3rd month
i = 12;
x_tas = AMST_calibration_bias_corrected_2300_extraced{j}(:,i);
SIA_max = SIA_max_param_1850_extraced_2300(j,i);

calibrated_SIA = ( ((SIA_max+sm_shift) .* (1 + exp((x_tas(1,:)-b)))) - (a.*(x_tas - x_tas(1,:))) )  ./ ( 1 + exp(x_tas-b) ); plot(x_tas, calibrated_SIA); hold on






















%% GOF of SIA vs AMST 2300


clc
Acmip6 = cellfun(@(v) movmean(v, 20), siarean_final, 'UniformOutput', false);
Aemulated = cellfun(@(v) movmean(v, 20), calibrated_SIA_store_2300, 'UniformOutput', false);
A = cellfun(@(v1,v2) (v1 - v2).^2, Acmip6, Aemulated, 'un', 0);
A2 = cellfun(@(v) sum(v, 1, 'omitnan'), A,  'un', 0);
gof = cellfun(@(v) v ./ 451, A2, 'un', 0)';
gof = cell2mat(gof);
gofmean = mean(gof);
gofmean([3,5,12])

rsquared_save = [];
for j = 1:5
    for i = 1:12
        mdl = fitlm( Acmip6{j}(:,i) , Aemulated{j}(:,i));
        rsquared_save{j,i} = mdl.Rsquared.Ordinary; 
    end
end
rsquared_save = cell2mat(rsquared_save);
rsquared_savemean = mean(rsquared_save);
rsquared_savemean([3,5,12])




%% GOF of SIA vs GMST 2300


clc
Acmip6 = cellfun(@(v) movmean(v, 20), siarean_final, 'UniformOutput', false);
Aemulated = cellfun(@(v) movmean(v, 20), calibrated_SIA_store_2300, 'UniformOutput', false);
A = cellfun(@(v1,v2) (v1 - v2).^2, Acmip6, Aemulated, 'un', 0);
A2 = cellfun(@(v) sum(v, 1, 'omitnan'), A,  'un', 0);
gof = cellfun(@(v) v ./ 451, A2, 'un', 0)';
gof = cell2mat(gof);
gofmean = mean(gof);
gofmean([3,5,12])

rsquared_save = []; 
for j = 1:5
    for i = 1:12
        mdl = fitlm( Acmip6{j}(:,i) , Aemulated{j}(:,i));
        rsquared_save{j,i} = mdl.Rsquared.Ordinary; 
    end
end
rsquared_save = cell2mat(rsquared_save);
rsquared_savemean = mean(rsquared_save);
rsquared_savemean([3,5,12])







%% GOF of SIA vs AMST 2300: when foreced with CMIP6 tas


clc
Acmip6 = cellfun(@(v) movmean(v, 20), siarean_final, 'UniformOutput', false);
Aemulated = cellfun(@(v) movmean(v, 20), calibrated_SIA_store_2300_test, 'UniformOutput', false);
A = cellfun(@(v1,v2) (v1 - v2).^2, Acmip6, Aemulated, 'un', 0);
A2 = cellfun(@(v) sum(v, 1, 'omitnan'), A,  'un', 0);
gof = cellfun(@(v) v ./ 451, A2, 'un', 0)';
gof = cell2mat(gof);
gofmean = mean(gof);
gofmean([3,5,12])

rsquared_save = [];
for j = 1:5
    for i = 1:12
        mdl = fitlm( Acmip6{j}(:,i) , Aemulated{j}(:,i));
        rsquared_save{j,i} = mdl.Rsquared.Ordinary; 
    end
end
rsquared_save = cell2mat(rsquared_save);
rsquared_savemean = mean(rsquared_save);
rsquared_savemean([3,5,12])






%% Calculate sensitivity in each year between Arctic warming and global warming

refs = [1:20:451];
refs = [refs, 451];
clc
sens_yearly_emul_amst = [];
sens_yearly_emul_gmst = [];

sens_yearly_cmip6_amst = [];
sens_yearly_cmip6_gmst = [];
for j = 1:5
    for i = 3

        for ii = 1:length(refs)-1

            % Emul
            diff_amst_emul = AMST_calibration_bias_corrected_2300_extraced{j}(:,i); 
            diff_sia_emul = calibrated_SIA_store_2300{j}(:,i); 
            sens_yearly_emul_amst{j,ii} = polyfit( diff_amst_emul(refs(ii):refs(ii+1)) , diff_sia_emul(refs(ii):refs(ii+1)), 1);  
    
            diff_gmst_emul = movmean(tas_global_2300_anom_ruby_extracted{j}, 30);
            sens_yearly_emul_gmst{j,ii} = polyfit( diff_gmst_emul(refs(ii):refs(ii+1)) , diff_sia_emul(refs(ii):refs(ii+1)), 1); 
    
    
    
            % Emul
            diff_amst_cmip6 = movmean(tas_amst_2300_extraced{j}(:,i), 30); 
            diff_sia_cmip6 = movmean(siarean_final{j}(:,i), 30); 
            sens_yearly_cmip6_amst{j,ii} = polyfit( diff_amst_cmip6(refs(ii):refs(ii+1)) , diff_sia_cmip6(refs(ii):refs(ii+1)), 1); 
    
            diff_gmst_cmip6 = movmean(tas_global_2300_anom_ruby_extracted{j}, 30);
            sens_yearly_cmip6_gmst{j,ii} = polyfit( diff_gmst_cmip6(refs(ii):refs(ii+1)) , diff_sia_cmip6(refs(ii):refs(ii+1)), 1); 

        end
    end
end

sens_yearly_cmip6_gmst = cellfun(@(v) v(1), sens_yearly_cmip6_gmst, 'UniformOutput', true);
sens_yearly_cmip6_amst = cellfun(@(v) v(1), sens_yearly_cmip6_amst, 'UniformOutput', true);
sens_yearly_emul_gmst = cellfun(@(v) v(1), sens_yearly_emul_gmst, 'UniformOutput', true);
sens_yearly_emul_amst = cellfun(@(v) v(1), sens_yearly_emul_amst, 'UniformOutput', true);



gof_save_amst = [];
rsquared_save_amst = [];
for j = 1:5
    % GOF of SIA vs AMST Sensitivity 2300
    clc
    Acmip6 = sens_yearly_cmip6_amst(j,:); 
    Aemulated = sens_yearly_emul_amst(j,:); 
    A = (Acmip6 - Aemulated).^2;
    A2 = sum(A, 'omitnan');
    gof = A2 ./ 451;
    gof_save_amst{j} = gof;
    
    
    mdl = fitlm( Acmip6 , Aemulated); 
    rsquared_save_amst(j) = mdl.Rsquared.Ordinary; 

end
gofmean_amst = mean(gof)
rsquared_savemean_amst = mean(rsquared_save_amst)


gof_save_gmst = [];
rsquared_save_gmst = [];
for i = 3
    % GOF of SIA vs GMST Sensitivity 2300
    Acmip6 = sens_yearly_cmip6_gmst(j,:); 
    Aemulated = sens_yearly_emul_gmst(j,:); 
    A = (Acmip6 - Aemulated).^2;
    A2 = sum(A, 'omitnan');
    gof = A2 ./ 451;
    gof_save_gmst{j} = gof;
    

    mdl = fitlm( Acmip6 , Aemulated); 
    rsquared_save_gmst(j) = mdl.Rsquared.Ordinary; 
    
end
gofmean_gmst = mean(gof)
rsquared_savemean_gmst = mean(rsquared_save_gmst)



