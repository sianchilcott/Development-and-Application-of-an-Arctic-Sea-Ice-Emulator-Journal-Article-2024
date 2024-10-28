%% Function for Sea Ice Area (SIA) Parameterisation
% Section 2.5

function obj = SIA_CMIP6_Calibration_Publication(SIA_max,x_tas,SIA,dec_previous_year_tas,p)

% Calibration Parmaters
a = p(1);  
b = p(2);
sm_shift = p(3);
w1_store = p(4);



% Get the December temperature of the year before to weight January
dec_tas_test = x_tas(:,12);
dec_tas_test = circshift(dec_tas_test,1);
x_tas_shift = [dec_tas_test, x_tas];
x_tas_shift(1) = dec_previous_year_tas;


% Weight the temperatures (discussed in thesis Chapter 2 Section 2.3.5)
for tt = 1:12 

    if tt == 1
        hold_x_tas = ((x_tas(:,tt) .* w1_store) + (x_tas_shift(:,1) .* (1-w1_store))) ./ (w1_store+(1-w1_store));
    else
        hold_x_tas = ((x_tas(:,tt) .* w1_store) + (x_tas(:,tt-1) .* (1-w1_store))) ./ (w1_store+(1-w1_store));
    end

    x_tas_hold(:,tt) = hold_x_tas;
    x_tas_diff(:,tt) = hold_x_tas - x_tas(:,tt);
end

x_tas = x_tas_hold; 

 
% Calculate SIA
y_sia = ( ((SIA_max+sm_shift) .* (1 + exp((x_tas(1,:)-b)))) - (a.*(x_tas - x_tas(1,:))) ) ./ (1 + exp(x_tas-b));


% Calculate the residual sum of the squared di↵erences (RSS) to find the global miuimum of each parameter
A = (SIA - y_sia).^2;
obj = sum(A);
obj = mean(obj);

end




