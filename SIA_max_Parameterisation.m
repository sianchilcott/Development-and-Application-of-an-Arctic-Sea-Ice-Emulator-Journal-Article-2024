%% Function for SIA_max 
% Section 2.5.1

function obj = SIA_max_Calibration_Publication(ppp, x, SIA_max_cmip6)

f = ppp(1);
a = ppp(2);
e = ppp(3);
d = ppp(4);
g = ppp(5);


% SIA_max parameterisation
y_SIA_max = f .* (-exp(sin((x.^g) .* a - e))) + d;


% Calculate the residual sum of the squared di↵erences (RSS) to find the global miuimum of each parameter
A = (SIA_max_cmip6 - y_SIA_max).^2;
obj = sum(A);
obj = mean(obj);

end

