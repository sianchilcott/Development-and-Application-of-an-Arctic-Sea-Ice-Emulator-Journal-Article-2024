%% Function for SIA_max 
% Chapter 2 Section 2.3.4

function obj = SIA_max_Calibration_Publication(ppp, x, SIA_max_cmip6)

f = ppp(1);
a = ppp(2);
e = ppp(3);
d = ppp(4);
g = ppp(5);


% y_SIA_max = f .* -exp(sin(x .* a - e)) + sin(x.^g) + (f*d);
% y_SIA_max = f .* (-exp(sin(x .* a - e))) + sin(x.^g) + (f*d);
% y_SIA_max = f .* (-exp(sin((x.^g) .* a - e))) + (f*d );
y_SIA_max = f .* (-exp(sin((x.^g) .* a - e))) + d;


A = (SIA_max_cmip6 - y_SIA_max).^2;
obj = sum(A);
obj = mean(obj);

end


%%  2.0

% function obj = SIA_max_Calibration_Publication(ppp, x, SIA_max_cmip6)
% 
% f = ppp(1);
% a = ppp(2);
% d = ppp(3);
% g = ppp(4);
% 
% 
% y_SIA_max = f .* (-exp(sin(x .* a))) + sin(x.^g) + (f*d);
% 
% 
% A = (SIA_max_cmip6 - y_SIA_max).^2;
% obj = sum(A);
% obj = mean(obj);
% 
% end
