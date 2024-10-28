%% Function for our Arctic Seasonal Temperature Parameterisation
% Section 2.4

function obj = AMST_parameterisation_publication2(ppp,x,y,AAMST,e)

f1 = ppp(1);
f2 = ppp(2);
g1 = ppp(3);
g2 = ppp(4);
a1 = ppp(5);
a2 = ppp(6);
a3 = ppp(7);
        

f = (f1.*AAMST) + f2;
g = (g1.*AAMST) + g2;
a = cos((AAMST.* a1) +a2) +a3;   

        
y_tas = f.*(cos(x .* g - e .* exp(cos(x.^(a))))+((AAMST - mean(f.*(cos(x .* g - e .* exp(cos(x.^(a))))), 2, 'omitnan')) ./ f));

        
% Residual sum of the squared di↵erences (RSS) to find the global minimum
A = (y - y_tas).^2;
obj = sum(A);
obj = mean(obj);

end
