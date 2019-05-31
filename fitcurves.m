function [best_coeffs,best_fun,AICs,SS,best_residuals] = fitcurves(xdata,ydata,fits,LBC,UBC,IVC,options,varargin)

% FITCURVES test fit of curves to data
%   [C,F,AS,SS,R] = FITCURVES(X,Y,FITS,LB,UB,IV,OPTS) fits the set of functions requested
%   in FITS to find the best function that describes Y = F(X). The best function
%   is the one with the lowest AIC score.
%
%   [...] = FITCURVES(...,'S') will return the best function as the lowest
%   sum-of-squares, ignoring the number of free parameters - useful when
%   all you want is the very best fit!
%   
%   The set of functions is defined by the array FITS, whose numeric entries 
%   request the following:
%       1 = linear      (2 par)         a + b*x
%       2 = exponential (2 par)         a*exp(b*x)  
%       3 = exponential (3 par)         a + b*exp(c*x)
%       4 = double exponential (4 par)  a*exp(b*x) + c*exp(d*x)
%       5 = double exponential (5 par)  a + b*exp(c*x) + d*exp(e*x)
%       6 = triple exponential (6 par)  a*exp(b*x) + c*exp(d*x) + e*exp(f*x)
%       7 = power law  (2 par)          a*x^b
%       8 = power law (3 par)           a + b*x^c
%       9 = inverse second order (3 par)a + b/x + c/x^2
%      10 = Rayleigh (2 par)            a*x*exp(-b*sqrt(x))
%      11 = Gamma function (3 par)  1/(b^a * Gamma(a)) * x^(a-1) * exp(x/b) * c [x>0; set c=1 to get gamma PDF]    
%      12 = piece-wise linear#1 (5 par) 1: a+ b*x (x <= c) 2: d + e*x (x >= c)
%      13 = truncated power law (3 par) a*x^b * exp(x*c)    
%      14 = Pareto function (1 par)     1 - 1/(x^a)     
%      15 = modded power law I (2 par)  a*(1-x)^-b)
%      16 = modded power law II (2 par) a*(1+x)^b
%      17 = piece-wise power law (5 par) 1: a*x^b (x <= c) 2: d*x^e (x > c)
%      18 = inverse cubic polynomical (4 par) a/x^3 + b/x^2 + c/x + d  
%      19 = quadratic (3 par)           a + bx + cx^2 
%      20 = exponential rise-to-max (2par) a*(1-exp(-b*x))
%      21 = exponential rise-to-max (3par) a + b*(1-exp(-c*x))
%      22 = dbl exponential max (4par)  a*(1-exp(-b*x)) + c*(1-exp(-d*x))
%      23 = inverse first-order (2 par) a + b/x 
%      24 = stretched exponential decay (3 par)  a*exp(-x/b)^c
%      25 = cubic (4par):               a + bx + cx^2 + dx^3 
%      26 = quartic (5 par):            a + bx + cx^2 + dx^3 + e*x^4 
%      27 = hyperbola, double (4par):   a*x/(b+x) + c*x/(d+x)
%      28 = hyperbola, double (5par):   a*x/(b+x) + c*x/(d+x) + e*x
%      29 = quintic (6 par):            a + bx + cx^2 + dx^3 + e*x^4 + f*x^5
%      30 = exp+linear (4 par):         a + b*exp(c*x) + d*x
%      31 = rational (4 par):           (a + bx)/(1 + cx + dx^2)
%      32 = rational (5 par):           (a + bx + cx^2)/(1 + dx +ex^2)
%      33 = Gaussian (3 par):           a*exp(-(x-b)^2/2c^2) 
%      34 = double Gaussian (6 par):    a*exp(-(x-b)^2/2c^2) + d*exp(-(x-e)^2/2f^2)
%      35 = log Gaussian (3 par):       a*exp(-(ln(x)-b)^2/2c^2)     
%   
%   The range of values over which the coefficients for the requested functions 
%   are tested can be limited by lower (LB) and upper (UB) bound values.
%   Similarly, the initial values for the search can be given in the cell
%   array IV. The cell arrays LB, UB, and IV must have the same number of elements as FITS.
%   In each cell, there must be as many values as there are coefficients.
%   Any entry can be omitted by putting [] (including the entire array). OPTS is the
%   standard MatLab OPTIONS structure; set [] to omit (e.g. if using flag 'S').
%
%   Returns: C a cell array of the coefficients of the best-fit function; F
%   the inline function that was the best-fit; AS the array of AIC scores for every tested function; SS the
%   sum of squares all fitted functions; R the residuals for the *best*
%   fit function.
%
%   Note#1: TO FIX: random allocation of initial values always in [0 1] if
%   not specified; should be within [LB, UB] if these specified!!
%
%   Note#2: To get e.g. exponential decay, set bounds of LB = {[-Inf -Inf]}, UB =
%   {[Inf 0]}; 
%
%   Mark Humphries 22/6/2009

if any(fits > 35)
    error('Non-existent curve chosen')
end

rtn = 'AIC';
if nargin > 7 & strfind(varargin{1},'S')
    rtn = 'SS';
end

% linear fit
linfun = inline('x(1) + x(2) .* xdata','x','xdata');
% exponential fit (2 par)
expAfun = inline('x(1) .* exp(x(2).*xdata)','x','xdata');
% exponential fit (3 par)
expBfun = inline('x(1) + x(2) .* exp(x(3).*xdata)','x','xdata');
% double exponential fit (4 par)
dblexpAfun = inline('x(1) .* exp(x(2).*xdata) + x(3) .* exp(x(4).*xdata)','x','xdata');
% double exponential fit (5 par)
dblexpBfun = inline('x(1) + x(2) .* exp(x(3).*xdata) + x(4) .* exp(x(5).*xdata)','x','xdata');
% triple exponential fit
triexpfun = inline('x(1) .* exp(x(2).*xdata) + x(3) .* exp(x(4).*xdata) + x(5) .* exp(x(6).*xdata)','x','xdata');
% power law (2 par)
pwrAfun = inline('x(1) .* xdata.^x(2)','x','xdata');
% power law (3 par)
pwrBfun = inline('x(1) + x(2) .* xdata.^x(3)','x','xdata');
% inverse 2nd order fit
inv2fun = inline('x(1) + x(2)./xdata + x(3)./xdata.^2','x','xdata');
% Rayleigh function
ralfun = inline('x(1) .*xdata .* exp(-x(2).* sqrt(xdata))','x','xdata');
% Gamma function
gamdistfun = inline('1/(x(2).^x(1) .* gamma(x(1))) .* xdata.^(x(1)-1) .* exp(-xdata ./ x(2)) .* x(3)','x','xdata');
% piecewise linear fit 
piecelinfun = inline('(x(1) + x(2) .* xdata) .* (xdata<=x(3)) + (x(4) + x(5) .* xdata) .* (xdata > x(3)) ','x','xdata');
% truncated power law
trunpwrfun = inline('x(1) .* xdata.^x(2) .* exp(x(3) .* xdata)','x','xdata');
% Pareto 
Paretofun = inline('1 - 1 ./ xdata.^x(1)','x','xdata');
% modified power law I
modpwrAfun = inline('x(1) .* (1 - xdata.^-x(2))','x','xdata');
% modified power law II
modpwrBfun = inline('x(1) .* (1 + xdata.^x(2))','x','xdata');
% piecewise power law
piecepwrfun = inline('(x(1) .* xdata .^ x(2)) .* (xdata<=x(3)) + (x(4) .* xdata .^ x(5)) .* (xdata > x(3)) ','x','xdata');
% inverse cubic
invcubicfun = inline('x(1) ./ xdata.^3 + x(2) ./ xdata.^2 + x(3) ./ xdata + x(4)','x','xdata');
% quadratic
quadfun = inline('x(1) + x(2) .* xdata + x(3).*xdata.^2','x','xdata');
% exponential rise-to-max (2par) a*(1-exp(-b*x))
expmaxAfun = inline('x(1) .* (1-exp(-x(2).*xdata))','x','xdata');
% exponential rise-to-max (3par) a + b*(1-exp(-c*x))
expmaxBfun = inline('x(1) + x(2).* (1-exp(-x(3).*xdata))','x','xdata');
% dbl exponential max (4par)  a*(1-exp(-b*x)) + c*(1-exp(-d*x))
dblexpmaxAfun = inline('x(1) .* (1-exp(-x(2).*xdata))+x(3) .* (1-exp(-x(4).*xdata))','x','xdata');
% inverse first-order
invfun = inline('x(1) + x(2)./xdata','x','xdata');
% stretched exponential decay (3 par)
strexpfun = inline('x(1).*exp(-(xdata./x(2)).^x(3))','x','xdata');
% cubic
cubicfun =  inline('x(1) + x(2) .* xdata + x(3).*xdata.^2 + x(4).*xdata.^3','x','xdata');
% quartic 
quarticfun =  inline('x(1) + x(2) .* xdata + x(3).*xdata.^2 + x(4).*xdata.^3 + x(5).*xdata.^4','x','xdata');
% hyperbola, double (4 par)
hyprdblAfun =  inline('x(1)*xdata ./ (x(2) + xdata) + x(3).*xdata ./(x(4)+ xdata)','x','xdata');
% hyperbola, double (5 par)
hyprdblBfun =  inline('x(1)*xdata ./ (x(2) + xdata) + x(3).*xdata ./(x(4)+ xdata) + x(5) .* xdata','x','xdata');
% quintic
quinticfun =  inline('x(1) + x(2) .* xdata + x(3).*xdata.^2 + x(4).*xdata.^3 + x(5).*xdata.^4 + x(6).*xdata.^5','x','xdata');
% exp and linear (4 par)
explinfun = inline('x(1) + x(2) .* exp(x(3).*xdata) + x(4).*xdata','x','xdata');
% rational (4 par)
rtnl4fun = inline('(x(1) + x(2) .* xdata) ./ (1 + x(3).*xdata + x(4).*xdata.^2)','x','xdata');
% rational (5 par)
rtnl5fun = inline('(x(1) + x(2) .* xdata + x(3).*xdata.^2) ./ (1 + x(4).*xdata + x(5).*xdata.^2)','x','xdata');
% Gaussian (3 par)
Gaussfun = inline('x(1).*exp(-(xdata-x(2)).^2./x(3).^2)','x','xdata');
% double Gaussian (6 par)
dblGaussfun = inline('x(1).*exp(-(xdata-x(2)).^2./x(3).^2) + x(4).*exp(-(xdata-x(5)).^2./x(6).^2)','x','xdata');
% Gaussian (3 par)
logGaussfun = inline('x(1).*exp(-(log(xdata)-x(2)).^2./x(3).^2)','x','xdata');


nfits = length(fits);
SS = zeros(nfits,1);
coeffs = cell(nfits,1);
residuals = cell(nfits,1);

if isempty(LBC) LBC = cell(nfits,1); end
if isempty(UBC) UBC = cell(nfits,1); end
if isempty(IVC) IVC = cell(nfits,1); end

for loop = 1:nfits
    if ~isempty(LBC{loop})
        LB = LBC{loop};
        UB = UBC{loop};
    else
        LB = [];
        UB = [];
    end
    if ~isempty(IVC{loop}) IV = IVC{loop}; end

    switch fits(loop)
        case 1
            if isempty(IVC{loop}) IV = rand(2,1); end
            
            [coeffs{loop},SS(loop),residuals{loop}] = lsqcurvefit(linfun,IV, xdata, ydata,LB,UB,options);   
        case 2
            if isempty(IVC{loop}) IV = rand(2,1); end

            [coeffs{loop},SS(loop),residuals{loop}] = lsqcurvefit(expAfun,IV, xdata, ydata,LB,UB,options);   
        case 3
            if isempty(IVC{loop}) IV = rand(3,1); end
            [coeffs{loop},SS(loop),residuals{loop}] = lsqcurvefit(expBfun,IV, xdata, ydata,LB,UB,options);

        case 4 
            if isempty(IVC{loop}) IV = rand(4,1); end
            [coeffs{loop},SS(loop),residuals{loop}] = lsqcurvefit(dblexpAfun,IV, xdata, ydata,LB,UB,options);   

        case 5
            if isempty(IVC{loop}) IV = rand(5,1); end
            [coeffs{loop},SS(loop),residuals{loop}] = lsqcurvefit(dblexpBfun,IV, xdata, ydata,LB,UB,options);   
        case 6
            if isempty(IVC{loop}) IV = rand(6,1); end
            [coeffs{loop},SS(loop),residuals{loop}] = lsqcurvefit(triexpfun,IV, xdata, ydata,LB,UB,options);   

        case 7
            if isempty(IVC{loop}) IV = rand(2,1); end
            [coeffs{loop},SS(loop),residuals{loop}] = lsqcurvefit(pwrAfun,IV, xdata, ydata,LB,UB,options);   

        case 8
            if isempty(IVC{loop}) IV = rand(3,1); end
            [coeffs{loop},SS(loop),residuals{loop}] = lsqcurvefit(pwrBfun,IV, xdata, ydata,LB,UB,options);   
        case 9
            if isempty(IVC{loop}) IV = rand(3,1); end
            [coeffs{loop},SS(loop),residuals{loop}] = lsqcurvefit(inv2fun,IV, xdata, ydata,LB,UB,options);   
        case 10
            if isempty(IVC{loop}) IV = rand(2,1); end
            [coeffs{loop},SS(loop),residuals{loop}] = lsqcurvefit(ralfun,IV, xdata, ydata,LB,UB,options);    
        case 11            
            if isempty(IVC{loop}) IV = rand(3,1); end
            [coeffs{loop},SS(loop),residuals{loop}] = lsqcurvefit(gamdistfun,IV, xdata, ydata,LB,UB,options);          
        case 12
            if isempty(IVC{loop}) IV = rand(5,1); end
            [coeffs{loop},SS(loop),residuals{loop}] = lsqcurvefit(piecelinfun,IV, xdata, ydata,LB,UB,options);
%             max_ix = find(ydata == max(ydata));
%             x1 = xdata(1:max_ix); y1 = ydata(1:max_ix);
%             x2 = xdata(max_ix:end); y2 = ydata(max_ix:end);
%             [coeffs1,SS1,residuals1] = lsqcurvefit(linfun,IV(1:2), x1, y1,LB,UB);          
%             [coeffs2,SS2,residuals2] = lsqcurvefit(linfun,IV(3:4), x2, y2,LB,UB);          
%             coeffs{loop} = [coeffs1 coeffs2];
%             SS(loop) = SS1 + SS2;
%             residuals{loop} = [residuals1; residuals2];
        case 13
            if isempty(IVC{loop}) IV = rand(3,1); end
            [coeffs{loop},SS(loop),residuals{loop}] = lsqcurvefit(trunpwrfun,IV, xdata, ydata,LB,UB,options);              
        case 14
            if isempty(IVC{loop}) IV = rand(1,1); end
            [coeffs{loop},SS(loop),residuals{loop}] = lsqcurvefit(Paretofun,IV, xdata, ydata,LB,UB,options);              
        case 15
            if isempty(IVC{loop}) IV = rand(2,1); end
            [coeffs{loop},SS(loop),residuals{loop}] = lsqcurvefit(modpwrAfun,IV, xdata, ydata,LB,UB,options);              
        case 16
            if isempty(IVC{loop}) IV = rand(2,1); end
            [coeffs{loop},SS(loop),residuals{loop}] = lsqcurvefit(modpwrBfun,IV, xdata, ydata,LB,UB,options);              
        case 17
            if isempty(IVC{loop}) IV = rand(5,1); end
            [coeffs{loop},SS(loop),residuals{loop}] = lsqcurvefit(piecepwrfun,IV, xdata, ydata,LB,UB,options);              
        case 18
            if isempty(IVC{loop}) IV = rand(4,1); end
            [coeffs{loop},SS(loop),residuals{loop}] = lsqcurvefit(invcubicfun,IV, xdata, ydata,LB,UB,options);              
        case 19
            if isempty(IVC{loop}) IV = rand(3,1); end
            [coeffs{loop},SS(loop),residuals{loop}] = lsqcurvefit(quadfun,IV, xdata, ydata,LB,UB,options);              
        case 20
            if isempty(IVC{loop}) IV = rand(2,1); end
            [coeffs{loop},SS(loop),residuals{loop}] = lsqcurvefit(expmaxAfun,IV, xdata, ydata,LB,UB,options);              
        case 21
            if isempty(IVC{loop}) IV = rand(3,1); end
            [coeffs{loop},SS(loop),residuals{loop}] = lsqcurvefit(expmaxBfun,IV, xdata, ydata,LB,UB,options);              
        case 22
            if isempty(IVC{loop}) IV = rand(4,1); end
            [coeffs{loop},SS(loop),residuals{loop}] = lsqcurvefit(dblexpmaxAfun,IV, xdata, ydata,LB,UB,options);              
        case 23
            if isempty(IVC{loop}) IV = rand(2,1); end
            [coeffs{loop},SS(loop),residuals{loop}] = lsqcurvefit(invfun,IV, xdata, ydata,LB,UB,options);                  
        case 24
            if isempty(IVC{loop}) IV = rand(3,1); end
            [coeffs{loop},SS(loop),residuals{loop}] = lsqcurvefit(strexpfun,IV, xdata, ydata,LB,UB,options);                      
        case 25
            if isempty(IVC{loop}) IV = rand(4,1); end
            [coeffs{loop},SS(loop),residuals{loop}] = lsqcurvefit(cubicfun,IV, xdata, ydata,LB,UB,options);                      
        case 26
            if isempty(IVC{loop}) IV = rand(5,1); end
            [coeffs{loop},SS(loop),residuals{loop}] = lsqcurvefit(quarticfun,IV, xdata, ydata,LB,UB,options);                      
        case 27
            if isempty(IVC{loop}) IV = rand(4,1); end
            [coeffs{loop},SS(loop),residuals{loop}] = lsqcurvefit(hyprdblAfun,IV, xdata, ydata,LB,UB,options);                      
        case 28
            if isempty(IVC{loop}) IV = rand(5,1); end
            [coeffs{loop},SS(loop),residuals{loop}] = lsqcurvefit(hyprdblBfun,IV, xdata, ydata,LB,UB,options);                      
        case 29
            if isempty(IVC{loop}) IV = rand(6,1); end
            [coeffs{loop},SS(loop),residuals{loop}] = lsqcurvefit(quinticfun,IV, xdata, ydata,LB,UB,options);                      
        case 30
            if isempty(IVC{loop}) IV = rand(4,1); end
            [coeffs{loop},SS(loop),residuals{loop}] = lsqcurvefit(explinfun,IV, xdata, ydata,LB,UB,options);                      
        case 31
            if isempty(IVC{loop}) IV = rand(4,1); end
            [coeffs{loop},SS(loop),residuals{loop}] = lsqcurvefit(rtnl4fun,IV, xdata, ydata,LB,UB,options);                      
        case 32
            if isempty(IVC{loop}) IV = rand(5,1); end
            [coeffs{loop},SS(loop),residuals{loop}] = lsqcurvefit(rtnl5fun,IV, xdata, ydata,LB,UB,options);                      
        case 33
            if isempty(IVC{loop}) IV = rand(3,1); end
            [coeffs{loop},SS(loop),residuals{loop}] = lsqcurvefit(Gaussfun,IV, xdata, ydata,LB,UB,options);                      
        case 34
            if isempty(IVC{loop}) IV = rand(6,1); end
            [coeffs{loop},SS(loop),residuals{loop}] = lsqcurvefit(dblGaussfun,IV, xdata, ydata,LB,UB,options);                      
        case 35
            if isempty(IVC{loop}) IV = rand(3,1); end
            [coeffs{loop},SS(loop),residuals{loop}] = lsqcurvefit(logGaussfun,IV, xdata, ydata,LB,UB,options);                      
   
    end
end

AICs = zeros(1,nfits);
for j = 1:nfits
    AICs(j) = AIC(SS(j),length(ydata),length(coeffs{j}));
end

switch rtn 
    case {'AIC'}
    %%% then find smallest AIC
    [AICsort,I] = sort(AICs);
    case {'SS'}
    %% then find smallest sum-of-squares
    [SSsort,I] = sort(SS);
    otherwise
        error('returned function criterion unknown')
end

best_model = fits(I(1));        
best_coeffs = coeffs{I(1)};
best_SS = SS(I(1));
best_residuals = residuals{I(1)};

switch best_model
    case 1
        best_fun = linfun;
    case 2
        best_fun = expAfun;
    case 3
        best_fun = expBfun;
    case 4
        best_fun = dblexpAfun;
    case 5
        best_fun = dblexpBfun;      
    case 6
        best_fun = triexpfun;
    case 7
        best_fun = pwrAfun;
    case 8
        best_fun = pwrBfun;
    case 9
        best_fun = inv2fun;
    case 10
        best_fun = ralfun;
    case 11
        best_fun = gamdistfun;    
    case 12
        best_fun = piecelinfun;    
    case 13
        best_fun = trunpwrfun;
    case 14
        best_fun = Paretofun;
    case 15
        best_fun = modpwrAfun;
    case 16
        best_fun = modpwrBfun;
    case 17
        best_fun = piecepwrfun;   
     case 18
        best_fun = invcubicfun;   
    case 19
        best_fun = quadfun;
    case 20
        best_fun = expmaxAfun;
    case 21
        best_fun = expmaxBfun;
    case 22
        best_fun = dblexpmaxAfun;       
    case 23
        best_fun = invfun;
    case 24
        best_fun = strexpfun;
    case 25
        best_fun = cubicfun;      
    case 26
        best_fun = quarticfun;   
    case 27
        best_fun = hyprdblAfun;   
    case 28
        best_fun = hyprdblBfun;      
    case 29
        best_fun = quinticfun;    
    case 30
        best_fun = explinfun;           
    case 31
        best_fun = rtnl4fun;
    case 32
        best_fun = rtnl5fun;       
    case 33
        best_fun = Gaussfun;
    case 34
        best_fun = dblGaussfun;
    case 35
        best_fun = logGaussfun;
                
end

