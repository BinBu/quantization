function pd = mypdffit(data,nbins,dist)
%MYPDFFIT Fit probability distribution to data.
%   MYPDFFIT(DATA,NBINS) plots a pdf of the values in the vector DATA,
%   along with a normal density function with parameters estimated from the
%   data.  NBINS is the number of bars in the histogram. With one input
%   argument, NBINS is set to the square root of the number of elements in
%   DATA. 
%
%   MYPDFFIT(DATA,NBINS,DIST) plots a pdf of the values in the vector DATA,
%   along with a density from the DIST distribution. DIST can take the following values:
%
%         'beta'                             Beta
%         'birnbaumsaunders'                 Birnbaum-Saunders
%         'exponential'                      Exponential
%         'extreme value' or 'ev'            Extreme value
%         'gamma'                            Gamma
%         'generalized extreme value' 'gev'  Generalized extreme value
%         'generalized pareto' or 'gp'       Generalized Pareto (threshold 0)
%         'inverse gaussian'                 Inverse Gaussian
%         'logistic'                         Logistic
%         'loglogistic'                      Log logistic
%         'lognormal'                        Lognormal
%         'negative binomial' or 'nbin'      Negative binomial
%         'nakagami'                         Nakagami
%         'normal'                           Normal
%         'poisson'                          Poisson
%         'rayleigh'                         Rayleigh
%         'rician'                           Rician
%         'tlocationscale'                   t location-scale
%         'weibull' or 'wbl'                 Weibull
%
%   PD = MYPDFFIT(...) returns an object PD representing
%   the fitted distribution. PD is an object in a class derived from the
%   ProbDist class.

%   Written by Bin Hu in MATLAB R2009a, learned from histfit.m (Copyright 1993-2008 The MathWorks, Inc.) in MATLAB R2009a. 
%   $Revision: 1.0.0.0 $  $Date: 2017/07/18 $

% validation verification and format conversion.
if ~isvector(data)
    error('stats:histfit:VectorRequired','DATA must be a vector.');
end

data = data(:);
data(isnan(data)) = [];
n = numel(data);

if nargin<2 || isempty(nbins)
    nbins = ceil(sqrt(n));
elseif ~isscalar(nbins) || ~isnumeric(nbins) || ~isfinite(nbins) ...
        || nbins~=round(nbins)
    error('stats:histfit:BadNumBins','NBINS must be a positive integer.')
end

% Do histogram calculations
[bincounts,bincenters]=hist(data,nbins);

% Plot the pdf of the values in the vector data
binwidth = bincenters(2) - bincenters(1); % Finds the width of each bin
y_data = bincounts / (n * binwidth); % Normalize the density
plot(bincenters,y_data,'b-','LineWidth',2);

% Fit distribution to data
if nargin<3 || isempty(dist)
    dist = 'normal';
end
try
    pd = fitdist(data,dist);
catch myException
    if isequal(myException.identifier,'stats:ProbDistUnivParam:fit:NRequired')
        % Binomial is not allowed because we have no N parameter
        error('stats:histfit:BadDistribution',...
            'Binomial distribution not allowed.')
    else
        % Pass along another other errors
        throw(myException)
    end
end

% Find range for plotting
q = icdf(pd,[0.0013499 0.99865]); % three-sigma range for normal distribution
x = linspace(q(1),q(2));
if ~pd.Support.iscontinuous
    % For discrete distribution use only integers
    x = round(x);
    x(diff(x)==0) = [];
end

% Overlay the density
y_fitted = pdf(pd,x);
np = get(gca,'NextPlot');
set(gca,'NextPlot','add')
plot(x,y_fitted,'r-','LineWidth',2);
set(gca,'NextPlot',np)
xmin1 = min(bincenters);
xmax1 = max(bincenters);
xmin2 = min(x);
xmax2 = max(x);
xlim([max(xmin1, xmin2), min(xmax1, xmax2)]);