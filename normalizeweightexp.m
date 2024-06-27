function w = normalizeweightexp(e,varargin)
%NORMALIZEWEIGHTEXP Normalize weights that include exponentials
%   w = NORMALIZEWEIGHTEXP(e,c,method) performs numerically robust
%   normalization of weights c.*exp(e) such that their row sum equals one,
%   i.e. w = c.*exp(e)./sum(c.*exp(e),2). If no coefficients 'c' are given
%   they are set to one. The implemented methods are:
%       1) (default) The greatest exponent 'e' that has nonzero coefficient
%       'c' is found and canceled out of the fraction. Then one of the
%       summands of denominator has value of 'c' that corresponds to that
%       maximum exponent.
%
%       2) The greatest product c.*exp(e) is found and
%       canceled out of the fraction. Then one of the summands of
%       denominator has value 1. This second option is computationally more
%       demanding but can hypotetically improve numerical robustness when
%       the greatest exponent 'e' has nonzero coefficient 'c' that is
%       really small compared to other coefficients. It seems that this
%       second methods is rarely needed and the first one is enough to get
%       numerically robust normalization.


% Check the maximum number of optional input arguments
numVarArgMax = 2;
if length(varargin) > numVarArgMax
    exception = MException('MyToolbox:normalizeweightexp:tooManyOptInArgs',...
        'The function ''%s'' accepts at most %i optional input arguments.',mfilename,numVarArgMax);
    throw(exception);
end

% Set default values for optional input arguments
optArg = {1,1};

% Flag all nonempty input arguments in varargin
flagNotEmpty = cellfun(@isnotempty,varargin);

% Overwrite default values by the ones specified in varargin
optArg(flagNotEmpty) = varargin(flagNotEmpty);

% Place optional input arguments in memorable variable names
[c,method] = optArg{:};

% Get dimensions of exponents and coeficients
sizC = size(c);
sizE = size(e);

% Check input arguments
if numel(sizE)~=2 || any(e(:)>0)
    exception = MException('MyToolbox:normalizeweightexp:wrongInArgs',...
        'The first argument must be a vector of nonpositive exponents');
    throw(exception);
end

isSameSize = isequal(sizC,sizE);
if numel(sizC)~=2 || any(c(:)<0) || ~(isSameSize || isscalar(c)) || (isscalar(c) && c<=0) ||  (isSameSize && any(sum(c,2)==0))
    exception = MException('MyToolbox:normalizeweightexp:wrongInArgs',...
        'The second argument must be a positive scalar or conforming matrix of nonnegative weights, each row has at least one positive weight');
    throw(exception);
end

switch method
    case 1
        % Find maximum exponent for each row
        if isscalar(c)
            maxE = max(e,[],2);
        else
            %TODO - do an implementation without FOR LOOP
            maxE  = zeros(sizC(1),1);
            for i = 1:sizC(1)
                maxE(i) = max(e(i,c(i,:)>0));
            end
        end
        
        % Compute non-normalized weights
        tmp = exp(e-maxE);
        % Check for possile infs, if they occur at positions where c is
        % zero it is ok an set them to zero
        flagInf = isinf(tmp(:));
        flagZero = c(:)==0;
        if all(flagZero | ~flagInf)
            tmp(flagInf) = 0;
        end
        
        w = c.*tmp;
        
    case 2
        error('not implemented for matrix case')
        if  sizC(1)==1
            % Computation of ratios if c and e are row vectors
            ratios = (c'./c).*exp(e'-e); % This works only in new versions of MATLAB
        else
            % Computation of ratios if c and e are column vectors
            ratios = (c./c').*exp(e-e'); % This works only in new versions of MATLAB
        end
        
        % flag column with where all ratios are less or equal to one
        idxCol = all(ratios<=1);
        
        % Put selected terms into auxiliary variables, to treat the case
        % when there is not the unique maximum
        maxE = e(idxCol);
        maxC = c(idxCol);
        
        % Compute non-normalized weights
        w = (c/maxC(1)).*exp(e-maxE(1)); % Use the first maximum found
        
        
    otherwise
        error('unknown method')
end

% Common part of the code for both methods

% Compute the sum of the weights
sumW = sum(w,2);

if any(sumW<=0)
    exception = MException('MyToolbox:normalizeweightexp:wrongInArgs',...
        'Numerical problem while normalizing weights');
    throw(exception);
end

% Normalize weights
w = w./sumW;


end