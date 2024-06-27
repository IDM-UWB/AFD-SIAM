function [xp,Pxxp] = kfp(u,xf,Pxxf,A,B,SigmaGw,meanGw)
%KFP Predictive step of the Kalman filter
%   [xp,Pxxp] = kfp(u,xf,Pxxf,A,B,Sigmaw,meanw) computes predictive mean
%   value 'xp' and predictive covariance matrix 'Pxxp', based on the input
%   'u', filtering mean value 'xf', filtering covariance matrix 'Pxxf',
%   state matrix 'A', input matrix 'B', noise covariance matrix 'Sigmaw',
%   and noise mean value 'meanw'

% If mean value of state noise is not given, zero is assumed
if nargin<7
    nx = size(xf,1);
    meanGw = zeros(nx,1);
end

isemptyU = isempty(u);
isemptyB = isempty(B);
if (isemptyB && ~isemptyU) || (~isemptyB && isemptyU)
    error('Input matrix B and input u have to be both either empty or nonempty')
end

if isemptyU
    xp = A*xf + meanGw;
else
    xp = A*xf + B*u + meanGw;
end
Pxxp = A*Pxxf*A' + SigmaGw;


end