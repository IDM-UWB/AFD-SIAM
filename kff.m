function [xf,Pxxf,yp,Pyyp,K,e] = kff(y,xp,Pxxp,C,Sigmav,meanv)
%KFF Filtering step of the Kalman filter
%   [xf,Pxxf,yp,Pyyp,K,e] = kff(y,xp,Pxxp,C,Sigmav,meanv) computes
%   filtering mean value 'xf', filtering covariance matrix 'Pxxf',
%   predictive mean value 'yp', predictive covariance matrix 'Pyyp', Kalman
%   gain 'K', and innovation 'e' based on the measurement 'y', predictive
%   mean value 'xp', predictive covariance matrix 'Pxxp', measurement gain
%   matrix 'C', noise covariance matrix 'Sigmav', and noise mean value
%   'meanv'

% Get the dimension of output and state
[ny,nx] = size(C);

% If mean value of measurement noise is not given, zero is assumed
if nargin<6
    meanv = zeros(ny,1);
end

% Compute the predictive mean value and covariance matrix of output
yp = C*xp + meanv;
Pyyp = C*Pxxp*C' + Sigmav;

% Compute the Kalman gain
K = Pxxp*C'/Pyyp;

% Compute innovation
e = y - yp;

% Compute the filtering mean value and covariance matrix of state
xf = xp + K*e;
tmp = eye(nx) - K*C;
Pxxf = tmp*Pxxp*tmp' + K*Sigmav*K';

% Basic form - computationally efficient but suffers from numerical
% problems such as loss of symmetry or positive semidefinitiness
% Pxxf = (eye(nx)-K*C)*Pxxp;

% Symmetrical form - keeps symmetry but still prone to numerical problems -
% I've seen an example of t-invariant model where gain and covariance
% matrices are almost in steady state and suddenly start oscillating with
% increasing amplitude
% Pxxf = Pxxp - K*Pyyp*K';

% Joseph form - keeps symmetry, positive semidefinnitness and is valid for
% non-optimal gain
% tmp = eye(nx) - K*C;
% Pxxf = tmp*Pxxp*tmp' + K*Sigmav*K';


end