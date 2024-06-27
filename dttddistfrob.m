function distFrob = dttddistfrob(Fttd,Gttd)
%DTTDDISTFROB returns the Frobenius norm of difference of two TTD tensors
%   distFrob = dttddistfrob(Fttd,Gttd) returns the Frobenius norm
%   'distFrob' of difference between tensors 'Fttd' and 'Gttd' that are
%   specified by their tensor train decompositions
%
%   Example
%   Fttd{1}(1,:,:) = [1;3];
%   Fttd{2}(:,:,1) = [5 8];
%   Gttd{1}(1,:,:) = [9;1];
%   Gttd{2}(:,:,1) = [8 3];
%   distFrob = dttddistfrob(Fttd,Gttd)
%   ans =
%       73.0753036257804
%
%   See also DTTDDOT, DTTDNORMFROB
%
%   Copyright 2023 IDM (Ivo Punèocháø)

% No check of input arguments is performed here because the code is a
% simple and all checks are done by the dot product function anyway

% Compute the square of Frobenius norm using the dot products
% |X-Y|_{F}^{2} = <X,X> - 2*<X,Y> + <Y,Y>
distFrob2 = dttddot(Fttd,Fttd)-2*dttddot(Fttd,Gttd)+dttddot(Gttd,Gttd);

% If numerical errors make the result negative, set it to zero
if distFrob2<0
    distFrob2 = 0;
    warning('Numerical issues make the norm negative')
end

% Compute the Frobenius norm
distFrob = sqrt(distFrob2);


end