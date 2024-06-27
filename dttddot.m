function z = dttddot(Fttd,Gttd)
%DTTDDOT returns the dot product of two TTD tensors
%   z = dttddot(Fttd,Gttd) computes the dot product 'z' of tensor 'Fttd'
%   and 'Gttd' that are given by their tensor train decompositions. The dot
%   product of two d-dimensional tensors is defined as z =
%   \sum_{i_1=1}^{n_1}...\sum_{i_d=1}^{n_d}
%   Fttd(i_1,i_2,...,i_d)*Gttd(i_1,i_2,...,i_d)
%
%   Example
%   Fttd{1}(1,:,:) = [1;3];
%   Fttd{2}(:,:,1) = [5 8];
%   Gttd{1}(1,:,:) = [9;1];
%   Gttd{2}(:,:,1) = [8 3];
%   dttddot(Fttd,Gttd)
%   ans =
%       768
%
%   See also DTTDDISTFROB
%
%   Copyright 2023 IDM (Ivo Punèocháø)

% Get sizes of input args
sizFttd = size(Fttd);
sizGttd = size(Gttd);

% Check dimensions
if ~iscell(Fttd) || numel(sizFttd)~=2 || sizFttd(1)~=1 || size(Fttd{1},1)~=1 || size(Fttd{end},3)~=1
    exception = MException('MyToolbox:dttddot:wrongInArg',...
        'The first input argument must be a tensor in TTD given as a row cell array. Also the first and last core tensors must have outer ranks equal to one.');
    throw(exception);
end

if ~iscell(Gttd) || numel(sizGttd)~=2 || sizGttd(1)~=1 || size(Gttd{1},1)~=1 || size(Gttd{end},3)~=1
    exception = MException('MyToolbox:dttddot:wrongInArg',...
        'The second input argument must be a tensor in TTD given as a row cell array. Also the first and last core tensors must have outer ranks equal to one.');
    throw(exception);
end

if sizFttd(2)~=sizGttd(2)
    exception = MException('MyToolbox:dttddot:wrongInArg',...
        'The dimensions of the both tensors must be the same.');
    throw(exception);
end


% The original code by Petr Tichavsky was slightly adapted to remove
% special treatment of the first and last term of the tensor train

% Initialize the value of the dot product
z = 1;

% Perform left to right sweep through the tensor
for j = 1:sizFttd(2)
    % Get sizes of the i-th tensor cores
    [rfp, nf, rf] = size(Fttd{j});
    [rgp, ng, rg] = size(Gttd{j});
    
    % Check that the sizes of tensors at the i-th dimension are the same
    if nf~=ng
        exception = MException('MyToolbox:dttddot:wrongInArg',...
            'The sizes of the tensors differ in the %i-th dimension.');
        throw(exception);
    end
    
    % Permute and reshape tensor cores to get r_{i-1}*r_{i} by n_{i} matrices
    Fttdm = reshape(permute(Fttd{j},[1,3,2]),rfp*rf, nf);
    Gttdm = reshape(permute(Gttd{j},[1,3,2]),rgp*rg, nf);
    
    % Perform product of elements and summation over i-th dimension
    % (Xttdm*Yttdm'). Permute and reshape the result to a rx_{i-1}*ry_{i-1}
    % by rx_{i}*ry_{i} matrix
    tmp = reshape(permute(reshape(Fttdm*Gttdm',rfp,rf,rgp,rg),[1,3,2,4]),rfp*rgp,rf*rg);

    % Perform summation over ranks during by the left to right sweep
    z = z*tmp;
end

% The following code is a naive implementation that uses several FOR cycles
%
% % Compute the first dimesion separately
% z = zeros(size(Xttd{1},3),size(Yttd{1},3));
% for alpha1x = 1:size(Xttd{1},3)
%     for alpha1y = 1:size(Yttd{1},3)
%         z(alpha1x,alpha1y) = reshape(Xttd{1}(:,:,alpha1x),1,size(Xttd{1},2))*reshape(Yttd{1}(:,:,alpha1y),size(Yttd{1},2),1);
%     end
% end
%
% % Computation proceeds with the rest of dimensions
% for d = 2:sizXttd(2)
%     newz = zeros(size(Xttd{d},3),size(Yttd{d},3));
%     for alphax = 1:size(Xttd{d},3)
%         for alphay = 1:size(Yttd{d},3)
%             for alphaxp = 1:size(Xttd{d-1},3)
%                 for alphayp = 1:size(Yttd{d-1},3)
%                     tmp = 0;
%                     for i = 1:size(Xttd{d},2)
%                         tmp = tmp + Xttd{d}(alphaxp,i,alphax)*Yttd{d}(alphayp,i,alphay);
%                     end
%                     newz(alphax,alphay) = newz(alphax,alphay) + z(alphaxp,alphayp)*tmp;
%                 end
%             end
%         end
%     end
%     z = newz;
% end


end