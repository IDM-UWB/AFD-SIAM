function [x,P,pmu] = hyperstate2sufstat(xi,nx,nModels,idxvec2mat,reducedProbabilityVector)
%HYPERSTATE2SUFSTAT Summary of this function goes here
%   Detailed explanation goes here

sizXi = size(xi);

if numel(sizXi)~=2 || sizXi(2)~=1
    exception = MException('MyToolbox:hyperstate2sufstat:wtongInArgs',...
        'The function ''%s'' requires that the first input argument is a columns vector.',mfilename');
    throw(exception);
end

nxix = nx*nModels;
nxiP = nModels*nx*(nx+1)/2;

x = reshape(xi(1:nxix),nx,nModels);
P = xi(nxix+idxvec2mat);
pmu = xi(nxix+nxiP+1:end);
if reducedProbabilityVector
    pmu(end+1) = 1 - sum(pmu);
end


end