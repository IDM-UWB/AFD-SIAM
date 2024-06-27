function xi = sufstat2hyperstate(x,P,pmu,idxtril2vec,reducedProbabilityVector)
%SUFSTAT2HYPERSTATE Summary of this function goes here
%   Detailed explanation goes here

if reducedProbabilityVector
    xi = [x(:); squeeze(P(idxtril2vec)); pmu(1:end-1)];
else
    xi = [x(:); squeeze(P(idxtril2vec)); pmu];
end


end