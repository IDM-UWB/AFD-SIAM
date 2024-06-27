function [Vrhs,uopt] = getvp(subidx,V,inputs,eta,mL_h,mV_h)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
nu = size(inputs,2);
Q = zeros(nu,1);
for ui = 1:nu
    Q(ui) = mL_h(subidx,inputs(:,ui)) + eta*mV_h(subidx,inputs(:,ui),V);
end
[Vrhs,idx] = min(Q,[],1);

uopt = inputs(:,idx);

end
