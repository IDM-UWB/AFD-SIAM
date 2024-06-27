legendText = {};
figure
load bellman21
gridElements = prod(cellfun(@length,thetaGridLvl));
plot([1:10],dVfrobnorm_ci/gridElements,'b')
hold on
plot([1:10],dVfrobnorm_li/gridElements,'b*')
legendText = [legendText,sprintf('constant_li %1.2g elements',gridElements)];
legendText = [legendText,sprintf('linear_li %1.2g elements',gridElements)];

load bellman31
gridElements = prod(cellfun(@length,thetaGridLvl));
plot([1:10],dVfrobnorm_ci/gridElements,'g')
plot([1:10],dVfrobnorm_li/gridElements,'g*')
legendText = [legendText,sprintf('constant_li %1.2g elements',gridElements)];
legendText = [legendText,sprintf('linear_li %1.2g elements',gridElements)];

gridElements = 6.5757e+11;
dVfrobnorm41_ci = [3.542603e-07,3.803603e-08,2.622642e-08,1.759744e-08,1.785800e-08,1.577347e-08,3.243087e-07,2.304252e-07];
dVfrobnorm41_li = [3.542603e-07,3.756256e-08,1.895792e-08,8.126750e-09,3.853163e-09,3.422740e-09,1.889424e-09,1.289244e-09];
plot([1:8],dVfrobnorm41_ci,'r')
plot([1:8],dVfrobnorm41_li,'r*')
legendText = [legendText,sprintf('constant_li %1.2g elements',gridElements)];
legendText = [legendText,sprintf('linear_li %1.2g elements',gridElements)];
legend(legendText,'Interpreter','latex')


load bellman21_ci
gridElements = prod(cellfun(@length,thetaGridLvl));
plot([1:10],dVfrobnorm_ci/gridElements,'r')
hold on
plot([1:10],dVfrobnorm_li/gridElements,'r*')
legendText = [legendText,sprintf('constant_ci %1.2g elements',gridElements)];
legendText = [legendText,sprintf('linear_ci %1.2g elements',gridElements)];

