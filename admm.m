function [Xkp1] = admm(X0,M,b,Lambda,ADMM)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Steps = 1;
[m,n] = size(M);
Ykp1 = X0;
Zkp1 = X0;
ADMM.Rho = 1;
GradCheck = 1;
while GradCheck>ADMM.Tol && Steps < ADMM.Steps
    Yk = Ykp1; Zk = Zkp1;
    Xkp1 = (inv(M+ADMM.Rho*eye(m,n)))*(b+ADMM.Rho*Zk - Yk);
    ThreshArg = Xkp1 + Yk/ADMM.Rho;
    SoftThresh = sign(ThreshArg).*max(0, abs(ThreshArg) - Lambda/ADMM.Rho);
    Zkp1 = SoftThresh;
    Ykp1 = Yk + ADMM.Rho*(Xkp1 - Zkp1);
    Steps = Steps + 1;
    GradCheck = norm(M*Xkp1 - b,2);
    fprintf(' admm log10(Gradient) = %d \n',log10(GradCheck));
    fprintf(' admm Iter = %d \n',Steps);
end


end

