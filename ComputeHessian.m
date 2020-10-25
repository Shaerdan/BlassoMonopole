function [M,b] = ComputeHessian(Measurement,PhiComponent,Conditioning)
%UNTITLED2 Summary of this function goes here
%   M:   Hessian matrix of f_des, M = integral (AQ'*AQ) dQ;
%   b:   Measurement component of the gradient of f_des;
%   b =  integral(AQ*Measurement) dQ;
% SourceNum = length(FL2Norm);
[M,SourceNum] = size(PhiComponent);
AQ = zeros(M,SourceNum);

%% Compute d(0.5*||Estimates - Measurement||_2^2)/d(Intensity):


for i = 1:SourceNum
    AQ(:,i) = PhiComponent(:,i);
end
if (Conditioning == 1)
    ConditionTerm = 0.5*(norm(Measurement,2))^2;
else
    ConditionTerm = 1;
end
M = ConditionTerm*(AQ'*AQ);
b = ConditionTerm*(AQ'*Measurement);
end

