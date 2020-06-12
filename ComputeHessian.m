function [M,b] = ComputeHessian(Measurement,PhiComponent,FL2Norm,Mesh)
%UNTITLED2 Summary of this function goes here
%   M:   Hessian matrix of f_des, M = integral (AQ'*AQ) dQ;
%   b:   Measurement component of the gradient of f_des; 
%   b =  integral(AQ*Measurement) dQ;
dx = Mesh.ThetaQLine';
dy = Mesh.PsiQLine';
SIN1 = sin(Mesh.ThetaQ);
SourceNum = length(FL2Norm);
AQ = zeros(size(PhiComponent));
M = zeros(SourceNum,SourceNum);
b = zeros(SourceNum,1);
%% Compute d(0.5*||Estimates - Measurement||_2^2)/d(Intensity):
for i = 1:SourceNum
   AQ(:,:,i) = PhiComponent(:,:,i)/FL2Norm(i);
end
for i = 1:SourceNum
    for j = 1:SourceNum
    M(i,j) = trapz(dy,trapz(dx,AQ(:,:,i).*AQ(:,:,j).*SIN1,2));
    end
    b(i) = trapz(dy,trapz(dx,AQ(:,:,i).*Measurement.*SIN1,2));
end
end

