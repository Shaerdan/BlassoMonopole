function [D] = CalcDiscrepancy(X,Measurement,PhiComponent,FL2Norm,Mesh)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
SourceNum = length(X);
IntensityNormalized = X./FL2Norm';
Estimates = sum(bsxfun(@times,PhiComponent,reshape(IntensityNormalized,1,1,SourceNum)),3);
dx = Mesh.ThetaQ(1,:)';
dy = Mesh.PsiQ(:,1);
D = Estimates - Measurement;
end

