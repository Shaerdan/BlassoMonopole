function [J,GradientPseudo] = ObjectiveFuncNonLinearLasso(X,Measurement,Lambda, Mesh,FL2Norm)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
SourceNum = length(X)/4;
IntensityNormalized = X(1:SourceNum)./FL2Norm;
Locations = X(SourceNum +1:end);

[PhiComponent,CosGamma,DistPQ] = ComputePotentialComponent(Locations,Mesh);

Estimation=sum(bsxfun(@times,PhiComponent,reshape(IntensityNormalized,1,1,SourceNum)),3);
Descrepency = Estimation - Measurement;
LeastSquareTerm = norm(Descrepency,2);
RegularizationTerm = norm(IntensityNormalized,1);

J = 0.5*LeastSquareTerm^2 + Lambda*RegularizationTerm;
[GradientPseudo] = GradPseudo(X,Descrepency,PhiComponent,CosGamma,DistPQ,Mesh,FL2Norm, Lambda);
% grad_s

end

