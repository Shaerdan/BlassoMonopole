function [J,GradientPseudo] = ObjectiveFuncNonLinearLBFGSB(X,Measurement,Lambda, Mesh)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
SourceNum = length(X)/4;
Locations = X(SourceNum +1:end);
FL2Norm = L2NormF(Locations,Mesh);
IntensityNormalized = X(1:SourceNum)./FL2Norm';

[PhiComponent,CosGamma,DistPQ] = ComputePotentialComponent(Locations,Mesh);

Estimation=sum(bsxfun(@times,PhiComponent,reshape(IntensityNormalized,1,1,SourceNum)),3);
Descrepency = Estimation - Measurement;
LeastSquareTerm = norm(Descrepency,2);
RegularizationTerm = norm(IntensityNormalized,1);

J = 0.5*LeastSquareTerm^2 + Lambda*RegularizationTerm;
GradientPseudo = GradPseudoLBFGSB(X,Descrepency,PhiComponent,CosGamma,DistPQ,Mesh,FL2Norm, Lambda);
GradientPseudo = GradientPseudo';
end

