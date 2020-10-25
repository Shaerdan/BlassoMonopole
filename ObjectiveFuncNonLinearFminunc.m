function [J,GradientPseudo] = ObjectiveFuncNonLinearFminunc(X,Measurement,Lambda, Mesh,SwitchNorm)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

SourceNum = length(X)/4;
Locations = X(SourceNum +1:end);
I = X(1:SourceNum);
[PhiComponent] = ComputePotentialComponent(Locations,Mesh,SwitchNorm);

Estimation=PhiComponent*I;
Descrepency = Estimation - Measurement;
LeastSquareTerm = norm(Descrepency,2);
% RegularizationTerm = norm(X(1:SourceNum),1);
RegularizationTerm = norm(I,1);
J = 0.5*LeastSquareTerm^2 + Lambda*RegularizationTerm;
[GradientPseudo] = GradPseudoLBFGSB(X,Measurement,Mesh,Lambda,SwitchNorm);
% [GradientPseudo] = GradPseudoFminunc(X,Measurement,Mesh,Lambda,SwitchNorm);

end

