function [J] = GradTest(X,Measurement,Mesh,Lambda)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
SIN1 = sin(Mesh.ThetaQ);
dx = Mesh.ThetaQLine;
dy = Mesh.PsiQLine;

SourceNum = length(X)/4;
Locations = X(SourceNum +1:end);
FL2Norm = L2NormF(Locations,Mesh); 

[PhiComponent] = ComputePotentialComponent(Locations,Mesh);

Estimation=sum(bsxfun(@times,PhiComponent,reshape(X(1:SourceNum)./FL2Norm',1,1,SourceNum)),3);
Descrepency = Estimation - Measurement;
IntegrandLeastSquare = (Descrepency.^2).*SIN1;
LeastSquareTerm = trapz(dy,trapz(dx,IntegrandLeastSquare,2));
RegularizationTerm = norm(X(1:SourceNum),1);

J = 0.5*LeastSquareTerm + Lambda*RegularizationTerm;

end