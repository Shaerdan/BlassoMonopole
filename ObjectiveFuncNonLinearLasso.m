function [J,grad_s] = ObjectiveFuncNonLinearLasso(X,Measurement,Lambda, Mesh,F_L2Norm)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
SourceNum = length(X)/4;
Intensity = X(1:SourceNum)./F_L2Norm;
Locations = X(SourceNum +1:end);

[PhiComponent,CosGamma,L_in] = ComputePotentialComponent(Locations,Mesh);
% for i=1:n1
%     PhiComponent(:,:,i) = PhiComponent(:,:,i)/F_L2Norm(i);
% end


Estimation=sum(bsxfun(@times,PhiComponent,reshape(Intensity,1,1,SourceNum)),3);

Descrepency = Estimation - Measurement;
LeastSquareTerm = norm(Descrepency,2);
RegularizationTerm = norm(Intensity,1);

J = 0.5*LeastSquareTerm^2 + Lambda*RegularizationTerm;
[grad_s] = GradPseudo(X,Descrepency,PhiComponent,CosGamma,L_in,Mesh,F_L2Norm, Lambda);
% grad_s

end

