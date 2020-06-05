function [Grad,Descrepency] = GradFBS(I,Measurement,Mesh,PhiComponent,FL2Norm)

%%	Compute the gradient of ||Estimates - Measurement||_2^2;
%   I --- Source intensity,
%   Estimates --- Estimation of the data,
%   Estimates = Sum_i( IntensityNormalized(i)*PhiComponent(i)).

%% Initialize:
SourceNum = length(I);
IntensityNormalized = I./FL2Norm;
Estimates = sum(bsxfun(@times,PhiComponent,reshape(IntensityNormalized,1,1,SourceNum)),3);
dx = Mesh.ThetaLine';
dy = Mesh.PsiLine';
Descrepency = Estimates - Measurement;
Grad = zeros(1,SourceNum);

%% Compute d(||Estimates - Measurement||_2^2)/d(Intensity):
for i = 1:SourceNum
   Diff_Estimates_Intensity = PhiComponent(:,:,i)/FL2Norm(i);
   Integrand = Descrepency.*Diff_Estimates_Intensity.*sin(Mesh.ThetaQ);
   Grad(i) = trapz(dy,trapz(dx,Integrand,2));
end

end

