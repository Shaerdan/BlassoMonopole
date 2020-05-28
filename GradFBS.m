function [Grad,Descrepency] = GradFBS(Intensity,Measurement,Mesh,PhiComponent,F_L2Norm)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
SourceNum = length(Intensity);
IntensityNormalized = Intensity./F_L2Norm;
Estimates = sum(bsxfun(@times,PhiComponent,reshape(IntensityNormalized,1,1,SourceNum)),3);
dx = Mesh.ThetaQ(1,:)';
dy = Mesh.PsiQ(:,1);
Descrepency = Estimates - Measurement;
Grad = zeros(1,SourceNum);

for i = 1:SourceNum
   dEdS = PhiComponent(:,:,i)/F_L2Norm(i);
   Integrand = Descrepency.*dEdS.*sin(Mesh.ThetaQ);
   Grad(i) = trapz(dy,trapz(dx,Integrand,2));
end

end

