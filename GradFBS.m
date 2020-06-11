function Grad = GradFBS(Discrepancy,PhiComponent,Mesh,FL2Norm, SourceNum)

%%	Compute the gradient of ||Estimates - Measurement||_2^2;
%   I --- Source intensity,
%   Estimates --- Estimation of the data,
%   Estimates = Sum_i( IntensityNormalized(i)*PhiComponent(i)).

%% Initialize:
dx = Mesh.ThetaQLine';
dy = Mesh.PsiQLine';
SIN1 = sin(Mesh.ThetaQ);
Grad = zeros(1,SourceNum);

%% Compute d(0.5*||Estimates - Measurement||_2^2)/d(Intensity):
for i = 1:SourceNum
   Diff_Estimates_Intensity = PhiComponent(:,:,i)/FL2Norm(i);
   Integrand = Discrepancy.*Diff_Estimates_Intensity.*SIN1;
   Grad(i) = trapz(dy,trapz(dx,Integrand,2));
end

end

