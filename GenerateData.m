function [Data] = GenerateData(Intensity,X,Mesh,level,SwitchNorm)

%   Generate data from known(imposed) sources.
%   Data.Measurement is generated from analytical solution:
%   Data.Measurement = Sum_i( Intensity(i)*PhiComponent(i) );
%   EXPLAIN NORMALIZATION HERE.
%   ComputePotentialComponent(): return an analytical computation of
%   PhiComponent.
%   L2NormF(X,Mesh): return the normalization of the forward integral
%   kernel, using the source locations X.

NumSource = length(Intensity);
PhiComponent = ComputePotentialComponent(X,Mesh,SwitchNorm);
Data.Measurement = PhiComponent*Intensity';
Data.Noise = level*randn(size(Data.Measurement));
Data.Measurement = Data.Measurement + Data.Noise;
% [n1,n2] = size(Mesh.ThetaQ);
figure(1)
surf(Mesh.ThetaQ,Mesh.PsiQ,reshape(Data.Measurement,size(Mesh.ThetaQ)));
xlabel('\theta')
ylabel('\psi')
zlabel('\phi^{d}')
legend('Exact Data')
end

