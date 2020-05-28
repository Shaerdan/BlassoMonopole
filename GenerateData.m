function [Data] = GenerateData(S,X,Mesh,level)

%   Generate data from known(imposed) sources.
%   Data.Measurement is generated from analytical solution:  
%   Data.Measurement = Sum_i( Intensity(i)*PhiComponent(i) ); 
%   EXPLAIN NORMALIZATION HERE.
%   ComputePotentialComponent(): return an analytical computation of
%   PhiComponent.
%   L2NormF(X,Mesh): return the normalization of the forward integral
%   kernel, using the source locations X.

ns = length(S);
PhiComponent = ComputePotentialComponent(X,Mesh);
FL2Norm = L2NormF(X,Mesh);
S_norm = S./FL2Norm';
Data.Measurement = sum(bsxfun(@times,PhiComponent,reshape(S_norm,1,1,ns)),3);
Data.Noise = level*randn(size(Data.Measurement));
Data.Measurement = Data.Measurement + Data.Noise;
figure(1)
surf(Data.Measurement)

end

