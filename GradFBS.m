function Grad = GradFBS(I,M,b)

%%	Compute the gradient of ||Estimates - Measurement||_2^2;
%   I --- Source intensity,
%   Estimates --- Estimation of the data,
%   Estimates = Sum_i( IntensityNormalized(i)*PhiComponent(i)).

%% Initialize:

Grad = M*I - b;


end

