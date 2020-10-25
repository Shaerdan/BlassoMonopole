function [P] = EvaluateP0(X,I,Measurement,Mesh,FL2Norm)

%   Compute P = Measurement - Estimation.

%   X --- estimated source locations vector.
%   I --- estimated source intensity vector.
%   Estimation --- computed analytically from:
%                    Sum_i(I(i)*PhiComponent(X(i)));

    I_Normalized = I./FL2Norm;
    PhiComponent = ComputePotentialComponent(X,Mesh);
    Estimation=PhiComponent*I_Normalized';
    P = -Estimation + Measurement;
%     figure(40)
%     surf(Estimation);


end

