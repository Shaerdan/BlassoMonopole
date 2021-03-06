function [P] = EvaluateP(X,I,Measurement,Mesh,SwitchNorm)

%   Compute P = Measurement - Estimation.

%   X --- estimated source locations vector.
%   I --- estimated source intensity vector.
%   Estimation --- computed analytically from:
%                    Sum_i(I(i)*PhiComponent(X(i)));

SourceNum = length(I);

if (~isempty(I))
    PhiComponent = ComputePotentialComponent(X,Mesh,SwitchNorm);
    Estimation = PhiComponent*I';
    P = -Estimation + Measurement;
%     figure(40)
%     surf(Estimation);
else
    P = Measurement;
end

end

