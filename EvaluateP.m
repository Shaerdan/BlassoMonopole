function [P] = EvaluateP(X,I,Measurement,Mesh,FL2Norm)

%   Compute P = Measurement - Estimation.

%   X --- estimated source locations vector.
%   I --- estimated source intensity vector.
%   PhiEstimated --- computed analytically from:
%                    Sum_i(I(i)*PhiComponent(X(i)));

SourceNum = length(I);

if (~isempty(I))
    I_Normalized = I./FL2Norm;
    PhiComponent = ComputePotentialComponent(X,Mesh);
    PhiEstimated=sum(bsxfun(@times,PhiComponent,reshape(I_Normalized,1,1,SourceNum)),3);
    P = -PhiEstimated + Measurement;
%     figure(40)
%     surf(PhiEstimated);
else
    P = Measurement;
end

end

