function [P] = EvaluateP(X,Intensity,Measurement,Mesh,FL2Norm)

%   Compute P = Measurement - Estimation.

%   X --- estimated source locations vector.
%   Intensity --- estimated source intensity vector.
%   PhiEstimated --- computed analytically from:
%                    Sum_i(Intensity(i)*PhiComponent(X(i)));

SourceNum = length(Intensity);

if (~isempty(Intensity))
    Intensity_Normalized = Intensity./FL2Norm';
    PhiComponent = ComputePotentialComponent(X,Mesh);
    PhiEstimated=sum(bsxfun(@times,PhiComponent,reshape(Intensity_Normalized,1,1,SourceNum)),3);
    P = -PhiEstimated + Measurement;
    figure(40)
    surf(PhiEstimated);
else
    P = Measurement;
end

end

