function [XHat] = FBSSolver(X,M,b,Measurement,FBS,PhiComponent,Mesh,Lambda,FL2Norm)
%   FBS iteration:  X_{k+1/2} = X_{k} - Tau*Gradient(||Estimate - ...
%   Measurement||_2^2);  X_{k+1} = BackProjection(Tau,Lambda,X_{k+1/2});
%   X   --- The unknown parameters;
%


% Initialization:
XHat = X;
SourceNum = length(X);
GradCheck = 1;
IterationCount =1;
DiscrepancyResidualCheck = zeros(1,FBS.Iteration);
StepSize = 1;
dx = Mesh.ThetaQLine';
dy = Mesh.PsiQLine';
SIN1 = sin(Mesh.ThetaQ);
% FBS loop (termination condition is: Gradient <= tolerence && Stepsize <= tolerence && InterationCount = Maxiteration):
while (GradCheck > FBS.GradTol && IterationCount < FBS.Iteration && StepSize > FBS.StepTol)
    Estimates = sum(bsxfun(@times,PhiComponent,reshape(XHat./FL2Norm,1,1,SourceNum)),3);
    Discrepancy = Estimates - Measurement; 
    DiscrepancyGradient = GradFBS(XHat,M,b);
    IntegrandLeastSquare = (Discrepancy.^2)*SIN1;
    DiscrepancyResidualCheck(IterationCount+1) = trapz(dy,trapz(dx,IntegrandLeastSquare,2))+...
        Lambda*norm(XHat,1);    
    StepSize = abs(DiscrepancyResidualCheck(IterationCount+1) - DiscrepancyResidualCheck(IterationCount));
    
    XHat = XHat - FBS.Tau*DiscrepancyGradient;
    
    SoftThresholding = sign(XHat).*max(0, abs(XHat) - FBS.Tau*Lambda);
    
    XHat = SoftThresholding;
    
    DiscrepancyGradientCheck = norm(DiscrepancyGradient,2);
    
    if (mod(IterationCount,FBS.DisplayFrequency)==0)
        fprintf(' lfbs log10(Discrepency) = %d \n',log10(DiscrepancyResidualCheck(IterationCount)));
        fprintf(' lfbs log10(DiscrepancyGradientient) = %d \n',log10(DiscrepancyGradientCheck));
        fprintf(' lfbs Iter = %d \n',IterationCount);
    end
    
    IterationCount = IterationCount +1;
end

end

