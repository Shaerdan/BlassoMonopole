function [XHat] = FBSTest(X,M,b,Measurement,FBS,A,Lambda)
%   FBS iteration:  X_{k+1/2} = X_{k} - Tau*Gradient(||Estimate - ...
%   Measurement||_2^2);  X_{k+1} = BackProjection(Tau,Lambda,X_{k+1/2});
%   X   --- The unknown parameters;
%


% Initialization:
XHat = X;
GradCheck = 1;
IterationCount =1;
DiscrepancyResidualCheck = zeros(1,FBS.Iteration);
StepSize = 1;
% FBS loop (termination condition is: Gradient <= tolerence && Stepsize <= tolerence && InterationCount = Maxiteration):
while (GradCheck > FBS.GradTol && IterationCount < FBS.Iteration && StepSize > FBS.StepTol)
    Estimates = A*XHat;
    Discrepancy = Estimates - Measurement;
    DiscrepancyResidualCheck(IterationCount+1) = 0.5*norm(Discrepancy,2)^2+...
        Lambda*norm(XHat,1);
    StepSize = abs(DiscrepancyResidualCheck(IterationCount+1) - DiscrepancyResidualCheck(IterationCount));
    DiscrepancyGradient = GradFBS(XHat,M,b);
    GradCheck = norm(DiscrepancyGradient,2);
    XHat = XHat - FBS.Tau*DiscrepancyGradient;
    
    SoftThresholding = sign(XHat).*max(0, abs(XHat) - FBS.Tau*Lambda);
    
    XHat = SoftThresholding;
        
    if (mod(IterationCount,FBS.DisplayFrequency)==0)
        fprintf(' lfbs log10(Discrepency) = %d \n',log10(DiscrepancyResidualCheck(IterationCount)));
        fprintf(' lfbs log10(DiscrepancyGradientient) = %d \n',log10(GradCheck));
        fprintf(' lfbs Iter = %d \n',IterationCount);
    end
    
    IterationCount = IterationCount +1;
end

end

