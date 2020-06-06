function [XHat] = FBSSolver(X,Measurement,FBS,PhiComponent,Mesh,Lambda,FL2Norm)
%   FBS iteration:  X_{k+1/2} = X_{k} - Tau*Gradient(||Estimate - ...
%   Measurement||_2^2);  X_{k+1} = BackProjection(Tau,Lambda,X_{k+1/2});
%   X   --- The unknown parameters;
%   


% Initialization:
XHat = X;
GradCheck = 1;
IterationCount =1;
DescrepencyResidualCheck = zeros(1,FBS.Iteration);
StepSize = 1;
% FBS loop (termination condition is: Gradient <= tolerence && Stepsize <= tolerence && InterationCount = Maxiteration):
while (GradCheck > FBS.GradTol && IterationCount < FBS.Iteration && StepSize > FBS.StepTol)
    
    [DescrepencyGradient,Descrepency] = GradFBS(XHat,Measurement,Mesh,PhiComponent,FL2Norm);
        
    DescrepencyResidualCheck(IterationCount+1) = norm(Descrepency,2);
    
    StepSize = abs(DescrepencyResidualCheck(IterationCount+1) - DescrepencyResidualCheck(IterationCount));
    
    XHat = XHat - FBS.Tau*DescrepencyGradient;
    
    SoftThresholding = sign(XHat).*max(0, abs(XHat) - Lambda);
    
    XHat = SoftThresholding;
    
    DescrepencyGradientCheck = norm(DescrepencyGradient,2);
    
    if (mod(IterationCount,FBS.DisplayFrequency)==0)
        fprintf(' lfbs log10(Discrepency) = %d \n',log10(DescrepencyResidualCheck(IterationCount)));
        fprintf(' lfbs log10(DescrepencyGradientient) = %d \n',log10(DescrepencyGradientCheck));
        fprintf(' lfbs Iter = %d \n',IterationCount);
    end
    
    IterationCount = IterationCount +1;
end


end

