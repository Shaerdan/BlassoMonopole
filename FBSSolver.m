function [XHat] = FBSSolver(X,Measurement,Iter,PhiComponent,Mesh,GradTolerence,StepSizeTolerence,...
                      Tau,Lambda,DisplayFreq,F_L2Norm)
%   FBS iteration:  X_{k+1/2} = X_{k} - Tau*Gradient(||Estimate - ...
%   Measurement||_2^2);  X_{k+1} = BackProjection(Tau,Lambda,X_{k+1/2});
%   X   --- The unknown parameters;
%   


% Initialization:
XHat = X;
GradCheck = 1;
IterationCount =1;
DescrepencyResidualCheck = zeros(1,Iter);
StepSize = 1;
% FBS loop (termination condition is: Gradient <= tolerence && Stepsize <= tolerence && InterationCount = Maxiteration):
while (GradCheck > GradTolerence && IterationCount < Iter && StepSize > StepSizeTolerence)
    
    [DescrepencyGradient,Descrepency] = GradFBS(XHat,Measurement,Mesh,PhiComponent,F_L2Norm);
        
    DescrepencyResidualCheck(IterationCount+1) = norm(Descrepency,2);
    
    StepSize = abs(DescrepencyResidualCheck(IterationCount+1) - DescrepencyResidualCheck(IterationCount));
    
    XHat = XHat - Tau*DescrepencyGradient;
    XHat(XHat == 0) = 1e-14;
    
    BackProjection = XHat.*max(0, 1-Tau*Lambda./(abs(XHat)));
    
    XHat = BackProjection;
    
    DescrepencyGradientCheck = norm(DescrepencyGradient,2);
    
    if (mod(IterationCount,DisplayFreq)==0)
        fprintf(' lfbs log10(Discrepency) = %d \n',log10(DescrepencyResidualCheck(IterationCount)));
        fprintf(' lfbs log10(DescrepencyGradientient) = %d \n',log10(DescrepencyGradientCheck));
        fprintf(' lfbs Iter = %d \n',IterationCount);
    end
    
    IterationCount = IterationCount +1;
end


end

