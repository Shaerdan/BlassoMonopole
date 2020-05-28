function [XHat] = FBSSolver(X,Measurement,Iter,PhiComponent,theta,psi,GradTolerence,StepSizeTolerence,...
                      Tau,Lambda,DisplayFreq,F_L2Norm)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
% x   --- initial value of the unknown parameters


% Initialization:
XHat = X;
GradCheck = 1;
IterationCount =1;
DescrepencyResidualCheck = zeros(1,Iter);
StepSize = 1;
% FBS loop (termination condition is: Gradient <= tolerence && Stepsize <= tolerence && InterationCount = Maxiteration):
while (GradCheck > GradTolerence && IterationCount < Iter && StepSize > StepSizeTolerence)
    
    [DescrepencyGradient,] = GradFBS(XHat',Measurement,theta,psi,PhiComponent,F_L2Norm,Lambda);
    
    DescrepencyGradient = DescrepencyGradient';
    
    DescrepencyResidualCheck(IterationCount+1) = norm(D,2);
    
    StepSize = abs(DescrepencyResidualCheck(IterationCount+1) - DescrepencyResidualCheck(IterationCount));
    
    XHat = XHat - Tau*DescrepencyGradient;
    XHat(XHat == 0) = 1e-14;
    
    prox_g = XHat.*max(0, 1-Tau*Lambda./(abs(XHat)));
    
    XHat = prox_g;
    
    DescrepencyGradientCheck = norm(DescrepencyGradient,2);
    
    if (mod(IterationCount,DisplayFreq)==0)
        fprintf(' lfbs log10(discrepency) = %d \n',log10(DescrepencyResidualCheck(IterationCount)));
        fprintf(' lfbs log10(DescrepencyGradientient) = %d \n',log10(DescrepencyGradientCheck));
        fprintf(' lfbs Iter = %d \n',IterationCount);
    end
    
    IterationCount = IterationCount +1;
end


end

