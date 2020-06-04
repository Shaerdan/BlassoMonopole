function [ZNew] = RelaxedFBSSolver(X,Measurement,FBS,PhiComponent,Mesh,Lambda,FL2Norm)
%   FBS iteration:  X_{k+1/2} = X_{k} - FBS.Tau*Gradient(||Estimate - ...
%   Measurement||_2^2);  X_{k+1} = BackProjection(FBS.Tau,Lambda,X_{k+1/2});
%   X   --- The unknown parameters;
%   


% Initialization:
ZOld = X;
GradCheck = 1;
IterationCount =1;
DescrepencyResidualCheck = zeros(1,FBS.Iteration);
StepSize = 1;
% FBS loop (termination condition is: Gradient <= tolerence && Stepsize <= tolerence && InterationCount = Maxiteration):
while (GradCheck > FBS.GradTol && IterationCount < FBS.Iteration && StepSize > FBS.StepTol)
    
    [DescrepencyGradient,Descrepency] = GradFBS(ZOld,Measurement,Mesh,PhiComponent,FL2Norm);
        
    DescrepencyResidualCheck(IterationCount+1) = norm(Descrepency,2);
    
    StepSize = abs(DescrepencyResidualCheck(IterationCount+1) - DescrepencyResidualCheck(IterationCount));
    
    ZNew = ZOld - FBS.Tau*DescrepencyGradient;
    ZNew(ZNew == 0) = 1e-14;
    
    SoftThresholding = ZNew.*max(0, 1 - FBS.Tau*Lambda./(abs(ZNew)));
    
    ZNew = SoftThresholding;
    ZNew = ZNew +FBS.Mu*(ZNew - ZOld);
    ZOld = ZNew; 
    DescrepencyGradientCheck = norm(DescrepencyGradient,2);
    
    if (mod(IterationCount,FBS.DisplayFrequency)==0)
        fprintf(' lfbs log10(Discrepency) = %d \n',log10(DescrepencyResidualCheck(IterationCount)));
        fprintf(' lfbs log10(DescrepencyGradientient) = %d \n',log10(DescrepencyGradientCheck));
        fprintf(' lfbs FBS.Iteration = %d \n',IterationCount);
    end
    
    IterationCount = IterationCount +1;
end


end

