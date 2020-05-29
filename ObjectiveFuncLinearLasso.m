function [ObjectiveF,Grad] = ObjectiveFuncLinearLasso(X,Measurement,PhiComponent,Mesh,F_L2Norm,Lambda)
%UNTITLED6 Summary of this function goes here
%   Compute objective function and the gradient of the descrepency term.
%   X   --- The unknown parameters

[Grad,D] = GradFBS(X,Measurement,Mesh,PhiComponent,F_L2Norm);

DescrepencyTerm = norm(D,2);
RegularizationTerm = Lambda*norm(X,1);

ObjectiveF = 0.5*DescrepencyTerm^2 + RegularizationTerm;

end


