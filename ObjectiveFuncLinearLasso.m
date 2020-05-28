function [f,Grad] = ObjectiveFuncLinearLasso(x,Measurement,PhiComponent,Mesh,F_L2Norm,Lambda)
%UNTITLED6 Summary of this function goes here
%   Compute objective function and the gradient of the descrepency term.
%   x   --- initial value of the unknown parameters

[Grad,D] = GradFBS(x,Measurement,Mesh,PhiComponent,F_L2Norm);
Grad = Grad';

DescrepencyTerm = norm(D,2);
RegularizationTerm = Lambda*norm(x,1);

f = 0.5*DescrepencyTerm^2 + RegularizationTerm;

end


