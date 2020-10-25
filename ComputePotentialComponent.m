function [PhiComponent,CosGamma,DistancePQ] = ComputePotentialComponent(X,Mesh,SwitchNorm)

%   Compute the analytical value of PhiComponent:

% Initialization:
global LogSwitch

SourceNum = length(X)/3;
LocationRadius = X(1:SourceNum);
LocationTheta = X(SourceNum+1:2*SourceNum);
LocationPsi = X(2*SourceNum+1:3*SourceNum);


% Compute PhiComponent(ThetaQ,PsiQ,i):
for i = 1:SourceNum
    CosGamma(:,i) = cos(LocationTheta(i))*cos(Mesh.ThetaQ(:))+...
        sin(LocationTheta(i))*sin(Mesh.ThetaQ(:)).*cos(Mesh.PsiQ(:)-LocationPsi(i));
    DistancePQ(:,i) = sqrt(1+LocationRadius(i)^2-2*LocationRadius(i)* CosGamma(:,i))  ;
    PhiComponent(:,i) = ((2./DistancePQ(:,i)) - LogSwitch*log(1-LocationRadius(i)*CosGamma(:,i) + DistancePQ(:,i) ));
    if (SwitchNorm == 1)
        PhiComponent(:,i) = PhiComponent(:,i)/norm(PhiComponent(:,i));
    end
end

end

