function [Eta] = ObjEtaMatrix(X,Lambda,Mesh,Solution,Data,FL2Norm,SwitchNorm)
%%   Compute the objective function -(1/2)*Eta;
%   Eta -- ref Equation (45), Page 6 of the report:
%          Eta = Eta(Mesh.S) = integral (
%          F_Normalized(Mesh.Q,Mesh.S)*P(Mesh.Q)dMesh.Q ),
%          dMesh.Q = sin(Mesh.ThetaQ) * dMesh.ThetaQ * dMesh.PsiQ.
%   F  --- The integral kernel, ref Equation (45), Page 6 of the report.
%   FL2Norm --- L2Norm of F;
%%

dx = Mesh.ThetaQLine;
dy = Mesh.PsiQLine;
SIN1 = sin(Mesh.ThetaQ);
Radius = [Solution.Radius, X(1)];
Theta = [Solution.Theta, X(2)];
Psi = [Solution.Psi, X(3)];
Location = [Radius, Theta, Psi];
Intensity = [Solution.Intensity, 1];
P   = EvaluateP(Location, Intensity,Data.Measurement,Mesh,FL2Norm);

CosGamma = cos(X(2)).*cos(Mesh.ThetaQ)+sin(X(2)).*...
    sin(Mesh.ThetaQ).*cos(X(3)- Mesh.PsiQ);
DistanceSQ  = sqrt(1 + X(1)^2 - 2* X(1)*CosGamma);
Phi      = (2./DistanceSQ) - log(1-X(1)*CosGamma + DistanceSQ);

figure(22)
 K = convhull(100);
 patch('vertices',points,'faces',K,'cdata',P); % draw it
 axis equal tight vis3d % set axis
 view(3) 

if (SwitchNorm == 1)
    Phi = Phi/norm(Phi,2);
end
PhiStarP = trapz(dy,trapz(dx,Phi.*P.*SIN1,2));
Eta = abs(PhiStarP)/Lambda;
Eta = -0.5*Eta^2;
% Eta = -0.5*norm(Eta,inf)^2;

end

