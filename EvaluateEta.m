function [Eta] = EvaluateEta(Lambda,Mesh,P,SwitchNorm)

%   Evaluate Eta over a descrete volumetric mesh:

%   Eta -- ref Equation (45), Page 6 of the report:
%          Eta = Eta(Mesh.S) = integral (
%          F_Normalized(Mesh.Q,Mesh.S)*P(Mesh.Q)dMesh.Q ),
%          dMesh.Q = sin(Mesh.ThetaQ) * dMesh.ThetaQ * dMesh.PsiQ.
%   F  --- The integral kernel, ref Equation (45), Page 6 of the report.
%   FL2Norm --- L2Norm of F;
dx = Mesh.ThetaQLine;
dy = Mesh.PsiQLine;
Eta = zeros(Mesh.LengthRadiusS,Mesh.LengthThetaS,Mesh.LengthPsiS);
for i=1:Mesh.LengthRadiusS
    for j=1:Mesh.LengthThetaS
        for k=1:Mesh.LengthPsiS
            CosGamma = cos(Mesh.ThetaS(j)).*cos(Mesh.ThetaQ)+sin(Mesh.ThetaS(j)).*...
                sin(Mesh.ThetaQ).*cos(Mesh.PsiS(k)- Mesh.PsiQ);
            DistanceSQ  = sqrt(1 + Mesh.RadiusS(i).^2 - 2* Mesh.RadiusS(i).*CosGamma);
            F      = (2./DistanceSQ) - log(1-Mesh.RadiusS(i).*CosGamma + DistanceSQ);
            if (SwitchNorm == 1)
                Integrand0 = (F.^2).*sin(Mesh.ThetaQ);
                FL2Norm = sqrt(trapz(dy, trapz(dx,Integrand0,2)));
            else
                FL2Norm = 1;
            end
            F_Normalized  = F./FL2Norm;
            Integrand1 = F_Normalized.*P.*sin(Mesh.ThetaQ);
            Eta(i,j,k) = trapz(dy, trapz(dx,Integrand1,2));
        end
    end
end

Eta = abs(Eta)/Lambda;
end

