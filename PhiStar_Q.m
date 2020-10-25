function [PhiSQ] = PhiStar_Q(Mesh,Q,SwitchNorm)

%   Given Q(x), compute Phi^*(Q)
%%
%   PhiSQ -- Phi^*(Q)
%   F  --- The integral kernel, ref Equation (45), Page 6 of the report.
%   FL2Norm --- L2Norm of F;
%%
dx = Mesh.ThetaQLine;
dy = Mesh.PsiQLine;
PhiSQ = zeros(Mesh.LengthRadiusS,Mesh.LengthThetaS,Mesh.LengthPsiS);
for i=1:Mesh.LengthRadiusS
    for j=1:Mesh.LengthThetaS
        for k=1:Mesh.LengthPsiS
            CosGamma = cos(Mesh.ThetaS(j)).*cos(Mesh.ThetaQ)+sin(Mesh.ThetaS(j)).*...
                sin(Mesh.ThetaQ).*cos(Mesh.PsiS(k)- Mesh.PsiQ);
            DistanceSQ  = sqrt(1 + Mesh.RadiusS(i).^2 - 2* Mesh.RadiusS(i).*CosGamma);
            F      = (2./DistanceSQ) - log(1-Mesh.RadiusS(i).*CosGamma + DistanceSQ);
            %             if (SwitchNorm == 1)
            %                 Integrand0 = (F.^2).*sin(Mesh.ThetaQ);
            %                 FL2Norm = sqrt(trapz(dy, trapz(dx,Integrand0,2)));
            %             else
            %                 FL2Norm = 1;
            %             end
            %             F_Normalized  = F./FL2Norm;
            %             Integrand1 = F_Normalized.*Q.*sin(Mesh.ThetaQ);
            Integrand1 = F.*Q.*sin(Mesh.ThetaQ);
            PhiSQ(i,j,k) = trapz(dx, trapz(dy,Integrand1,2));
        end
    end
end
end

