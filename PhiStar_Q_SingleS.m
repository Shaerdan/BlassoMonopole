function [PhiSQ] = PhiStar_Q_SingleS(Mesh,S,Q,SwitchNorm)

%   Given Q(x), compute Phi^*(Q)
%%
%   PhiSQ -- Phi^*(Q)
%   F  --- The integral kernel, ref Equation (45), Page 6 of the report.
%   FL2Norm --- L2Norm of F;
%%
dx = Mesh.ThetaQLine;
dy = Mesh.PsiQLine;
PhiSQ = zeros(Mesh.LengthRadiusS,Mesh.LengthThetaS,Mesh.LengthPsiS);

            CosGamma = cos(S.ThetaS).*cos(Mesh.ThetaQ)+sin(S.ThetaS).*...
                sin(Mesh.ThetaQ).*cos(S.PsiS- Mesh.PsiQ);
            DistanceSQ  = sqrt(1 + S.RadiusS.^2 - 2* S.RadiusS.*CosGamma);
            F      = - log(1-S.RadiusS.*CosGamma + DistanceSQ);
            %             if (SwitchNorm == 1)
            %                 Integrand0 = (F.^2).*sin(Mesh.ThetaQ);
            %                 FL2Norm = sqrt(trapz(dy, trapz(dx,Integrand0,2)));
            %             else
            %                 FL2Norm = 1;
            %             end
            %             F_Normalized  = F./FL2Norm;
            %             Integrand1 = F_Normalized.*Q.*sin(Mesh.ThetaQ);
            Integrand1 = F.*Q.*sin(Mesh.ThetaQ);
%             surf(Integrand1)
            PhiSQ = trapz(dy', trapz(dx',Integrand1,2));

end

