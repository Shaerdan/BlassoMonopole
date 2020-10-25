function [GradPseudo,Jac] = GradPseudoLBFGSB(X,Measurement,Mesh,Lambda,SwitchNorm,ConditionTerm)

%%  Compute pseudo gradient of 0.5*||Estimates(X) - Measurement||_2^2 + ...
%%  ... Lambda* ||X||_1;
%  GradPseudo = d/dX (0.5*||Estimates - Measurement||_2^2) + ...
%  ... SubGradient of Lambda*||X||_1;
%  I --- Source intensity;

%% Initialize:
global LogSwitch
SourceNum = length(X)/4;
I = X(1:SourceNum);
Locations = X(SourceNum+1:end);
SourceRadius = X(SourceNum +1:2*SourceNum);
SourceTheta = X(2*SourceNum +1:3*SourceNum);
SourcePsi = X(3*SourceNum +1:4*SourceNum);
SIN1 = sin(Mesh.ThetaQ(:));
COS1 = cos(Mesh.ThetaQ(:));
GradPseudo = zeros(1,SourceNum);
GradSub = GradPseudo;
Jac  = zeros(length(SIN1),4*SourceNum); PDotQ = Jac; Di_DistPQ = Jac;
Di_PDotQ = Jac; Di_PhiPi = Jac; Di_CosGamma = Jac;
[PhiPi,CosGamma,DistPQ] = ComputePotentialComponent(Locations,Mesh,SwitchNorm);
FL2Norm = L2NormF(Locations,Mesh);
Estimation=PhiPi*I;
Descrepency = Estimation - Measurement;
%% Compute the Pseudo gradient of the objective function:
PhiPiNorm = FL2Norm(1:SourceNum);
%% dObjectiveFunction/dIntensity
GradSub(1:SourceNum) = Lambda*sign(I);
% GradSub(1:SourceNum) = 0;
GradSub(SourceNum+1:4*SourceNum) = 0;
for i = 1:SourceNum
    Jac(:,i)  = PhiPi(:,i);
    %% dObjectiveFunction/dRadius
    PDotQ(:,i) = SourceRadius(i)*CosGamma(:,i);
    Di_DistPQ(:,i) = (SourceRadius(i)-CosGamma(:,i))./DistPQ(:,i);
    Di_PDotQ(:,i)  = CosGamma(:,i);
    Di_PhiPi(:,i)  = (-2./(DistPQ(:,i).^2)).*Di_DistPQ(:,i) ...
        - LogSwitch*(-Di_PDotQ(:,i) + Di_DistPQ(:,i))./(1- PDotQ(:,i) + DistPQ(:,i));
    Di_PhiPiNorm(i)    = (PhiPi(:,i))'*Di_PhiPi(:,i);
    Jac(:,i+SourceNum) = I(i)*(PhiPiNorm(i)*Di_PhiPi(:,i)/(PhiPiNorm(i))^2 - PhiPi(:,i).*Di_PhiPiNorm(i)/PhiPiNorm(i));
    %         Jac(:,i+1)  = Di_Estimate'*Descrepency;
    %         GradPseudo(i+SourceNum) = Jac(i+SourceNum) ;
    %% dObjectiveFunction/dTheta
    Di_CosGamma(:,i) = -sin(SourceTheta(i))*COS1'+cos(SourceTheta(i))*SIN1'...
        .*(cos(Mesh.PsiQ(:) - SourcePsi(i)))';
    Di_DistPQ(:,i) = -(SourceRadius(i)*Di_CosGamma(:,i)./DistPQ(:,i));
    Di_PDotQ(:,i)  = SourceRadius(i)*Di_CosGamma(:,i);
    Di_PhiPi(:,i)  = (-2./(DistPQ(:,i).^2)).*Di_DistPQ(:,i) ...
        - LogSwitch*(-Di_PDotQ(:,i) + Di_DistPQ(:,i))./(1- PDotQ(:,i) + DistPQ(:,i));
    Di_PhiPiNorm(i)    = (PhiPi(:,i))'*Di_PhiPi(:,i);
    Jac(:,i+2*SourceNum) = I(i)*(PhiPiNorm(i)*Di_PhiPi(:,i)/(PhiPiNorm(i))^2 - PhiPi(:,i).*Di_PhiPiNorm(i)/PhiPiNorm(i));
    %         Jac(:,i+2)  = Di_Estimate;
    %         GradPseudo(i+2*SourceNum) = Jac(i+2*SourceNum) ;
    
    %% dObjectiveFunction/dPsi
    Di_CosGamma(:,i) =  -sin(SourceTheta(i)).*SIN1'...
        .*(sin(SourcePsi(i) - Mesh.PsiQ(:)))';
    Di_DistPQ(:,i) = -(SourceRadius(i)*Di_CosGamma(:,i)./DistPQ(:,i));
    Di_PDotQ(:,i)  = SourceRadius(i)*Di_CosGamma(:,i);
    Di_PhiPi(:,i)  = (-2./(DistPQ(:,i).^2)).*Di_DistPQ(:,i) ...
        - LogSwitch*(-Di_PDotQ(:,i) + Di_DistPQ(:,i))./(1- PDotQ(:,i) + DistPQ(:,i));
    Di_PhiPiNorm(i)    = (PhiPi(:,i))'*Di_PhiPi(:,i);
    Jac(:,i+3*SourceNum) = I(i)*(PhiPiNorm(i)*Di_PhiPi(:,i)/(PhiPiNorm(i))^2 - PhiPi(:,i).*Di_PhiPiNorm(i)/PhiPiNorm(i));
    %         Jac(:,i+3)  = Di_Estimate'*Descrepency;
    %         GradPseudo(i+3*SourceNum) = Jac(i+3*SourceNum) ;
end


GradDescrepency = ConditionTerm*Jac'*Descrepency;
GradPseudo = GradDescrepency + GradSub';

% GradPseudo = GradPseudo';
end