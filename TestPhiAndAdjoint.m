clc;
clear all;
close all;

%% Setup Parameters:
%%% Regularization
SwitchNorm = 0;
%%% Mesh Control
MeshPointSSize.Radius = 10; MeshPointSSize.Theta = 10; MeshPointSSize.Psi = 10;
MeshPointQSize.Theta = 500; MeshPointQSize.Psi = 500; Cap.Radius = 0.9999;
Cap.Theta = pi; Cap.Psi = 2*pi;
%%  Normalization switch: 0 -- disable  1 -- enable.
Switch.Normalize = 0;
TestType = 'constant'; TestPoints = 'single';
%% Create the mesh for point Q in R^3 and point S in R^2 (Q,S as defined
% in the report).
% Mesh:
% Mesh.RadiusQ -- Radius component of Q mesh.
% Mesh.ThetaQ -- Theta component of Q mesh.
% Mesh.PsiQ -- Psi component of Q mesh.
% Mesh.ThetaS -- Theta component of S mesh.
% Mesh.PsiS -- Psi component of S mesh.
Mesh = GenerateMesh(MeshPointSSize,MeshPointQSize,Cap);
[Mesh.STheta,Mesh.SPsi,Mesh.SRadius] = meshgrid(Mesh.ThetaS,Mesh.PsiS,Mesh.RadiusS);
switch TestType
    case 'random'
        Q = 10*randn(size(Mesh.ThetaQ));
        P = 1*randn(size(Mesh.STheta));
    case 'constant'
        Q = ones(size(Mesh.ThetaQ));
        P = ones(size(Mesh.STheta));
end
switch TestPoints
    case 'multi'
        PhiSQ = PhiStar_Q(Mesh,Q,SwitchNorm);
        SIN1 = sin(Mesh.STheta);
        R1 = Mesh.SRadius.^2;
        Integrand1 = P.*PhiSQ.*R1.*SIN1;
        dx1 = Mesh.ThetaS'; dy1 = Mesh.PsiS'; dz1 = Mesh.RadiusS';
        Inner_PhiSQ_P_L2Omega = trapz(dx1,trapz(dy1,trapz(dz1,Integrand1,3),2),1);
        PhiP = Phi_P(Mesh,P);
        dx2 = Mesh.ThetaQLine'; dy2 = Mesh.PsiQLine';
        Integrand2 = PhiP.*Q.*sin(Mesh.ThetaQ);
        Inner_Q_PhiP_L2PartialOmega = trapz(dx2,trapz(dy2,Integrand2,2));

    case 'single'
        S.ThetaS = pi*rand;
        S.PsiS = 2*pi*rand;
        S.RadiusS = 0.999*rand;   
%          S.ThetaS =0;
%          S.PsiS = 0;
%          S.RadiusS = 0.99;
        PhiSQ = PhiStar_Q_SingleS(Mesh,S,Q,SwitchNorm);
        fun = @(x,y) 2*((1+S.RadiusS^2 - 2*S.RadiusS*(cos(x))).^(-1/2)).*sin(x) ...
            - log(1-S.RadiusS*cos(x) + (1+S.RadiusS^2 - 2*S.RadiusS*(cos(x))).^(1/2)).*sin(x);
        format long
        q = integral2(fun,0,pi,0,2*pi);
        2*4*pi-q

end