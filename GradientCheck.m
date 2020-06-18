
clc;
clear all;
close all;
%% Setup Parameters:
%%% Regularization
Lambda = 0;
format longE

%%% Mesh Control
MeshPointSSize.Radius = 10; MeshPointSSize.Theta = 10; MeshPointSSize.Psi = 10;
MeshPointQSize.Theta = 100; MeshPointQSize.Psi = 100; Cap.Radius = 0.9;
Cap.Theta = 0.5*pi; Cap.Psi = 0.5*pi;

%% Souce Control:
RadiusReal = [0.8 0.8]; SourceNum = 2; % Real source depths and source number setup.
IntensityReal = [1 -1]; % Real source intensity setup.
ThetaReal = [.2 1.8]; PsiReal = [.2 1.8]; % Real source angle setup.
NoiseLevel = 0;
LocationReal = zeros(3*SourceNum,1);   % Initialize real source location assembly vector.
LocationReal (1:SourceNum) = RadiusReal; % Store radius component to the assembly vector.
LocationReal (SourceNum+1:2*SourceNum) = ThetaReal; % Store theta component ....
LocationReal (2*SourceNum+1:end) = PsiReal; % Store psi component ...

%% Create the mesh for point Q in R^3 and point S in R^2 (Q,S as defined
% in the report).
% Mesh:
% Mesh.RadiusQ -- Radius component of Q mesh.
% Mesh.ThetaQ -- Theta component of Q mesh.
% Mesh.PsiQ -- Psi component of Q mesh.
% Mesh.ThetaS -- Theta component of S mesh.
% Mesh.PsiS -- Psi component of S mesh.
Mesh = GenerateMesh(MeshPointSSize,MeshPointQSize,Cap);

%% Create measurement data.

% Data: Data.Measurement, Data.Noise
Data = GenerateData(IntensityReal,LocationReal,Mesh,NoiseLevel);
S_exact = [IntensityReal'; LocationReal];
% S1 = randn(4*SourceNum,1);
S1 = S_exact+(1e-2)*randn(4*SourceNum,1);
delta = 1e-10;
f = @(X)GradTest(X,Data.Measurement,Mesh,Lambda);
S2 = S1;
S0 = S1;
for i = 1:4*SourceNum
    S2(i) = S1(i) + 0.5*delta;
    S0(i) = S1(i) - 0.5*delta;
    y2 = f(S2); y1 = f(S1); y0 =f(S0);
    gradcentral(i) = (y2-y0)/(delta);
    S2 = S1;
    S0 = S1;
end

S2 = S1;
S0 = S1;
S3 = S1;
Sm1= S1;
for i = 1:4*SourceNum
    S2(i) = S1(i) + delta;
    S0(i) = S1(i) - delta;
    S3(i) = S1(i) + 2*delta;
    Sm1(i)= S1(i) - 2*delta; 
    y3=f(S3); y2 = f(S2); y0 =f(S0);ym1 = f(Sm1);
    gradhigher(i) = (-y3 + 8*y2 - 8*y0 +ym1)/(12*delta);
    S2 = S1;
    S0 = S1;
    S3 = S1;
    Sm1 = S1;
end

dydx_analytic = GradPseudoFminunc(S1,Data.Measurement,Mesh,Lambda);