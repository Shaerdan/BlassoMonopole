clc;
clear all;
close all;

%% Setup Parameters:
%%% Regularization
Lambda = 1e-16;

%%% Mesh Control
MeshPointSSize.Radius = 10; MeshPointSSize.Theta = 20; MeshPointSSize.Psi = 20;
MeshPointQSize.Theta = 100; MeshPointQSize.Psi = 100; Cap.Radius = 0.9;
Cap.Theta = pi; Cap.Psi = 2*pi;

%% Souce Control:
RadiusReal = 0.8; SourceNum = 2; % Real source depths and source number setup.
IntensityReal = [1 -1]; % Real source intensity setup.
ThetaReal = [0.2, 0.3]; PsiReal = [0.2, 0.3]; % Real source angle setup.
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

%% FBS solver control:
FBS.Tau = 1e-2;
FBS.GradTol = 1e-10; FBS.StepTol = 1e-10;
FBS.Iteration = 2e5;
FBS.DisplayFrequency = 100;
FBS.Mu = 0.4;

%% FBS solver:
% Initialize the linear minimization problem:
InitialIntensity = randn(SourceNum,1);
PhiComponent = ComputePotentialComponent(LocationReal,Mesh);
FL2Norm = L2NormF(LocationReal,Mesh);
LipschitzConst = ComputeLipschitzConstant(Mesh,PhiComponent,FL2Norm,...
    SourceNum);
FBS.tau = 1/LipschitzConst;
[M,b]= ComputeHessian(Data.Measurement,PhiComponent,FL2Norm,Mesh);
IntensitySoln = FBSSolver(InitialIntensity,M,b,Data.Measurement,...
    FBS,PhiComponent,Mesh,Lambda,FL2Norm');