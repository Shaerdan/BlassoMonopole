clc;
clear all;
close all;

%% Setup Parameters:
%%% Regularization
Lambda = 1e-16;

%%% Mesh Control
MeshPointSSize.Radius = 10; MeshPointSSize.Theta = 20; MeshPointSSize.Psi = 20;
MeshPointQSize.Theta = 1000; MeshPointQSize.Psi = 1000; Cap.Radius = 0.9;
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
FBS.GradTol = 1e-10; FBS.StepTol = 1e-10;
FBS.Iteration = 1000;
FBS.DisplayFrequency = 1;
FBS.Mu = 0.4;

%% ADMM solver control:
ADMM.Tol = 1e-8; ADMM.Steps = 1000; ADMM.Rho = 1;


casestudy = 'monopole';
switch casestudy
    %% Test for source reconstruction problem:
    % Initialize the linear minimization problem:
    case 'monopole'
        InitialIntensity = 100*randn(SourceNum,1);
        PhiComponent = ComputePotentialComponent(LocationReal,Mesh);
        FL2Norm = L2NormF(LocationReal,Mesh);
        [M,b]= ComputeHessian(Data.Measurement,PhiComponent,FL2Norm,Mesh);
        LipschitzConst = norm(M,2);
        FBS.Tau = 1/LipschitzConst;
        IntensitySoln = FBSSolver(InitialIntensity,M,b,Data.Measurement,...
            FBS,PhiComponent,Mesh,Lambda,FL2Norm');
        IntensitySoln1 = admm(InitialIntensity,M,b,Lambda,ADMM);
        
        %% General Test: solve 0.5*||AX-b||_{2}^{2} + ||X||_{1};
    case 'general'
        m=10;n=400;
        x = randn(m,1);
        A = randn(n,m);
        x0 = randn(m,1);
        Data.Measurement = A*x;
        figure(1)
        plot(Data.Measurement)
        FL2Norm = ones(m,1);
        M = A'*A;
        b = A'*Data.Measurement;
        chk0 = M*randn(m,1) - b;
        LipschitzConst = norm(M,2);
        FBS.Tau = 1/LipschitzConst;
        Soln = FBSTest(x0,M,b,Data.Measurement,FBS,A,Lambda);
        figure(2)
        plot(x); hold on; plot(Soln,'-*');
        res=norm(x-Soln,2);
        Soln1 =admm(x0,M,b,Lambda,ADMM) ;
end