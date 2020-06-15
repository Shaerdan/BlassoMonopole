clc;
clear all;
close all;
%% Setup Parameters:
%%% Regularization
Lambda = 1e-7;

%%% Mesh Control
MeshPointSSize.Radius = 10; MeshPointSSize.Theta = 10; MeshPointSSize.Psi = 10;
MeshPointQSize.Theta = 100; MeshPointQSize.Psi = 100; Cap.Radius = 0.9;
Cap.Theta = 0.5*pi; Cap.Psi = 0.5*pi;

%% Souce Control:
RadiusReal = [0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7]; SourceNum = 8; % Real source depths and source number setup.
IntensityReal = [1 -1 1 -1 1 -1 1 -1]; % Real source intensity setup.
ThetaReal = [.2 .25 0.7 0.75 1.2 1.25 1.5 1.55]; PsiReal = [.2 .25 0.7 0.75 1.2 1.25 1.5 1.55]; % Real source angle setup.
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

%% Choose method:
method = 'fmincon';       % 'fmincon' or 'lbfgsc'

%% fminunc solver control:
Opts = optimoptions(@fminunc,'PlotFcns',{@optimplotfval,@optimplotx},...
    'MaxIterations',10000,'MaxFunctionEvaluations',20000,...
    'OptimalityTolerance',1e-16,'StepTolerance',1e-16);


%% Initializing %%%%%%

InitialSolution = randn(4*SourceNum,1);

IntensityMax = 5;
LowerBoundNonLinear = [-IntensityMax*ones(1,SourceNum), ...
    repmat(0.001,1,SourceNum),zeros(1,SourceNum),...
    zeros(1,SourceNum)];
UpperBoundNonLinear = [IntensityMax*ones(1,SourceNum),...
    repmat(Cap.Radius,1,SourceNum),repmat(pi,1,SourceNum),...
    repmat(pi,1,SourceNum)];

switch method
    case 'fmincon'
        problem = createOptimProblem('fmincon','x0',InitialSolution,...
            'ub',UpperBoundNonLinear,'lb',LowerBoundNonLinear, ...
            'objective',@(X) ObjectiveFuncNonLinearFminunc(X,Data.Measurement,...
            Lambda, Mesh),...
            'options', Opts);
        [SolutionArgmin,fval2,flag,stepcount] = fmincon(problem);
        fprintf('flag of fmincon %d\n',flag);
        Solution.Intensity = SolutionArgmin(1:SourceNum)';
        Solution.Location = SolutionArgmin(SourceNum+1:end)';
    case 'lbfgsc'
        OptLBFGSB = struct('x0',InitialSolution,'m',20,'factr',1e4,...
            'pgtol',1e-10,'maxIts',10000,'maxTotalIts',500000,'printEvery',1);
        [SolutionArgmin,f_val2,info2] = lbfgsb( @(X) ...
            ObjectiveFuncNonLinearLBFGSB(X,Data.Measurement,Lambda, Mesh),...
            LowerBoundNonLinear', UpperBoundNonLinear', OptLBFGSB );
        Solution.Intensity = SolutionArgmin(1:SourceNum)';
        %             if (ii>4)
        %                 Solution.Intensity(abs(Solution.Intensity)<= Cropping) = 0;
        %             end
        Solution.Location = SolutionArgmin(SourceNum+1:end)';
end

IntensitySoln = Solution.Intensity;
RadiusSoln = Solution.Location(1:SourceNum);
ThetaSoln = Solution.Location(SourceNum+1:2*SourceNum);
PsiSoln = Solution.Location(2*SourceNum+1:3*SourceNum);
figure(6)
CartesianXSoln = RadiusSoln.*sin(ThetaSoln).*cos(PsiSoln);
CartesianYSoln = RadiusSoln.*sin(ThetaSoln).*sin(PsiSoln);
CartesianZSoln = RadiusSoln.*cos(ThetaSoln);
CartesianXReal = RadiusReal.*sin(ThetaReal).*cos(PsiReal);
CartesianYReal = RadiusReal.*sin(ThetaReal).*sin(PsiReal);
CartesianZReal = RadiusReal.*cos(ThetaReal);
[XSphere,YSphere,ZSphere] = sphere;
h = surf(XSphere, YSphere, ZSphere,'FaceColor','none');
set(h, 'FaceAlpha', 0.2);
hold on;
for i = 1:SourceNum
    scatter3(CartesianXReal(i),CartesianYReal(i),CartesianZReal(i),100,'r','linewidth',1.5);
end
for i = 1:SourceNum
    scatter3(CartesianXSoln(i),CartesianYSoln(i),CartesianZSoln(i),100,'k','*','linewidth',1.5);
    hold on;
end
hold off;


