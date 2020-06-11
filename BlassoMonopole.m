clc;
clear all;
close all;
%% Setup Parameters:
%%% Regularization
Lambda = 1e-4;

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

%% Choose solvers:
% Linear solver:
% 'FBS' -- linear proximal forward backward splitting;
% 'fminunc' -- matlab solver with a default quasi-newton algorithm
% and bfgs approximation of the Hessian.
LinearSolver = 'fbs';               % choose 'fbs' or 'relaxedfbs' or 'fista' or 'fmincon';
% Set the nonlinear solver:
% 'lbfgsc' -- limited memory bfgs solver.
% 'fminunc' -- matlab solver with a default quasi-newton algorithm
%  and bfgs approximation of the Hessian.
NonLinearSolver = 'lbfgsc';        % choose 'lbfgsc' or 'fmincon';

%% FBS solver control:
FBS.Tau = 1e-2;
FBS.GradTol = 1e-10; FBS.StepTol = 1e-16;
FBS.Iteration = 2e5;
FBS.DisplayFrequency = 100;
FBS.Mu = 0.4;

%% FISTA solver constrol:
FISTAopts.lambda = Lambda;
FISTAopts.max_iter = 1e5;
FISTAopts.tol =1e-10;
FISTAopts.verbose = 'true';
FISTAopts.L = 0;

%% fminunc solver control:
Opts = optimoptions(@fminunc,'PlotFcns',{@optimplotfval,@optimplotx},...
    'MaxIterations',10000,'MaxFunctionEvaluations',20000,...
    'OptimalityTolerance',1e-10,'StepTolerance',1e-16);


%% Initializing the Sliding Frank Wolfe iteration %%%%%%
GlobalIteration = 100;

Solution.Intensity= ([]); Solution.Location = ([]);
% Solution.Radius .Theta .Psi below is used to update the Location assembly vector.
Solution.Radius = ([]); Solution.Theta = ([]); Solution.Psi = ([]);
% Initialize EtaMax.
EtaMax = zeros(1,GlobalIteration);


%% Sliding Frank Wolfe Algorithm Interation:
for ii = 1:GlobalIteration
    
    FL2Norm = L2NormF(Solution.Location,Mesh);
    
    %%  Solve (r,theta,psi) = Argmax_(r',theta',psi') {Eta(r',theta',psi')}
    %          Eta = Eta(Mesh.S) = integral (
    %          F_Normalized(Mesh.Q,Mesh.S)*P(Mesh.Q)dMesh.Q ),
    %          dMesh.Q = sin(Mesh.ThetaQ) * dMesh.ThetaQ * dMesh.PsiQ.
    P   = EvaluateP(Solution.Location,Solution.Intensity,Data.Measurement,Mesh,FL2Norm);
    %%% MARKED FOR DEBUGGING %%%
    tic
    Eta = EvaluateEta(Lambda,Mesh,P);
    toc
    [EtaMax(ii),Index] = max(Eta(:));
    [IndexRadius,IndexTheta,IndexPsi] = ind2sub(size(Eta),Index);
    
    %%% Debug line, output EtaMax:
    fprintf('Eta= %f', EtaMax(ii));
    pause(1);
    %%% MARKED FOR DEBUGGING %%%
    
    % obtain the new source location from the index of EtaMax.
    RadiusUpdate = Mesh.RadiusS(IndexRadius);
    ThetaUpdate = Mesh.ThetaS(IndexTheta);
    PsiUpdate = Mesh.PsiS(IndexPsi);
    %%% Refine Steps: solve argmin[-(1/2)*Eta^2]:
    LowerBoundRefine = [0, 0, 0]';
    UpperBoundRefine = [Cap.Radius, Cap.Theta, Cap.Psi]';
    InitialSolution = [RadiusUpdate,ThetaUpdate,PsiUpdate]';
    %     OptLBFGSB = struct('x0',InitialSolution,'m',20,'factr',1e0,...
    %         'pgtol',1e-30,'maxIts',10000,'maxTotalIts',500000,'printEvery',1);
    %     [SolutionLoc,f_val2,info2] = lbfgsb( @(X) ObjEta(X,Lambda,Mesh,P),...
    %         LowerBoundRefine,UpperBoundRefine,OptLBFGSB );
    problem = createOptimProblem('fmincon','x0',InitialSolution,...
        'ub',UpperBoundRefine,'lb',LowerBoundRefine, ...
        'objective',@(X) ObjEta(X,Lambda, Mesh,P),...
        'options', Opts);
    [SolutionLoc,EtaMin] = fmincon(problem);
    RadiusUpdate = SolutionLoc(1);
    ThetaUpdate = SolutionLoc(2);
    PsiUpdate = SolutionLoc(3);
    figure(2)
    bar([RadiusUpdate, ThetaUpdate, PsiUpdate]);
    
    %%% Termination condition for the global interation
    if (EtaMax(ii) <=1)
        break
    end
    
    % Update the source location assembly vector:
    Solution.Radius = [Solution.Radius,RadiusUpdate];
    Solution.Theta = [Solution.Theta, ThetaUpdate];
    Solution.Psi = [Solution.Psi, PsiUpdate];
    Solution.Location = [Solution.Radius, Solution.Theta, Solution.Psi];
    % End of (r,theta,psi) = Argmax_(r',theta',psi'){Eta(r',theta',psi')} Procedure
    
    %% Update the source intensity by solving the linear lasso problem
    
    % Initialize the linear minimization problem:
    IntensityUpdate = 0;
    InitialIntensity = [Solution.Intensity, IntensityUpdate];
    PhiComponent = ComputePotentialComponent(Solution.Location,Mesh);
    FL2Norm = L2NormF(Solution.Location,Mesh);
    SourceNumUpdate = ii;
    % Choose solvers:
    switch LinearSolver
        case 'fbs'
            LipschitzConst = ComputeLipschitzConstant(Mesh,PhiComponent,FL2Norm,...
                SourceNumUpdate);
            FBS.tau = 1/LipschitzConst;
            Solution.Intensity = FBSSolver(InitialIntensity,Data.Measurement,...
                FBS,PhiComponent,Mesh,Lambda,FL2Norm);
        case 'fista'
            LipschitzConst = ComputeLipschitzConstant(Mesh,PhiComponent,FL2Norm,...
                SourceNumUpdate);
%             LipschitzConst = 10*LipschitzConst;
            Solution.Intensity = fista_general(@(X) GradFISTA(X,Data.Measurement,...
                Mesh,PhiComponent,FL2Norm),@(X) proj_l1(X,Lambda), InitialIntensity',...
                LipschitzConst,FISTAopts, @(X) CalcDiscrepancy(X,Data.Measurement,...
                PhiComponent,FL2Norm,Mesh));
            Solution.Intensity = Solution.Intensity';
        case 'relaxedfbs'
            LipschitzConst = ComputeLipschitzConstant(Mesh,PhiComponent,FL2Norm,...
                SourceNumUpdate);
            FBS.tau = 1/LipschitzConst;
            Solution.Intensity = RelaxedFBSSolver(InitialIntensity,Data.Measurement,...
                FBS,PhiComponent,Mesh,Lambda,FL2Norm);
        case 'fmincon'
            IntensityMax = 5;
            lb0 = -IntensityMax*ones(1,SourceNumUpdate)';    % Lower bound = -S_max* RealIntensity
            ub0 = IntensityMax*ones(1,SourceNumUpdate)';     % Upper bound =  S_max* RealIntensity
            problem = createOptimProblem('fmincon','x0',InitialIntensity,'ub',ub0,'lb',lb0, ...
                'objective',@(X) ObjectiveFuncLinearLasso(X,Data.Measurement,PhiComponent,...
                Mesh,FL2Norm,Lambda),'options',Opts);
            [Solution.Intensity] = fmincon(problem);
    end
    
    % Debug line, output intensities:
    figure(3)
    bar(Solution.Intensity);
    % End of source intensity update
    
    %% Update source intensity and locations by solving anonlinear Lasso problem:
    % Initialize the nonlinear minimization problem:
    InitialSolution = [Solution.Intensity,Solution.Location];
    SourceNumUpdate = length(Solution.Intensity);
    % The box constraint is set to be: -IntensityMax < Intensity <
    % IntensityMax, 0 < Radius < 0.999, 0 < Theta < pi, 0 < Psi < 2*pi:
    IntensityMax = 5;
    LowerBoundNonLinear = [-IntensityMax*ones(1,SourceNumUpdate), ...
        repmat(0.001,1,SourceNumUpdate),zeros(1,SourceNumUpdate),...
        zeros(1,SourceNumUpdate)];
    UpperBoundNonLinear = [IntensityMax*ones(1,SourceNumUpdate),...
        repmat(Cap.Radius,1,SourceNumUpdate),repmat(pi,1,SourceNumUpdate),...
        repmat(2*pi,1,SourceNumUpdate)];
    
    % Choose nonlinear solver:
    switch NonLinearSolver
        case 'lbfgsc'
            OptLBFGSB = struct('x0',InitialSolution','m',10,'factr',1e0,...
                'pgtol',1e-30,'maxIts',10000,'maxTotalIts',500000,'printEvery',1);
            [SolutionArgmin,f_val2,info2] = lbfgsb( @(X) ...
                ObjectiveFuncNonLinearLBFGSB(X,Data.Measurement,Lambda, Mesh),...
                LowerBoundNonLinear', UpperBoundNonLinear', OptLBFGSB );
            Solution.Intensity = SolutionArgmin(1:SourceNumUpdate)';
            Solution.Location = SolutionArgmin(SourceNumUpdate+1:end)';
        case 'fmincon'
            problem = createOptimProblem('fmincon','x0',InitialSolution,...
                'ub',UpperBoundNonLinear,'lb',LowerBoundNonLinear, ...
                'objective',@(X) ObjectiveFuncNonLinearFminunc(X,Data.Measurement,...
                Lambda, Mesh),...
                'options', Opts);
            [SolutionArgmin,fval2,flag,stepcount] = fmincon(problem);
            fprintf('flag of fmincon %d\n',flag);
            Solution.Intensity = SolutionArgmin(1:SourceNumUpdate);
            Solution.Location = SolutionArgmin(SourceNumUpdate+1:end);
    end
    
    % Update Solution Location Components:
    Solution.Radius = Solution.Location(1:ii);
    Solution.Theta  = Solution.Location(ii+1:2*ii);
    Solution.Psi    = Solution.Location(2*ii+1:3*ii);
    % Debug line, output the source intensities and locations:
    %     figure(4)
    %     bar(Solution.Intensity);
    %     figure(5)
    %     bar(Solution.Location);
    % End of source locations and intensities update.
    
    %% Display The Solution:
    IntensitySoln = Solution.Intensity;
    RadiusSoln = Solution.Location(1:SourceNumUpdate);
    ThetaSoln = Solution.Location(SourceNumUpdate+1:2*SourceNumUpdate);
    PsiSoln = Solution.Location(2*SourceNumUpdate+1:3*SourceNumUpdate);
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
    for i = 1:SourceNumUpdate
        scatter3(CartesianXSoln(i),CartesianYSoln(i),CartesianZSoln(i),100,'k','*','linewidth',1.5);
        hold on;
    end
    hold off;
    
end





