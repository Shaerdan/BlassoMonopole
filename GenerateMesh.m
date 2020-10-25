function [Mesh] = GenerateMesh(MeshPointSSize,MeshPointQSize,Cap)
%UNTITLED Summary of this function goes here
%   Generate Mesh
%   Mesh.ThetaQ:  a 2d matrix meshing the theta commponent of the point Q \in R^2;
%   Mesh.PsiQ:    a 2d matrix meshing the psi commponent of the point Q \in R^2;
%   Mesh.RadiusS:   a 3d array(tensor) meshing of the radius commponent of the point S \in R^3;
%   Mesh.ThetaS:    a 3d array(tensor) meshing of the theta commponent of the point S \in R^3;
%   Mesh.PsiS:    a 3d array(tensor) meshing of psi commponent of the point S \in R^3;

SIncreSizeRadius = (Cap.RadiusUpper - Cap.RadiusLow)/MeshPointSSize.Radius;
SIncreSizeTheta = Cap.Theta/MeshPointSSize.Theta;
SIncreSizePsi = Cap.Psi/MeshPointSSize.Psi;
QIncreSizeTheta = pi/MeshPointQSize.Theta;
QIncreSizePsi = 2*pi/MeshPointQSize.Psi;
Mesh.ThetaQLine =  0:QIncreSizeTheta:pi;
Mesh.PsiQLine = 0:QIncreSizePsi:2*pi;
Mesh.RadiusS = Cap.RadiusLow:SIncreSizeRadius:Cap.RadiusUpper; 
Mesh.ThetaS = 0:SIncreSizeTheta:Cap.Theta;
Mesh.PsiS  =  0:SIncreSizePsi:Cap.Psi;
% Mesh.RadiusS = Cap.Radius*rand(MeshPointSSize.Radius,1); 
% Mesh.ThetaS = Cap.Theta*rand(MeshPointSSize.Theta,1);
% Mesh.PsiS  =  Cap.Psi*rand(MeshPointSSize.Psi,1);
[Mesh.ThetaQ,Mesh.PsiQ] = meshgrid(Mesh.ThetaQLine,Mesh.PsiQLine);
[Mesh.ThetaSpike,Mesh.PsiSpike,Mesh.RadiusSpike] = meshgrid(Mesh.ThetaS,Mesh.PsiS,Mesh.RadiusS);
Mesh.LengthRadiusS = length(Mesh.RadiusS); 
Mesh.LengthThetaS = length(Mesh.ThetaS); 
Mesh.LengthPsiS = length(Mesh.PsiS);

end

