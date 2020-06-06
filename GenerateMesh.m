function [Mesh] = GenerateMesh(MeshPointSSize,MeshPointQSize,RadiusCap)
%UNTITLED Summary of this function goes here
%   Generate Mesh
%   Mesh.ThetaQ:  a 2d matrix meshing the theta commponent of the point Q \in R^2;
%   Mesh.PsiQ:    a 2d matrix meshing the psi commponent of the point Q \in R^2;
%   Mesh.RadiusS:   a 3d array(tensor) meshing of the radius commponent of the point S \in R^3;
%   Mesh.ThetaS:    a 3d array(tensor) meshing of the theta commponent of the point S \in R^3;
%   Mesh.PsiS:    a 3d array(tensor) meshing of psi commponent of the point S \in R^3;

SIncreSizeRadius = 1/MeshPointSSize.Radius;
SIncreSizeTheta = pi/MeshPointSSize.Theta;
SIncreSizePsi = 2*pi/MeshPointSSize.Psi;
QIncreSizeTheta = pi/MeshPointQSize.Theta;
QIncreSizePsi = 2*pi/MeshPointQSize.Psi;
Mesh.ThetaQLine =  0:QIncreSizeTheta:pi;
Mesh.PsiQLine = 0:QIncreSizePsi:2*pi;
Mesh.RadiusS = 0:SIncreSizeRadius:RadiusCap; 
Mesh.ThetaS = 0:SIncreSizeTheta:pi;
Mesh.PsiS  =  0:SIncreSizePsi:2*pi;
[Mesh.ThetaQ,Mesh.PsiQ] = meshgrid(Mesh.ThetaQLine,Mesh.PsiQLine);
Mesh.LengthRadiusS = length(Mesh.RadiusS); 
Mesh.LengthThetaS = length(Mesh.ThetaS); 
Mesh.LengthPsiS = length(Mesh.PsiS);

end

