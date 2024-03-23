close all; clearvars;

% INPUT DATA

% Mesh discretization and subvolume dimensions
nx = 100;
ny = 100;

% Optimization parameters
volfrac = 0.50;
penal   = 3;
rfil    = [];  % undefined filter radius
ft      = [];  % undefined filter type

% RUN TOPOLOGY OPTIMIZATION CODE (TopMatFVT)
TopMatFVT(nx,ny,volfrac,penal,rfil,ft);

% RUN topX code
% rfil = 1;
% ft = 1;
% topX(nx,ny,volfrac,penal,rfil,ft);