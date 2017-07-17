% Main program for quantization.

% Written by Bin Hu
% Date: 2017/07/16

clear all
close all
% clc

% input
bit = 5;
codebook = [];
tol = eps;
plot_flag = 1;

% filename (.mat) of Lloyd quantization for the standard Gaussian
% distribution
lloydsfile = ['lloyds_', num2str(bit), 'bit'];
load(lloydsfile, 'codebook'); % comment if no such file exists

% Lloyd quantization for the standard Gaussian distribution
[partition, codebook, distor, rel_distor] = mylloyds(bit, codebook, tol, plot_flag);

% save the codebook and partition in the file (.mat), or overwrite it if it
% already exists
save(lloydsfile,'codebook','partition');
