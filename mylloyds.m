function [partition, codebook, distor, rel_distor] = mylloyds(bit, tol, plot_flag)
%MYLLOYDS Optimize quantization parameters for the standard Gaussian distribution by using the Lloyd algorithm.
%   [PARTITION, CODEBOOK] = MYLLOYDS(BIT) optimizes the
%   scalar quantization PARTITION and CODEBOOK based on the standard
%   Gaussian distribution using the Lloyd algorithm.
%   BIT is the bit of quantization and should be a scalar positive integer. The optimized
%   CODEBOOK is a vector of length equal to 2^BIT. PARTITION is a vector of
%   length equal to 2^BIT minus 1. The optimization will be terminated if the relative distortion
%   is less than 10^(-7).
%
%   [PARTITION, CODEBOOK] = MYLLOYDS(BIT, TOL) provides the
%   tolerance in the optimization.
%
%   [PARTITION, CODEBOOK, DISTORTION] = MYLLOYDS(...) outputs the final distortion
%   value in terminating the computation.
%
%   [PARTITION, CODEBOOK, DISTORTION, REL_DISTORTION] = MYLLOYDS(...) outputs the
%   relative distortion value in terminating the computation.

%   Written by Bin Hu, learned from lloyds.m (Copyright 1996-2005 The MathWorks, Inc.) in MATLAB.
%   $Revision: 1.0.0.0 $ $Date: 2017/07/12 $

% validation verification and format conversion.
error(nargchk(1,3,nargin, 'struct'));

if length(bit) < 1
    error('mylloyds:EmptyBIT','Bit parameter cannot be empty.');
elseif (min(size(bit)) > 1) || (max(size(bit)) > 1)
    error('mylloyds:NonScalarBIT','Bit parameter must be a scalar.');
elseif ~isreal(bit) || ~(bit == fix(bit)) || (bit < 1)
    error('mylloyds:InvalidBIT',['Bit parameter must be a positive ', ...
                        'integer.']);
else
    if nargin < 2 || isempty(tol)
        tol = 10^(-7);
    end;
end;

% initial half of codebook for the standard Gaussian distribution
len_codebook_half = 2^bit / 2;
codebook_half = linspace(0.03, 4, len_codebook_half);

% initial half of partition
partition_half = (codebook_half(2 : len_codebook_half) + codebook_half(1 : len_codebook_half-1)) / 2;

% distortion initialization, computation
distor = 0;
partition_half_ex = [0, partition_half, inf];
for ii = 1 : len_codebook_half
    distor = distor + (1/sqrt(2*pi)) * quadgk(@(x) (x - codebook_half(ii)).^2 .* exp(-x.^2/2), partition_half_ex(ii), partition_half_ex(ii+1));
end;
distor = distor * 2;
last_distor = 0;
ter_cond2 = eps;
if distor > ter_cond2
    rel_distor = abs(distor - last_distor)/distor;
else
    rel_distor = distor;
end;

while (rel_distor > tol) && (rel_distor > ter_cond2)
    % using the centroid condition, find half of the optimal codebook.
    for jj = 1 : len_codebook_half
        codebook_half(jj) = (1/sqrt(2*pi)) * (exp(-(partition_half_ex(jj))^2/2)-exp(-(partition_half_ex(jj+1))^2/2)) / (qfunc(partition_half_ex(jj))-qfunc(partition_half_ex(jj+1)));
    end;

    % compute half of partition
    partition_half = (codebook_half(2 : len_codebook_half) + codebook_half(1 : len_codebook_half-1)) / 2;

    % compute distortion
    distor = 0;
    partition_half_ex = [0, partition_half, inf];
    for ii = 1 : len_codebook_half
        distor = distor + (1/sqrt(2*pi)) * quadgk(@(x) (x - codebook_half(ii)).^2 .* exp(-x.^2/2), partition_half_ex(ii), partition_half_ex(ii+1));
    end;
    distor = distor * 2;
    if distor > ter_cond2
        rel_distor = abs(distor - last_distor)/distor;
    else
        rel_distor = distor;
    end;
end;

% output
partition = [-inf, -fliplr(partition_half), 0, partition_half, inf];
codebook = [-fliplr(codebook_half), codebook_half];

% plot
if (nargin > 2) && plot_flag
    x_len = 5;
    partition_plot = [-fliplr(partition_half), 0, partition_half];
    x = -x_len:0.01:x_len;
    x_codebook = [codebook(:)'; codebook(:)'; codebook(:)'];
    x_partition = [partition_plot(:)'; partition_plot(:)'; partition_plot(:)'];
    
    y_max = 0.4;
    len_codebook = length(codebook);
    y_codebook = [zeros(1, len_codebook); ones(1, len_codebook)*y_max; zeros(1, len_codebook)];
    y_codebook = y_codebook(:);
    y_partition = y_codebook(1:3*len_codebook-3);
    
    plot(x, normpdf(x, 0, 1), x_codebook(:), y_codebook, x_partition(:), y_partition)
    axis([-x_len, x_len, 0, y_max])
    title('Optimized quantization partition and codebook for the standard Gaussian distribution.')
    xlabel('Blue: standard Gaussian distribution; Green: codebook; Red: partition.');
end;
