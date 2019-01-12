function x = FDSynthesis(w)
% input size: [Nx(N-1) ; ((N-1)xN)'] = 2Nx(N-1)
% output size: NxN
n = size(w,2)+1;

x = zeros(n, n);
x(:, 2:end) = w(1:n, :);
x(:, 1:(end-1)) = x(:, 1:(end-1)) - w(1:n, :);
x(2:end, :) = x(2:end, :) + w((n+1):end, :).';
x(1:(end-1), :) = x(1:(end-1), :) - w((n+1):end, :).';

