function w = FDAnalysis(x)
% input size: NxN
% output size: [Nx(N-1) ; ((N-1)xN)'] = 2Nx(N-1)

w = [x(:, 2:end)-x(:, 1:(end-1)); (x(2:end, :)-x(1:(end-1), :)).'];

