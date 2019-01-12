function y = FourierSampling(x, Loc)
% input size: NxN

y = fft2(x);
y = y(Loc);
y = y(:);
