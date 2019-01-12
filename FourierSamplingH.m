function x = FourierSamplingH(y, N, Loc)

t = zeros(N^2, 1);
t(Loc) = y;
t = reshape(t, N, N);
x = ifft2(t);
