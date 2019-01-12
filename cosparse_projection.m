function [x_proj,res] = cosparse_projection(x, Omega, OmegaH, Lambda)

N = size(x,1);
Omega_x = Omega(x);
Omega_x_Lambda = zeros(size(Omega_x));
Omega_x_Lambda(Lambda) = Omega_x(Lambda);
OmegaLHOmegaL_op = @(b) OmegaLHOmegaL(b, Omega, OmegaH, Lambda, N);
[x_p_proj, res, iter] = cgsolve(OmegaLHOmegaL_op, vec(OmegaH(Omega_x_Lambda)), 1e-6, 1000, 0);
x_p_proj = reshape(x_p_proj,N,N);
x_proj = x - x_p_proj;


function w = OmegaLHOmegaL(x, Omega, OmegaH, L, N)
x = reshape(x,N,N);
v = Omega(x);
vt = zeros(size(v));
vt(L) = v(L);
w = OmegaH(vt);
w = vec(w);

