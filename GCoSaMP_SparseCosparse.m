function [best_X1,best_alpha_hat,best_X2,errTrace,best_iter] = GCoSaMP_SparseCosparse(y, A, At, D, Dt, Omega, Omegat, l, s, tol, maxIter, tol_early)

% Relaxation of the generalized CoSaMP for Synthesis-Analysis combined model
% Reference: "Generalizing CoSaMP to Signals from a Union of Low Dimensional Linear Subspaces"
% Authors: Tom Tirer and Raja Giryes
% Journal: Applied and Computational Harmonic Analysis (accepted 2018)

%%%%
% A = measurements operator (function handle)
% Omega = analysis operator (function handle)
% D = dictionary synthesis operator (function handle)
% tol, maxIter are used for termination
% tol_early is used for early termination
%%%%%%%

N_cg_iters_for_LS = 60;
cosparse_constraint_penalty_for_LS  = 1;
max_non_improvement_counter = 2; % before early termination
a_sparse = 2;
a_cosparse = 1;

temp = At(y);
rows = size(temp,1); cols = size(temp,2);
rowsAnalysis = size(Omega(temp),1); colsAnalysis = size(Omega(temp),2);
p = rowsAnalysis*colsAnalysis;
x_sparse = zeros(rows,cols);
x_cosparse = zeros(rows,cols);

r = y;
T = []; % support
Lambda = [1:p]; % cosupport
t = 0;
counter = 0;
err = 1e6;
best_err = inf;
errTrace = [];

while ((err>tol) && (t <= maxIter))
    
    t = t + 1;
    
    v = real(At(r));
    
    % find new support & cosupport elements
    Dt_v = Dt(v);
    alpha_siz = size(Dt_v);
    [~,indices] = sort(abs(Dt_v(:)),'descend');
    T_delta = indices(1:(a_sparse*s));
    Omega_v = Omega(v);
    [~,indices] = sort(abs(Omega_v(:)),'ascend');
    Lambda_delta = indices(1:floor(a_cosparse*l));

    % Merge supports, intersect cosupports
    T_tilde = union(T,T_delta);
    Lambda_tilde = intersect(Lambda,Lambda_delta);
    
    % LS estimation
    % Note that there is a relaxation (with lagrange_mult=1) of the cosparse subspace constraint
    A_LS = @(z) fullAt_fullA(z,A,At,D,Dt,Omega,Omegat,T_tilde,Lambda_tilde,rows,cols,alpha_siz,cosparse_constraint_penalty_for_LS);
    temp = vec(Dt(At(y)));
    b_LS = [temp(T_tilde); vec(At(y))];
    [x_LS, res, iter] = cgsolve(A_LS, b_LS, 1e-6, N_cg_iters_for_LS, 0);
    x_LS = real(x_LS); % for A that produce complex measurements (like FFT)
    
    alpha_tilde = zeros(alpha_siz);
    alpha_tilde(T_tilde) = x_LS(1:length(T_tilde));
    x2_tilde = reshape(x_LS(length(T_tilde)+1:end),rows,cols);
    
    % update support & cosupport
    [~,indices] = sort(abs(alpha_tilde(:)),'descend');
    T = indices(1:s);
    Omega_x2_tilde = Omega(x2_tilde);
    [~,indices] = sort(abs(Omega_x2_tilde(:)),'ascend');
    Lambda = indices(1:l);
    
    % calculate new estimates
    alpha = zeros(alpha_siz);
    alpha(T) = alpha_tilde(T);
    x_sparse = D(alpha); % simple thresholding for x_sparse
    % for x_cosparse projection onto null(OmegaL):
    if length(Lambda) < rows*cols
        OmegaL_x2 = zeros(rowsAnalysis,colsAnalysis);
        OmegaL_x2(Lambda) = Omega_x2_tilde(Lambda);
        [x2_proj, res_dagger] = dagger(Lambda,Omega,Omegat,OmegaL_x2);
        x_cosparse = x2_tilde - x2_proj;
    else
        [x_cosparse,res_dagger] = cosparse_projection(x2_tilde, Omega, Omegat, Lambda);
    end

    % calculate residual
    r = y - A(x_sparse + x_cosparse);
    
    errTemp = norm(r)/norm(y);
    fprintf('SACoSaMP Iter: %d, residual error: %f \n', t, errTemp);
    
    if ((err-errTemp)/err <= tol_early) %Early termination condition
        counter = counter+1;
    else
        counter = 0;
    end
    if errTemp < best_err
        best_err = errTemp;
        best_iter = t;
        best_alpha_hat = alpha;
        best_X1 = x_sparse;
        best_X2 = x_cosparse;
    end
    err = errTemp;
    errTrace = [errTrace,err];
    if counter==max_non_improvement_counter
        fprintf('Terminating.\n');
        break;
    end
end

