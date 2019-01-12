% Run script for Synthesis-Analysis CoSaMP (relaxation of the generalized CoSaMP)

% Reference: "Generalizing CoSaMP to Signals from a Union of Low Dimensional Linear Subspaces"
% Authors: Tom Tirer and Raja Giryes
% Journal: Applied and Computational Harmonic Analysis (accepted 2018)

% written in late 2016
% paper's experiments were performed in Matlab2013b (but it works in newer versions)

clear;
addpath('local_dct');

% load X2 - clean image
X2 = double(imread('modified_phantom.png'));
N = size(X2,1); % code assumes square images
n = N^2; % signal dimension
figure; subplot(2,3,1); imshow(uint8(X2)); title('clean image');

% operator for cosparse-analysis sub-model
Omega = @(t) FDAnalysis(t);
Omegat = @(u) FDSynthesis(u);
% dictionary for sparse-synthesis sub-model
D = @(alpha) D_dct(alpha,N);
Dt = @(x) Dt_dct(x,N);

% prepare X1 - structured noise
structured_noise_Ratio = 0.2;
s = 500; % sparsity level
alpha0 = Dt(randn(N,N));
[~,indices] = sort(abs(alpha0(:)),'descend');
alpha0(indices(s+1:end))=0;
%load('alpha0_phantom_paper.mat'); % to use X1 from the paper, else put it comment
X1 = D(alpha0)*structured_noise_Ratio*norm(X2,'fro')/norm(D(alpha0)+eps,'fro');

% create the combined (noisy) signal X1+X2
X = double(uint8(X1 + X2));
subplot(2,3,2); imshow(uint8(X)); title('noisy image');

% create Fourier observations
Loc = SamplingLocations(N, 25);
Loc = logical(Loc);
m = length(find(Loc~=0)); % observation dimension
subplot(2,3,3); imshow(fftshift(Loc)); title(['FFT samples, m/n=' num2str(m/n)]);
A = @(x) FourierSampling(x, Loc);
At = @(y) FourierSamplingH(y, N, Loc);
y = A(X);


naive_rec_X = real(At(y));
naive_rec_X = double(uint8(naive_rec_X));
subplot(2,3,6); imshow(uint8(naive_rec_X)); title(['naive rec., PSNR=' num2str(20*log10(255*sqrt(n)/norm(X-naive_rec_X,'fro')))]);


%% run SACoSaMP

l = 127064; % cosparsity (#zeros) for the analysis sub-model
s = 500; % sparsity (#nonzeros) for the synthesis sub-model

[X1_hat,alpha_hat,X2_hat,errTrace,best_iter] = GCoSaMP_SparseCosparse(y, A, At, D, Dt, Omega, Omegat, l, s, 1e-3, 30, 1e-2);

X_hat = X1_hat+X2_hat;
X_hat = double(uint8(X_hat));
X2_hat = double(uint8(X2_hat));
subplot(2,3,4); imshow(uint8(X2_hat)); title(['clean image rec., PSNR=' num2str(20*log10(255*sqrt(n)/norm(X2-X2_hat,'fro')))]);
subplot(2,3,5); imshow(uint8(X_hat)); title(['noisy image rec., PSNR=' num2str(20*log10(255*sqrt(n)/norm(X-X_hat,'fro')))]);



