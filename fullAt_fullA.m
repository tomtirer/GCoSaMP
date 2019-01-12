function z = fullAt_fullA(fullx,A,At,D,Dt,Omega,Omegat,T,L,rows,cols,alpha_siz,cosparse_constraint_penalty_for_LS)

alpha = fullx(1:length(T));
x2 = reshape(fullx(length(T)+1:end),rows,cols);

alpha_full = zeros(alpha_siz);
alpha_full(T) = alpha;

temp1 = D(alpha_full);
temp2 = At(A(temp1)) + At(A(x2));
temp3 = Dt(temp2);
upperPart = temp3(T);

temp4 = Omega(x2);
temp5 = zeros(size(temp4));
temp5(L) = temp4(L);
temp6 = temp2 + cosparse_constraint_penalty_for_LS * Omegat(temp5);
lowerPart = temp6(:);

z = [upperPart;lowerPart];

