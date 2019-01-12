function [z,res] = dagger(L,Omega,Omegat,b)

A_LS = @(z) OmegaL_OmegaLt(L,Omega,Omegat,z,size(b));
[bTemp, res, iter] = cgsolve(A_LS, b(L), 1e-7, 1000, 0);
zTemp = zeros(size(b));
zTemp(L) = bTemp;
z = Omegat(zTemp);

end

function z = OmegaL_OmegaLt(L,Omega,Omegat,bTemp,SIZ)

temp1 = zeros(SIZ);
temp1(L) = bTemp;
temp2 = Omegat(temp1);
temp3 = Omega(temp2);
z = temp3(L);
end

