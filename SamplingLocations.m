function Loc = SamplingLocations(n, k)
% returns sampling locations
% corresponding to Fourier coefficients of n x n image
% along k radial lines.

Loc = zeros(n, n);  
alpha = pi/k;
for j = 1:k
    tana = tan(alpha*j);
    if abs(tana) <= 1
        for x = 1:n
            y = floor(tana*(-(n+1)/2 + x) + (n+1)/2 + 0.5);
            Loc(x, y) = 1;
        end
    else
        for y = 1:n
            x = floor((-(n+1)/2 + y)/tana + (n+1)/2 + 0.5);
            Loc(x, y) = 1;
        end
    end
end

Loc = ifftshift(Loc);
