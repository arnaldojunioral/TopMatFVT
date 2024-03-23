% COMPUTE FILTER
function Hs = Filtering(nx,ny,rfil)

iH = ones(nx*ny*(2*(ceil(rfil)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;

% Compute filter values
for j1 = 1:ny
    for i1 = 1:nx
        e1 = (j1 - 1) * nx + i1;
        for j2 = max(j1-(ceil(rfil)-1),1):min(j1+(ceil(rfil)-1),ny)
            for i2 = max(i1-(ceil(rfil)-1), 1):min(i1+(ceil(rfil)-1),nx)
                e2 = (j2 - 1) * nx + i2;
                k = k + 1;
                iH(k) = e1;
                jH(k) = e2;
                sH(k) = max(0,rfil-sqrt((i1-i2)^2+(j1-j2)^2));
            end
        end
    end
end
% Create sparse matrix and compute row sums
H = sparse(iH,jH,sH);
Hs = H./sum(H,2);

%% This Matlab function is defined in topX code (Xia and Breitkopf, 2015)