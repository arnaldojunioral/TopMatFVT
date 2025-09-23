% TOPOLOGY OPTIMIZATION OF PERIODIC MATERIAL BY STANDARD FVT FORMULATION
function TopMatFVT(nx,ny,volfrac,penal,rfil,ft)

%______________________________________________________________USER-DEFINED
E0    = 1.0;                                                               % Young's modulus of solid material
nu    = 0.3;                                                               % Poisson's ratio
ctp   = 3;                                                                 % obj. func: 1 (shear modulus), 2 (bulk modulus) and 3 (Poisson's ratio)
R     = min(nx,ny)/6;                                                      % radius of circular material heterogeneity
mdl   = 'RAMP';                                                            % 'SIMP' or 'RAMP' method
eta   = 1/3;                                                               % damping factor
move  = 0.1;                                                               % move limit

%______________________________________________PERIODIC BOUNDARY CONDITIONS

% Subvolume faces
[i,j] = meshgrid(1:nx,1:ny);
q = (i+(j-1)*nx)';
faces = [q(:),q(:)+nx*ny+1,q(:)+nx,q(:)+nx*ny];
faces(end-nx+1:end,3) = faces(1:nx,1);
faces(nx:nx:end,2) = faces(1:nx:end-nx+1,4);

% Degrees of freedom
dof = zeros(nx*ny,8);
dof(:,2:2:end) = 2*faces;
dof(:,1:2:end) = 2*faces-1;
ndof = max(dof(:));

% Degrees of freedom: fixed and free
fixed = [dof(1,1:2),dof(end,3)];
free = setdiff(dof(:),fixed);

% Sparse mapping indices
iK = reshape(kron(dof,ones(8,1))',64*nx*ny,1);
jK = reshape(kron(dof,ones(1,8))',64*nx*ny,1);
iF = repmat((dof)',3,1);
jF = [ones(8,nx*ny);2*ones(8,nx*ny);3*ones(8,nx*ny)];

%____________________________________LOCAL STIFFNESS MATRIX AND LOAD VECTOR

% Constitutive matrix for unit Young's modulus
C0 = 1/(1-nu^2)*[1,nu,0;nu,1,0;0,0,(1-nu)/2];

% Normal vectors
N1 = [0,0,-1;0,-1 0]; N2 = [1,0,0;0,0,1];
N3 = [0,0,1;0,1,0];   N4 = [-1,0,0;0,0,-1];
N = [N1,zeros(2,9);zeros(2,3),N2,zeros(2,6);
    zeros(2,6),N3,zeros(2,3);zeros(2,9),N4];

% Kinematic matrices:
a = repmat(eye(2), 4, 1);
invA = [0,0,1,0,0,0,-1,0;-1,0,0,0,1,0,0,0;
    0,0,2,0,0,0,2,0;2,0,0,0,2,0,0,0;
    0,0,0,1,0,0,0,-1;0,-1,0,0,0,1,0,0;
    0,0,0,2,0,0,0,2;0,2,0,0,0,2,0,0];

% Strain/Displacement operator per local_faces: strain = E * Wi
E = [1,0,0,0,0,0,0,0;0,0,0,0,0,1,0,-3/2;0,1,0,-3/2,1,0,0,0;
    1,0,3/2,0,0,0,0,0;0,0,0,0,0,1,0,0;0,1,0,0,1,0,3/2,0;
    1,0,0,0,0,0,0,0;0,0,0,0,0,1,0,3/2;0,1,0,3/2,1,0,0,0;
    1,0,-3/2,0,0,0,0,0;0,0,0,0,0,1,0,0;0,1,0,0,1,0,-3/2,0];

% Operator that relates the local surfaces-averaged tractions x Wi
B = N*[C0,zeros(3,9);zeros(3,3),C0,zeros(3,6);
    zeros(3,6),C0,zeros(3,3);zeros(3,9),C0]*E;
sumB = (B(1:2,:)+B(5:6,:))+(B(3:4,:)+B(7:8,:));

% Matrices a_barra and A_barra
ab = (sumB*invA*a)\(sumB*invA);
Ab = invA*(eye(8)-a*ab);

% Local stiffness matrix
K0 = B*Ab;

% Local load vector
H0 = [N1;N2;N3;N4]*C0;

%_____________________________________________________MATERIAL PENALIZATION
Emin  = E0*(1e-9);                                                         % soft material to avoid singularity
x = InitialMaterialDesign(nx,ny,R,volfrac);
if (strcmp(mdl,'SIMP'))
    MatInt = @(p,x) {Emin+(E0-Emin)*x.^p, p*(E0-Emin)*x.^(p-1)};
elseif (strcmp(mdl,'RAMP'))
    MatInt = @(p,x) {Emin+(E0-Emin)*x./(1+p*(1-x)),...
        (1+p)*(E0-Emin)./(1+p*(1-x)).^2};
end

%_________________________________________________________________FILTERING
if ft == 1                                                                 % sensitivity filter
    Hs = Filtering(nx,ny,rfil);
    sensitivity = @(x,dfdx,dvdx){Hs*(x(:).*dfdx(:))./max(1e-3,x(:)),dvdx};
elseif ft == 2                                                             % density filter
    Hs = Filtering(nx,ny,rfil);
    sensitivity = @(x,dfdx,dvdx){Hs*(dfdx(:)),Hs*(dvdx(:))};
else
    Hs = [];
    sensitivity = @(x,dfdx,dvdx){dfdx,dvdx};                               % no-filtering
end

%_____________________________________________________PREALLOCATE VARIABLES
uf = zeros(ndof,3);                                                        % surface-averaged fluctuating displacements
u = zeros(nx*ny,8,3);                                                      % total surface-averaged displacements (macro + fluctuating)
C = zeros(3);                                                              % effective constitutive matrix
dCdx = cell(3);                                                            % sensitivity of local constitutive matrix
dvdx = ones(nx,ny);                                                        % sensitivity of the volume

% Compute macroscopic displacements
u0 = cell(nx*ny,1);
for j = 1:ny
    for i = 1:nx, q = i+(j-1)*nx;
        u0{q} = [1/2+(i-1),0,1/2*(j-1);0,(j-1),1/4+1/2*(i-1);   ...
            1+(i-1),0,1/4+1/2*(j-1);0,1/2+1*(j-1),1/2+1/2*(i-1);...
            1/2+1*(i-1),0,1/2+1/2*(j-1);0,1+(j-1),1/4+1/2*(i-1);...
            (i-1),0,1/4+1/2*(j-1);0,1/2+(j-1),1/2*(i-1)];
    end
end

%______________________________________________________OPTIMIZATION PROCESS
itemax = 0;
for p = 1:length(penal(:))                                                 % continuation method

    % Convergency criterium and iteration number
    change = 1.0; iter = 0;

    % Actualize design variable
    xPhys = x;

    % Print current penalty factor
    fprintf('\nPenalty factor: %1.2f\n',penal(p));

    % Iterative process
    while (change > 0.01), iter = iter+1;

        % MATERIAL PROPERTIES
        Mat = MatInt(penal(p),xPhys);
        E = Mat{1}; dEdx = Mat{2};

        % FINITE-VOLUME THEORY ANALYSIS

        % Material Interpolation
        sK = K0(:)*E(:)';
        sF = H0(:)*E(:)';

        % Assembly global stiffness matrix
        K = sparse(iK(:),jK(:),sK(:),ndof,ndof); K = (K+K')/2;

        % Assembly load vectors corresponding for three unit strain test
        F = -sparse(iF(:),jF(:),sF(:),ndof,3);

        % Compute fluctuating displacements for three unit strain test
        uf(free,:) = K(free,free)\F(free,:);

        % MATERIALS' HOMOGENIZATION

        % Compute total displacements
        for q = 1:(nx*ny)
            u(q,:,:) = u0{q}+uf(dof(q,:),:);
        end
        for i = 1:3, ui = u(:,:,i);
            for j = 1:3, uj = u(:,:,j);
                sumE = reshape(sum((ui*K0).*uj,2),nx,ny)/(nx*ny);
                C(i,j) = sum(sum(E.*sumE));
                dCdx{i,j} = dEdx.*sumE;
            end
        end

        % OBJECTIVE FUNCTION AND SENSITIVITY
        switch(ctp)
            case 1 % Shear modulus maximization
                f    = C(3,3);
                dfdx = dCdx{3,3};
            case 2 % Bulk modulus maximization
                f    = 1/4*(C(1,1)+C(2,2)+C(1,2)+C(2,1));
                dfdx = 1/4*(dCdx{1,1}+dCdx{2,2}+dCdx{1,2}+dCdx{2,1});
            case 3 % Poisson's ratio minimization
                f    = -(C(1,2)-(0.8^iter)*(C(1,1)+C(2,2)));
                dfdx = -(dCdx{1,2}-(0.8^iter)*(dCdx{1,1}+dCdx{2,2}));
        end

        % FILTER/MODIFICATION OF SENSITIVITIES
        sens = sensitivity(x,dfdx,dvdx);
        dfdx(:) = sens{1}; dvdx(:) = sens{2};

        % UPDATE OF DESIGN VARIABLES - OC METHOD
        xold = x; xlow = xold-move; xupp = xold+move;
        l1 = 0; l2 = 1e9;
        frac = max(0,dfdx./dvdx);
        while (l2-l1 > 1e-9)
            lmid = (l1+l2)/2;
            x = max(0,max(xlow,min(1,min(xupp,xold.*(frac/lmid).^eta))));
            if (mean(x(:)) > volfrac), l1 = lmid; else, l2 = lmid; end
        end
        if (ft == 2), xPhys = reshape(Hs*x(:),size(x));
        else, xPhys = x;
        end
        change = max(abs(x(:)-xold(:)));

        % PRINT RESULTS
        fprintf('It: %i\tObjec.: %1.4f\tVol.: %1.3f\tChange: %1.3f\n',...
            iter,f,mean(xPhys(:)),change);
    end
    itemax = itemax+iter;
end
fprintf('Objec.: %1.4f\n', f);
fprintf('\n'); disp('Optimized Homogenized Constitutive Matrix'); disp(C);
fprintf('itemax: %i\n',itemax);

% PLOT DESIGN
figure(1);
colormap(gray);imagesc(1-xPhys');clim([0 1]);axis equal;axis off;          % optimized topology
figure(2);
colormap(gray);imagesc(repmat(1-xPhys',3,3));clim([0 1]);axis equal off;   % 3x3 cells

%___________________________________________________INITIAL MATERIAL DESIGN
function x = InitialMaterialDesign(nx,ny,R,volfrac)

x = ones(nx,ny);

for j = 1:ny
    for i = 1:nx
        if sqrt((i-nx/2-0.5)^2+(j-ny/2-0.5)^2)< R
            x(i,j) = 1e-20;
        end
    end
end
% for j = 1:ny
%     for i = 1:nx
%         for k = 0.5
%             for m = [0,1]
%                 if sqrt(((i-0.5)-nx*m)^2+((j-0.5)-ny*k)^2) < R
%                     x(i,j) = 1e-20;
%                 end
%             end
%         end
%     end
% end

% density regularization for constraint solid material volume fraction
x = volfrac*x/mean(x(:));

%____________________________________________________________COMPUTE FILTER
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
