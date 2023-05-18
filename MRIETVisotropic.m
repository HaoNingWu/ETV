function [u,err2] = MRIETVisotropic(R,f, pm)
% |x|+|y| - alpha sqrt(x^2+y^2) + 0.5*beta*||Dx u-x+bx||^2 + 0.5*beta*||Dy u-y+by||^2
% s.t. RFu = f
%
% DCA and Split Bregman
%

[rows,cols] = size(f);

mu = 20; beta = 5; alpha = 0.5;
maxit = 1000;
u_orig = zeros(rows, cols);
u0 = zeros(rows,cols); tol = 1e-4;
tau = 1;

if isfield(pm,'mu'); mu = pm.mu; end
if isfield(pm,'beta'); beta = pm.beta; end
if isfield(pm,'alpha'); alpha = pm.alpha; end
if isfield(pm,'maxit'); maxit = pm.maxit; end
if isfield(pm,'u_orig'); u_orig = pm.u_orig; end
if isfield(pm,'u0'); u0 = pm.u0; end
if isfield(pm,'tol'); tol = pm.tol; end; % inner iteration tolerance
if isfield(pm,'tau'); tau = pm.tau; end


u = u0;
ux = Dx(u);
uy  =Dy(u);
ugrad = sqrt(abs(ux).^2+abs(uy).^2);
eps = 1e-8;

%DCA iteration counts
oit = 1;
stop = 0;
kkk = 1;
tstart = tic;


% Build Kernels
scale = sqrt(rows*cols);

uker = zeros(rows,cols);
uker(1,1) = 4;uker(1,2)=-1;uker(2,1)=-1;uker(rows,1)=-1;uker(1,cols)=-1;
uker = mu*(conj(R).*R)+beta*fft2(uker);

x = zeros(rows,cols);
y = zeros(rows,cols);
z = zeros(rows,cols);

ff = f;

F(1) = sum(sum(abs(ux)+abs(uy)-alpha*ugrad)) + mu/2*norm(R.*fft2(u)/scale-f)^2;



uold = ones(rows,cols);



while (oit <= 15 && norm(uold-u,'fro') > 1e-10)
    uold = u;
    
    % Reserve memory for the auxillary variables
    f0 = ff; f = f0; 
    
    bx = zeros(rows,cols);
    by = zeros(rows,cols);
    lambda = zeros(rows,cols);
    
    
    
    %  Do the reconstruction
    for outer = 1:1000

            % update u
            murf = ifft2(mu*R.*(f0+z-lambda))*scale;
            rhs = murf+beta*Dxt(x-bx)+beta*Dyt(y-by);
            u = (ifft2(fft2(rhs)./uker));



            
            % update x and y
            dx = Dx(u);
            dy = Dy(u);
            
            % anisotropic TV
            % x = shrink(dx+bx+alpha*ux/beta, 1/beta);
            % y = shrink(dy+by+alpha*uy/beta, 1/beta);
            s = sqrt(    (abs(dx + bx + alpha*ux/beta)).^2 + (abs(dy + by + alpha*uy/beta)).^2);
            x = max(s-1/beta,0) .* (dx + bx + alpha*ux/beta)./s;
            y = max(s-1/beta,0) .* (dy + by + alpha*uy/beta)./s;
             
            % update z
%            (Caution !) 
%            If the measurements are NOISY, then consider projection when updating z 
%            Uncomment the following and use project_L2 in THE ROXIMITY OPERATOR REPOSITORY

%             if norm((R.*fft2(u)/scale - f0) + lambda) <= tau
%                z = (R.*fft2(u)/scale - f0) + lambda;
%             else
%                z = z;
%             end
%             z = project_L2((R.*fft2(u)/scale - f0) + lambda, tau);
            z = 0;

            % update dual variables
            bx = bx+dx-x;
            by = by+dy-y;
            lambda = lambda + (R.*fft2(u)/scale -f0) - z;
            

        
    end
    
    u = abs(u);
    ux = Dx(u);
    uy = Dy(u);

   err2(oit) = norm(uold-u,'fro');
    
    F(oit+1) = sum(sum(abs(ux)+abs(uy)-alpha*ugrad))+ mu/2*norm(R.*fft2(u)/scale-f)^2;
    
    
    oit = oit + 1;
end


return;


function d = Dx(u)
[rows,cols] = size(u);
d = zeros(rows,cols);
d(:,2:cols) = u(:,2:cols)-u(:,1:cols-1);
d(:,1) = u(:,1)-u(:,cols);
return

function d = Dxt(u)
[rows,cols] = size(u);
d = zeros(rows,cols);
d(:,1:cols-1) = u(:,1:cols-1)-u(:,2:cols);
d(:,cols) = u(:,cols)-u(:,1);
return

function d = Dy(u)
[rows,cols] = size(u);
d = zeros(rows,cols);
d(2:rows,:) = u(2:rows,:)-u(1:rows-1,:);
d(1,:) = u(1,:)-u(rows,:);
return

function d = Dyt(u)
[rows,cols] = size(u);
d = zeros(rows,cols);
d(1:rows-1,:) = u(1:rows-1,:)-u(2:rows,:);
d(rows,:) = u(rows,:)-u(1,:);
return

function z = shrink(x,r)
z = sign(x).*max(abs(x)-r,0);
return;



