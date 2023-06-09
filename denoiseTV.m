function [u,error] = denoiseTV(f, pm)
% |x|+|y| + 0.5*lambda*||Dx u-x+bx||^2 + 0.5*lambda*||Dy u-y+by||^2 + 0.5*mu*||u-f||^2
%


[rows,cols] = size(f);

mu = 20; lambda = 1; nIter = 10000; tol = 1e-5;
u_orig = zeros(rows, cols);

if isfield(pm,'lambda'); lambda = pm.lambda; end
if isfield(pm,'mu'); mu = pm.mu; end
if isfield(pm,'nIter'); nIter = pm.nIter; end
if isfield(pm,'u_orig'); I = pm.u_orig; end
if isfield(pm,'tol'); tol = pm.tol; end; % stop condition



% Reserve memory for the auxillary variables
f0 = f;
u = zeros(rows,cols);
x = zeros(rows,cols);
y = zeros(rows,cols);
bx = zeros(rows,cols);
by = zeros(rows,cols);

% Build Kernels
uker = zeros(rows,cols);
uker(1,1) = 4;uker(1,2)=-1;uker(2,1)=-1;uker(rows,1)=-1;uker(1,cols)=-1;
uker = mu+lambda*fft2(uker);   


%  Do the reconstruction
for inner = 1:nIter
    % update u
    uold = u;
    rhs = mu*f+lambda*Dxt(x-bx)+lambda*Dyt(y-by); 
    u = real(ifft2(fft2(rhs)./uker));
    
    % update x and y
    dx = Dx(u);
    dy = Dy(u);
    
    
    %% isotropic TV
    % [x,y] = shrink2( dx+bx, dy+by,1/lambda);
    % anisotropic TV
    x = shrink(dx+bx, 1/lambda);
    y = shrink(dy+by, 1/lambda);
    
    % update bregman parameters
    bx = bx+dx-x;
    by = by+dy-y;
    % 
    if (norm(uold-u)/norm(u)<tol)
        break
    end
    
    error(inner) = norm(uold-u,'fro')/norm(uold,'fro');
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

function [xs,ys] = shrink2(x,y,lambda)

s = sqrt(x.*conj(x)+y.*conj(y));
ss = s-lambda;
ss = ss.*(ss>0);

s = s+(s<lambda);
ss = ss./s;

xs = ss.*x;
ys = ss.*y;

return;

function z = shrink(x,r)
z = sign(x).*max(abs(x)-r,0);
return;

