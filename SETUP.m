h = 1;
T = 200;
delt = 1/5;
J = round((b-a)/h);
dimJ = J+1;
n = (dimJ)^2;
N = round(T/delt);

% Set up the domain
[XX,YY] = meshgrid(1:dimJ,1:dimJ);
xx = XX(:); yy = YY(:);

L = sparse(n,n);
% Construct matrix L (without 1/h^2 factor)
L(1,1) = 3; L(1,2) = -3/2; L(J+1,J+1) = 6; L(J+1,J) = -3;
L = L+sparse(2:J,3:J+1,-1,n,n);
L = L+sparse(2:J,2:J,4,n,n);
L = L+sparse(2:J,1:J-1,-1,n,n);
L(1,J+2) = -3/2; L(J+1,2*J+2) = -3;
L = L+sparse(2:J,J+3:2*J+1,-2,n,n);
L(n-J,n-J) = 6; L(n-J,n-J+1) = -3;
L(n,n) = 3; L(n,n-1) = -3/2;
L = L+sparse(n-J+1:n-1,n-J+2:n,-1,n,n);
L = L+sparse(n-J+1:n-1,n-J+1:n-1,4,n,n);
L = L+sparse(n-J+1:n-1,n-J:n-2,-1,n,n);
L(n-J,n-(2*J+1)) = -3; L(n,n-dimJ) = -3/2;
L = L+sparse(n-J+1:n-1,n-2*J:n-(J+2),-2,n,n);
L = L+sparse(J+2:n-dimJ,2*J+3:n,-1,n,n);
L = L+sparse(J+2:n-dimJ,1:n-2*dimJ,-1,n,n);
L = L+sparse(J+2:n-dimJ,J+2:n-dimJ,4,n,n);
L = L+sparse(J+2:n-(J+2),J+3:n-dimJ,-1,n,n);
L = L+sparse(J+2:dimJ:n-(2*J+1),J+3:dimJ:n-2*J,-1,n,n);
L = L+sparse(2*J+2:dimJ:n-2*dimJ,2*J+3:dimJ:n-(2*J+1),1,n,n);
L = L+sparse(J+3:n-dimJ,J+2:n-(J+2),-1,n,n);
L = L+sparse(2*J+2:dimJ:n-dimJ,2*J+1:dimJ:n-(J+2),-1,n,n);
L = L+sparse(2*J+3:dimJ:n-(2*J+1),2*J+2:dimJ:n-2*dimJ,1,n,n);


u = zeros(n,NumSpp);
F = zeros(n,NumSpp);
U_grid = zeros(dimJ,dimJ,NumSpp);

CompHeight = 0.5;
