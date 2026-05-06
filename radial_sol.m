function [D1,Q1,D2,Q2,r] = radial_sol(Z,el)
%RADIAL_SOL Summary of this function goes here
%   Detailed explanation goes here

load('F3.mat');
load('z3.mat');


r = [1e-4:1e-4:12];

b = 1./2.*(3.*pi./4).^(2./3);
x = (Z.^(1./3)./b).*r;

chi = interp1(z3,F3(:,1),x);
dchi = interp1(z3,F3(:,2),x).*(Z.^(1./3)./b);


V = -2.*((Z-1).*chi + 1)./r;
%dV0 = (V(2) - V(1))./dr;
%dVn = (V(end) - V(end-1))./dr;
%dV = [dV0,(V(3:end) - V(1:end-2))./(2.*dr),dVn];
dVth = 2.*((Z-1).*chi+1)./r.^2 - 2.*(Z-1).*dchi./r;

alf = 1./(137.035999);

if (el~=0)
    ls = -(el+1)./2;
else
    ls = 0;
end

Hd = el.*(el+1)./(r.^2) + V + alf.^2.*dVth.*ls./(2.*r);

H = hamiltonian(r,Hd);

E0 = -Z.^2;

[D1] = bisection_eig_sym(H,E0,0,1e-6);

Q1 = zeros(length(r),length(D1));
for ii = 1 : length(D1)
    [q,d] = my_inv_pw_eig(H,D1(ii));
    Q1(:,ii) = q;
    D1(ii) = d;
end

ls = (el-1./2)./2;

Hd = el.*(el+1)./(r.^2) + V + alf.^2.*dVth.*ls./2;

H = hamiltonian(r,Hd);

E0 = -Z.^2;

[D2] = bisection_eig_sym(H,E0,0,1e-6);

Q2 = zeros(length(r),length(D2));
for ii = 1 : length(D2)
    [q,d] = my_inv_pw_eig(H,D2(ii));
    Q2(:,ii) = q;
    D2(ii) = d;
end



end



function [H] = hamiltonian(xx,diag)
%

%   

N = length(xx);

dx = xx(2) - xx(1);

ii = [1:N];
iix = [1:N-1];
jjx = [2:N];


t = ones(1,N).*(2./(dx.^2)) + diag;
a = ones(size(iix)).*(-1./(dx.^2));

jj_ = [ii,jjx,iix];
ii_ = [ii,iix,jjx];
t  = [t,a,a];

H = sparse(ii_,jj_,t,N,N);

end



