function [F] = ThF_pt(x)
%THF_PT Summary of this function goes here
%   Detailed explanation goes here



D = laplacian_1D(x);


ye = exp(-x);

f0 = 0.059.*x.*exp(-x);

f1 = f0;
for ii = 1 : 10
    df0 = sqrt(((f1 + ye).^3)./x) - ye;
    f1 = D\df0;
    ff(:,ii) = f1 + ye;
end



end


function [D] = laplacian_1D(xx)
%LAPLACIAN_1D генератор оператора лапласиана для 1D сетки

%   

N = length(xx);


dx = xx(2) - xx(1);



ii = [1:N];
iix = [1:N-1];
jjx = [2:N];


t = ones(1,N).*(2./dx^2);
a = ones(size(iix)).*(-1./dx^2);

jj = [ii,jjx,iix];
ii = [ii,iix,jjx];
t  = -[t,a,a];

D = sparse(ii,jj,t,N,N);

end
