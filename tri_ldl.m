function [l,d] = tri_ldl(A,z)
%TRI_LDL LDL разложение для трёхдиагональной матрицы симметричной A-I*z
%   

szA = size(A);

a = full(diag(A));
b = full(diag(A,1));

a = a- z;
d = a;
d(:) = 0;
l = b;
l(:)= 0;

d(1) = a(1);
l(1) = b(1)./d(1);

for ii = 2:szA(1)-1
    d(ii) = a(ii)-l(ii-1).*b(ii-1);%(b(ii-1).^2)./d(ii-1);
    l(ii) = b(ii)./d(ii);
end
ii = szA(1);
d(ii) = a(ii)-l(ii-1).*b(ii-1);

end

