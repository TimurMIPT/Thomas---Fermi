function [x] = tri_slau_solver(l,d,b)
%TRI_SLAU_SOLVER решение СЛАУ A*x = b;
% А = A';  
% d - главная диагональ матрицы D из разложения L*D*L' = A, 
% l - нижняя диагональ матрицы L из разложения L*D*L' = A
% L*D*L' * x = b
% L*z = b
% D*y = z
% L'*x = y

    ld = length(d);
    z = b;
    for ii = 2:ld % цикл прогонки вниз
       z(ii) = b(ii) - l(ii-1).*z(ii-1);
    end
    z = z./d;
    x = zeros(ld,1);
    x(ld) = z(ld);
    for ii = ld-1:-1:1 % цикл прогонки вверх
        x(ii) = z(ii)- l(ii).*x(ii+1);
    end
end