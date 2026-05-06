function [F,z] = ThF_potential(x,C)
%THF_POTENTIAL solver for Chi function of Thomas-Fermi potential

%x = [1e-10,0.0001:0.0001:40];
%[F,z] = ThF_potential(x,-1.616624655);

C1 = C + 2.*sqrt(x(1));

y0 = [1+C1.*x(1),C1];


[z,y] = ode45(@fun,x,y0,[x(1),x(end)]);

F = y;

end


function df = fun(x,y)
    df(1,1) =  y(2);
    df(2,1) =  sqrt((y(1).^3)./x);
end
