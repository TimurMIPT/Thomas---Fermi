function [F,z] = ThF_ptnl(x,C)
%THF_PTNL solver for Chi function of Thomas-Fermi potential

%x = [0:0.0001:60];
%[F3,z3] = ThF_ptnl(x,-1.5878657608815);


t = sqrt(x);
y0 = [1,C];

[z,y] = ode45(@fun,t,y0,t);

F = y;
z = z.^2;

end


function df = fun(t,y)
    df(1,1) =  2.*y(2).*t;
    df(2,1) =  2.*sqrt((y(1).^3));
end

