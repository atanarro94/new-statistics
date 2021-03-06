function [y,alpha] = hill_values(px)
xf = px*54/1.93;
x1 = xf;
y = zeros(1,length(xf));
ymax = 28;

%% Values of y
%Compute points first hill
for i=1:length(px)
    if xf(i)>=54.0/1.93*7.07
        x1(i)=9*54/1.93-x1(i);
    end
    if x1(i)<=9
        y(i) = min([ymax, 28 + 6.77507969851e-03*x1(i)^2 - 2.1245277758e-03*x1(i)^3]);
    elseif x1(i)>9 && x1(i)<=14
        y(i) = 25.07355893131+0.9754803562315*x1(i)-1.016116352781e-01*x1(i)^2+1.889794677828e-03*x1(i)^3;
    elseif x1(i)>14 && x1(i)<= 20
        y(i) = 2.579601052357E+01+8.206693007457E-01*x1(i)-9.055370274339e-02*x1(i)^2+1.626510569859e-03*x1(i)^3;
    elseif x1(i)>20 && x1(i)<=30
        y(i) = 4.046435022819E+01-1.379581654948E+00*x1(i)+1.945884504128E-02*x1(i)^2-2.070318932190E-04*x1(i)^3;
    elseif x1(i)>30 && x1(i)<=40
        y(i) = 1.792461334664E+01+8.743920332081E-01*x1(i)-5.567361123058E-02*x1(i)^2+6.277731764683E-04*x1(i)^3;
    elseif x1(i)>40 && x1(i)<=54 
        y(i) = max([0,5.639011190988E+01-2.010520359035E+00*x1(i)+1.644919857549E-02*x1(i)^2+2.674976141766E-05*x1(i)^3]);
    end
end
y = y/ymax;

%% Compute derivatives
for i=1:length(px)
    if xf(i)>=54.0/1.93*7.07
        x1(i)=9*54/1.93-xf(i);
        ft = -1;
    else
        ft = 1;
    end 
    if x1(i)<=9
        dydx(i) = ft*min([0,2*6.77507969851e-03*x1(i) - 3*2.1245277758e-03*x1(i)^2]);
    elseif x1(i)>9 && x1(i)<=14
        dydx(i) = ft*(0.9754803562315-2*1.016116352781e-01*x1(i)+3*1.889794677828e-03*x1(i)^2);
    elseif x1(i)>14 && x1(i)<= 20
        dydx(i) = ft*(8.206693007457E-01-2*9.055370274339e-02*x1(i)+3*1.626510569859e-03*x1(i)^2);
    elseif x1(i)>20 && x1(i)<=30
        dydx(i) = ft*(-1.379581654948E+00+2*1.945884504128E-02*x1(i)-3*2.070318932190E-04*x1(i)^2);
    elseif x1(i)>30 && x1(i)<=40
        dydx(i) = ft*(8.743920332081E-01-2*5.567361123058E-02*x1(i)+3*6.277731764683E-04*x1(i)^2);
    elseif x1(i)>40 && x1(i)<=54 
        dydx(i) = ft*(-2.010520359035E+00+2*1.644919857549E-02*x1(i)+3*2.674976141766E-05*x1(i)^2);
    else 
        dydx(i) = 0;
    end
end
alpha = atan(dydx*54/1.93/28);
