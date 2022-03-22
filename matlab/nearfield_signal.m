function [x,sig]=nearfield_signal(N);
% examples de champs (N= 1 ou 2, ...)
% sig: near_field signal, x: absica
% N: number of source: 1 or 2

u=1;             % length unit, e.g. 1 µm
x=[-5:0.01:5]*u; % signal zone

if N==1
sig=dirac(x); 
% dirac(vector): arrays of zero except x=0 =>infinity
elseif N==2
sig=dirac(x-1) +dirac(x+1);

elseif N==-1
    sig=dirac(x-x(1));
elseif N==0
    sig=dirac(x-x(end));
    
elseif N==10
    sig=cos(2*pi*(x-x(1)).^2);
elseif N==11
    sig=cos(2*pi*(x.^2));
elseif N==12
    sig=cos(2*pi*(x-x(end)).^2);
elsif N==20
    sig=cos(2*pi*(x-x(end)));
end


dirac_position=(sig==Inf);

sig(dirac_position)=1;

%plot(x,sig)

