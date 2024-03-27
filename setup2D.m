clear all
clc
format long g
x = 0:.01:2*pi;
y = 0:.01:2*pi;

v = [-sin(y), sin(x), 0];
B = [-sin(y), sin(2*x), 0];

kin_energy = 0;
mag_energy = 0;
for i=1:length(x)
    for j=1:length(y)
        ux = -2*sin(j);
        uy = 2*sin(i);
        u(i,j) = sqrt(ux*ux+uy*uy);
        bx = -2*sin(j);
        by = 2*sin(2*i);
        b(i,j) = sqrt(bx*bx+by*by);
        kin_energy =  kin_energy+u(i,j)*u(i,j);
        mag_energy =  mag_energy+b(i,j)*b(i,j);
    end
end
kin_energy_density = 0.5*kin_energy/(length(x))^2;
mag_energy_density = 0.5*mag_energy/(length(x))^2;

U_ref = rms(rms(u));
B_ref = rms(rms(b));

T = 2;
D = 100;
St = 1E-3;
Sl = 2*pi/D;
Sv = Sl/St;
U_reflb = U_ref/Sv;
nstep = T/St
Mach = U_reflb*sqrt(3)
V = St*U_ref/Sl
B = St*B_ref/Sl;

Re = 628;
ni = U_ref*2*pi./Re;
%ni = 0.01;
nilb = ni/Sl/Sl*St;
Relb = U_reflb*D./nilb;
tau = nilb*3+.5
taum = nilb*3+.5;


Scale_energy = Sl*Sl/(St*St);
Scale_enstrophy = St*St;



