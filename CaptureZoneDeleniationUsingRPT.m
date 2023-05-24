clc
clear all
close all
format short
syms paq;   % Unknown Aquifer Constant
syms phls1;

%% Define Parameters
por=0.35;                  % Porosity of the soil
k=25;                      % Permeability Outside the Traingular Domain         unit--->m/d
H=10;                      % Aquifer Thickness
b=0;                       % Base Elevation
h_ref=25;                  % Reference Point head 
angle=0;                   % Uniform FLow Direction
z_ref=10000+10000*i;         % Reference Point Coordinate
Q=12000;
grad=0.05;
zw=150+150*i;
 zl1=0+0i;
 zl2=1000+0i;
 hl1=25;

%% Equation at Reference Point
z=z_ref;                             
F_Well = Well_Fun(Q,z,zw);
F_Uniform_FLow = Uniform_Flow(k,H,grad,z,angle);
F_Linesink1=S_headlineSink(zl1,zl2,z);
Dis_Pot=Discharge_potential(b,h_ref,k,H);                                    %%  Discharge potential at Reference Point
% paq=Dis_Pot-F_Well-F_Uniform_FLow;
Eq_ref=paq+F_Well+F_Uniform_FLow+phls1*F_Linesink1-Dis_Pot;


%% Equation at the middle of line sink
z=(zl1+zl2)/2;
F_Well = Well_Fun(Q,z,zw);
F_Uniform_FLow = Uniform_Flow(k,H,grad,z,angle);
F_Linesink1=S_headlineSink(zl1,zl2,z);
Dis_Pot=Discharge_potential(b,hl1,k,H);
Eq1=paq+F_Well+F_Uniform_FLow+phls1*F_Linesink1-Dis_Pot;



sol=solve(Eq_ref,Eq1);
paq=double(sol.paq);
phls1=double(sol.phls1);



%% Contour Plot
py=0;                        % Arbitrary variable defined
for y=0:1:300;
    py=py+1;
    px=0;                    % Arbitrary variable defined
    for x=0:1:300;       
        px=px+1;
        z=x+i*y;   
        r=abs(z-zw);
        if r>3
        F_Well = Well_Fun(Q,z,zw);
        F_Uniform_FLow = Uniform_Flow(k,H,grad,z,angle);
        F_Linesink1=S_headlineSink(zl1,zl2,z);
        Dis_Pot=paq+F_Well+F_Uniform_FLow+phls1*F_Linesink1;
        Head(py,px)=Head_Conversion(k,H,Dis_Pot);   
      
        else
         Head(py,px)=NaN; 
        end;
    end       
end
Head;

[cc,hh]=contour(Head);
clabel(cc,hh);
xlabel('x(m)');
ylabel('y(m)');
grid on
title('Plot of Head Contour: Strack plot');

%%.................................................................................................................%%



%% Determination of pore velocity,

%% Method1 (Using Direct Formula):::--->

zb=150+100i;
zt=150+200i;
zl=100+150i;
zr=200+150i;
z=zb;
Q0=k*H*grad;
Z=(z-(0.5*(zl2+zl1)))/(0.5*(zl2-zl1)); 
Wls=(phls1/(2*pi))*log((Z-1)/(Z+1));
Wls_real=real(Wls);
Wls_img=imag(Wls);
Wls=Wls_real+i*Wls_img;
W=-Q0+Wls+(-(Q/(2*pi))*(1/(z-zw))); % Discharge Vector
Qx=real(W);               % X component of discharge vector   
Qy=imag(W);               % Y component of discharge vector
vxd=Qx/(por*H);                  
vyd=Qy/(por*H);
vd=vxd-i*vyd;



%% Method2 (Using hydraulic calculations):::--->
z=z;
F_Well1 = Well_Fun(Q,z,zw);
F_Uniform_FLow = Uniform_Flow(k,H,grad,z,angle);
F_Linesink1=S_headlineSink(zl1,zl2,z);
Dis_Pot1=paq+F_Well1+F_Uniform_FLow+phls1*F_Linesink1;
Head1=Head_Conversion(k,H,Dis_Pot1); 

zx1=z+0.001;
F_Wellx1 = Well_Fun(Q,zx1,zw);
F_Uniform_FLowx1 = Uniform_Flow(k,H,grad,zx1,angle);
F_Linesink1=S_headlineSink(zl1,zl2,zx1);
Dis_Potx1=paq+F_Wellx1+F_Uniform_FLowx1+phls1*F_Linesink1;
Headx1=Head_Conversion(k,H,Dis_Potx1); 
vxh=(-k*(Headx1-Head1)/(zx1-z))/por;

zy1=z+0.001i;
F_Welly1 = Well_Fun(Q,zy1,zw);
F_Uniform_FLowy1 = Uniform_Flow(k,H,grad,zy1,angle);
F_Linesink1=S_headlineSink(zl1,zl2,zy1);
Dis_Poty1=paq+F_Welly1+F_Uniform_FLowy1+phls1*F_Linesink1;
Heady1=Head_Conversion(k,H,Dis_Poty1); 
vyh=(-k*(Heady1-Head1)/(zy1-z))/por;
vh=vxh+vyh;

break_point=1;

%% Particle Tracking
R=20;                       % Say radious of circle in 20 m
delt=0.01;
%% lets create  points over the circle
theta=0:(pi/24):2*pi;
Np=length(theta); 
for ii=1:1:Np
            xp(ii)=R*cos(theta(ii));
            yp(ii)=R*sin(theta(ii));
end
xw=real(zw);
yw=imag(zw);
xw=repmat(xw,1,Np);
yw=repmat(yw,1,Np);
xwp=xw+xp;
ywp=yw+yp;  
zwp=xwp+i*ywp;                                                       %% Points over the point of Circle.
plot(xwp,ywp);                                                       %% Plot
hold on;
grid on;

for kk=1:1:Np          
zwp_C(kk,1)=zwp(kk);                                                 %% zwp converted to Column Vector  %% These are the points over the circle
end

for ii=1:1:Np
      xwp1=xwp(ii);
      ywp1=ywp(ii);
      zwp1=xwp1+i*ywp1;
      zwp_v(ii)=zwp1;
      for tt=1:1:10000   
      vwp1=+(Q0+((Q/(2*pi))*(1/(zwp1-zw))))/(por*H);    %% Porosity also Important Factor
      v_vect(ii,tt)=vwp1;
      vx=(real(vwp1));
      vy=-(imag(vwp1));
      xwp1=xwp1+delt*vx;
      ywp1=ywp1+delt*vy;
      zwp1=xwp1+i*ywp1;
      zwp1_v(ii,tt)=zwp1;
end
end
    
for ii=1:1:Np
      zwp1_v1=[zwp1_v(ii,:)];
      zwp1_v1_real=real(zwp1_v1);
      zwp1_v1_imag=imag(zwp1_v1);
      plot(zwp1_v1_real,zwp1_v1_imag);
      hold on;
end 


break_point=2;











