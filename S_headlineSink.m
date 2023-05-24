function [F_Head_line_Sink]=S_headlineSink(z1,z2,z)
zc=0.5*(z2+z1);
zL=0.5*(z2-z1);
zm=(z-zc)/zL;
F_Head_line_Sink=real((abs(z2-z1)/(4*pi))*((zm+1)*log(zm+1)-(zm-1)*log(zm-1)+2*log(zL)-2));
end
%S_headlineSink(0+0*i,0+1500*i,1500+2000*i)
%S_headlineSink(0+0*i,0+1500*i,0+750*i)


