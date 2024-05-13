function [Pi]=DPi(dP0,dP1,dP2,dP3,dP4,dX)
Pi=zeros(1,5);
Pi(1,1)=trapz(dX,dP0);
Pi(1,2)=trapz(dX,dP1);
Pi(1,3)=trapz(dX,dP2);
Pi(1,4)=trapz(dX,dP3);
Pi(1,5)=trapz(dX,dP4);
end