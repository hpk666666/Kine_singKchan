function Pp=Ve_int(u_e,Vt,Ve,data,n)
P=zeros(1,5);
for i=1:n+1
    P(1)=1/(i+u_e+4)*data(1,i)*(Vt-Ve)^(i+u_e+4)+P(1);
    P(2)=1/(i+u_e+3)*data(2,i)*(Vt-Ve)^(i+u_e+3)+P(2);
    P(3)=1/(i+u_e+2)*data(3,i)*(Vt-Ve)^(i+u_e+2)+P(3);
    P(4)=1/(i+u_e+1)*data(4,i)*(Vt-Ve)^(i+u_e+1)+P(4);
    P(5)=1/(i+u_e)*data(5,i)*(Vt-Ve)^(i+u_e)+P(5);
end
Pp=P(1)+P(2)+P(3)+P(4)+P(5);
end