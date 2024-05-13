function Pp=V_int_noise(V_up,V_down,V0,data,n)
P=zeros(1,5);
for i=1:n+1
    P(1)=1/(i)*data(1,i)*((V_up-V0)^(i)-(V_down-V0)^(i))+P(1);
    P(2)=1/(i)*data(2,i)*((V_up-V0)^(i)-(V_down-V0)^(i))+P(2);
    P(3)=1/(i)*data(3,i)*((V_up-V0)^(i)-(V_down-V0)^(i))+P(3);
    P(4)=1/(i)*data(4,i)*((V_up-V0)^(i)-(V_down-V0)^(i))+P(4);
    P(5)=1/(i)*data(5,i)*((V_up-V0)^(i)-(V_down-V0)^(i))+P(5);
end
Pp=P(1)+P(2)+P(3)+P(4)+P(5);
end