function Pp=VL_int(u,Vt,VL,data,n)
P=zeros(1,5);
for i=1:n+1
     P(1)=1/(i+u)*data(1,i)*(VL-Vt)^(i+u)+P(1);
     P(2)=1/(i+u)*data(2,i)*(VL-Vt)^(i+u)+P(2);
     P(3)=1/(i+u)*data(3,i)*(VL-Vt)^(i+u)+P(3);
     P(4)=1/(i+u)*data(4,i)*(VL-Vt)^(i+u)+P(4);
     P(5)=1/(i+u+1)*data(5,i)*(VL-Vt)^(i+u+1)+P(5);
end
Pp=P(1)+P(2)+P(3)+P(4)+P(5);
end