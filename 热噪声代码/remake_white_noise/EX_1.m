function [EX_1]=EX_1(VL,V_mid,Ve,data,n)
EX_1=zeros(1,5);
for i=1:n+1
    EX_1(1,1)=EX_1(1,1)+data(1,i)*1/((i-1)+2)*(VL-V_mid)^((i-1)+2)+data(1,i)*V_mid*1/((i-1)+1)*(VL-V_mid)^((i-1)+1)...
        -data(1,i)*1/((i-1)+2)*(Ve-V_mid)^((i-1)+2)-data(1,i)*V_mid*1/((i-1)+1)*(Ve-V_mid)^((i-1)+1);
    EX_1(1,2)=EX_1(1,2)+data(2,i)*1/((i-1)+2)*(VL-V_mid)^((i-1)+2)+data(2,i)*V_mid*1/((i-1)+1)*(VL-V_mid)^((i-1)+1)...
        -data(2,i)*1/((i-1)+2)*(Ve-V_mid)^((i-1)+2)-data(2,i)*V_mid*1/((i-1)+1)*(Ve-V_mid)^((i-1)+1);
    EX_1(1,3)=EX_1(1,3)+data(3,i)*1/((i-1)+2)*(VL-V_mid)^((i-1)+2)+data(3,i)*V_mid*1/((i-1)+1)*(VL-V_mid)^((i-1)+1)...
        -data(3,i)*1/((i-1)+2)*(Ve-V_mid)^((i-1)+2)-data(3,i)*V_mid*1/((i-1)+1)*(Ve-V_mid)^((i-1)+1);
    EX_1(1,4)=EX_1(1,4)+data(4,i)*1/((i-1)+2)*(VL-V_mid)^((i-1)+2)+data(4,i)*V_mid*1/((i-1)+1)*(VL-V_mid)^((i-1)+1)...
        -data(4,i)*1/((i-1)+2)*(Ve-V_mid)^((i-1)+2)-data(4,i)*V_mid*1/((i-1)+1)*(Ve-V_mid)^((i-1)+1);
    EX_1(1,5)=EX_1(1,5)+data(5,i)*1/((i-1)+2)*(VL-V_mid)^((i-1)+2)+data(5,i)*V_mid*1/((i-1)+1)*(VL-V_mid)^((i-1)+1)...
        -data(5,i)*1/((i-1)+2)*(Ve-V_mid)^((i-1)+2)-data(5,i)*V_mid*1/((i-1)+1)*(Ve-V_mid)^((i-1)+1);
end
end