ALL_D=cell(1,5);
ALL_S=cell(1,5);
D=zeros(1,n);
S=zeros(1,n);
for k=1:5
for i=1:n
gL=3*M(i);%单位是电导率 S/m ^2
ge=gL+gK;%单位是电导率 S
Cm=0.01*M(i);%单位是电容 F
E=ge/Cm;%单位是 pS/pF
L=gL/Cm;%单位是 pS/pF
[Pe_0,Pe_1,Pe_2,Pe_3,Pe_4]=Pe_Vt(ALLu_Ve{1,k}(i),ALLVt{1,k}(i),ALLVe{1,k}(i),ALLdata{1,k}{1,i},ddai);
[PL1_0,PL1_1,PL1_2,PL1_3,PL1_4]=PL_Vt(ALLu_VL{1,k}{1,i}(1),ALLVt{1,k}(i),ALLVL_star{1,k}(i),ALLdata_1{1,k}{1,i},ddai);
[PL2_0,PL2_1,PL2_2,PL2_3,PL2_4]=PL_Vt(ALLu_VL{1,k}{1,i}(2),ALLVt{1,k}(i),ALLVL_star{1,k}(i),ALLdata_2{1,k}{1,i},ddai);
[PL3_0,PL3_1,PL3_2,PL3_3,PL3_4]=PL_Vt(ALLu_VL{1,k}{1,i}(3),ALLVt{1,k}(i),ALLVL_star{1,k}(i),ALLdata_3{1,k}{1,i},ddai);
[PL4_0,PL4_1,PL4_2,PL4_3,PL4_4]=PL_Vt(ALLu_VL{1,k}{1,i}(4),ALLVt{1,k}(i),ALLVL_star{1,k}(i),ALLdata_4{1,k}{1,i},ddai);
Pe_int=Ve_int(ALLu_Ve{1,k}(i),ALLVt{1,k}(i),ALLVe{1,k}(i),ALLdata{1,k}{1,i},ddai);%这里利用牛顿莱布尼茨公式算(Pe(0,Vt)+...Pe(4,Vt))在【Ve,Vt】的积分
PL_int_1=VL_int(ALLu_VL{1,k}{1,i}(1),ALLVt{1,k}(i),ALLVL_star{1,k}(i),ALLdata_1{1,k}{1,i},ddai);%这里利用牛顿莱布尼茨公式算(PL1(0,Vt)+...PL1(4,Vt))在【Vt,VL】的积分
PL_int_2=VL_int(ALLu_VL{1,k}{1,i}(2),ALLVt{1,k}(i),ALLVL_star{1,k}(i),ALLdata_2{1,k}{1,i},ddai);
PL_int_3=VL_int(ALLu_VL{1,k}{1,i}(3),ALLVt{1,k}(i),ALLVL_star{1,k}(i),ALLdata_3{1,k}{1,i},ddai);
PL_int_4=VL_int(ALLu_VL{1,k}{1,i}(4),ALLVt{1,k}(i),ALLVL_star{1,k}(i),ALLdata_4{1,k}{1,i},ddai);
 A=[ -Pe_0, PL1_0,   PL2_0,   PL3_0,   PL4_0;
     -Pe_1, PL1_1,   PL2_1,   PL3_1,   PL4_1;
     -Pe_2, PL1_2,   PL2_2,   PL3_2,   PL4_2;
     -Pe_3, PL1_3,   PL2_3,   PL3_3,   PL4_3;
     -Pe_4, PL1_4,   PL2_4,   PL3_4,   PL4_4;
     Pe_int,PL_int_1,PL_int_2,PL_int_3,PL_int_4];
 B=[0;0;0;0;0;1];
 format rat
 K=A\B;
X1=linspace(ALLVe{1,k}(i),ALLVt{1,k}(i),100);
X2=linspace(ALLVt{1,k}(i),ALLVL_star{1,k}(i),300);
X=[X1,X2];
[Ve_P0,Ve_P1,Ve_P2,Ve_P3,Ve_P4]=leijia2_Ve(ALLu_Ve{1,k}(i),X1,ALLVe{1,k}(i),ALLdata{1,k}{1,i},ddai);
[VL1_P0,VL1_P1,VL1_P2,VL1_P3,VL1_P4]=leijia2_VL(ALLu_VL{1,k}{1,i}(1),X2,ALLVL_star{1,k}(i),ALLdata_1{1,k}{1,i},ddai);
[VL2_P0,VL2_P1,VL2_P2,VL2_P3,VL2_P4]=leijia2_VL(ALLu_VL{1,k}{1,i}(2),X2,ALLVL_star{1,k}(i),ALLdata_2{1,k}{1,i},ddai);
[VL3_P0,VL3_P1,VL3_P2,VL3_P3,VL3_P4]=leijia2_VL(ALLu_VL{1,k}{1,i}(3),X2,ALLVL_star{1,k}(i),ALLdata_3{1,k}{1,i},ddai);
[VL4_P0,VL4_P1,VL4_P2,VL4_P3,VL4_P4]=leijia2_VL(ALLu_VL{1,k}{1,i}(4),X2,ALLVL_star{1,k}(i),ALLdata_4{1,k}{1,i},ddai);
P0=[K(1)*Ve_P0,K(2)*VL1_P0+K(3)*VL2_P0+K(4)*VL3_P0+K(5)*VL4_P0];
P1=[K(1)*Ve_P1,K(2)*VL1_P1+K(3)*VL2_P1+K(4)*VL3_P1+K(5)*VL4_P1];
P2=[K(1)*Ve_P2,K(2)*VL1_P2+K(3)*VL2_P2+K(4)*VL3_P2+K(5)*VL4_P2];
P3=[K(1)*Ve_P3,K(2)*VL1_P3+K(3)*VL2_P3+K(4)*VL3_P3+K(5)*VL4_P3];
P4=[K(1)*Ve_P4,K(2)*VL1_P4+K(3)*VL2_P4+K(4)*VL3_P4+K(5)*VL4_P4];
%%%%%%%%%%%%%%%%%%%%%%%%%%%
dP0=P0(2:length(P0)-1);
dP1=P1(2:length(P0)-1);
dP2=P2(2:length(P0)-1);
dP3=P3(2:length(P0)-1);
dP4=P4(2:length(P0)-1);
dX=X(2:length(P0)-1);
[nP0,nP1,nP2,nP3,nP4]=newP(dP0,dP1,dP2,dP3,dP4);
S0=trapz(dX,nP0);
S1=trapz(dX,nP1);
S2=trapz(dX,nP2);
S3=trapz(dX,nP3);
S4=trapz(dX,nP4);
S(i)=-(S0+S1+S2+S3+S4);
[Pi]=DPi(dP0,dP1,dP2,dP3,dP4,dX);%?
[DP0,DP1,DP2,DP3,DP4]=DnewP(dP0,dP1,dP2,dP3,dP4,Pi);%?
D0=trapz(dX,DP0);
D1=trapz(dX,DP1);
D2=trapz(dX,DP2);
D3=trapz(dX,DP3);
D4=trapz(dX,DP4);
% [D0,D1,D2,D3,D4]=jifen(DP1,DP2,DP3,DP4,DP5,dX);
D(i)=-(D0+D1+D2+D3+D4);
end
ALL_D{1,k}=D;
ALL_S{1,k}=S;
end