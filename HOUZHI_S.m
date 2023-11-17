X1=linspace(Ve,Vt,100);
X2=linspace(Vt,VL_star,300);
X=[X1,X2];
% [Ve_P0,Ve_P1,Ve_P2,Ve_P3,Ve_P4]=leijiaS_Ve(u_Ve,X1,Ve,data,n,K);
% [VL1_P0,VL1_P1,VL1_P2,VL1_P3,VL1_P4]=leijiaS_VL(u_VL(1),X2,VL_star,data_1,n,K);
% [VL2_P0,VL2_P1,VL2_P2,VL2_P3,VL2_P4]=leijiaS_VL(u_VL(2),X2,VL_star,data_2,n,K);
% [VL3_P0,VL3_P1,VL3_P2,VL3_P3,VL3_P4]=leijiaS_VL(u_VL(3),X2,VL_star,data_3,n,K);
% [VL4_P0,VL4_P1,VL4_P2,VL4_P3,VL4_P4]=leijiaS_VL(u_VL(4),X2,VL_star,data_4,n,K);
% P0=[Ve_P0,VL1_P0+VL2_P0+VL3_P0+VL4_P0];
% P1=[Ve_P1,VL1_P1+VL2_P1+VL3_P1+VL4_P1];
% P2=[Ve_P2,VL1_P2+VL2_P2+VL3_P2+VL4_P2];
% P3=[Ve_P3,VL1_P3+VL2_P3+VL3_P3+VL4_P3];
% P4=[Ve_P4,VL1_P4+VL2_P4+VL3_P4+VL4_P4];
%*****************************
[Ve_P0,Ve_P1,Ve_P2,Ve_P3,Ve_P4]=leijia2_Ve(u_Ve,X1,Ve,data,n);
[VL1_P0,VL1_P1,VL1_P2,VL1_P3,VL1_P4]=leijia2_VL(u_VL(1),X2,VL_star,data_1,n);
[VL2_P0,VL2_P1,VL2_P2,VL2_P3,VL2_P4]=leijia2_VL(u_VL(2),X2,VL_star,data_2,n);
[VL3_P0,VL3_P1,VL3_P2,VL3_P3,VL3_P4]=leijia2_VL(u_VL(3),X2,VL_star,data_3,n);
[VL4_P0,VL4_P1,VL4_P2,VL4_P3,VL4_P4]=leijia2_VL(u_VL(4),X2,VL_star,data_4,n);
P0=[K(1)*Ve_P0,K(2)*VL1_P0+K(3)*VL2_P0+K(4)*VL3_P0+K(5)*VL4_P0];
P1=[K(1)*Ve_P1,K(2)*VL1_P1+K(3)*VL2_P1+K(4)*VL3_P1+K(5)*VL4_P1];
P2=[K(1)*Ve_P2,K(2)*VL1_P2+K(3)*VL2_P2+K(4)*VL3_P2+K(5)*VL4_P2];
P3=[K(1)*Ve_P3,K(2)*VL1_P3+K(3)*VL2_P3+K(4)*VL3_P3+K(5)*VL4_P3];
P4=[K(1)*Ve_P4,K(2)*VL1_P4+K(3)*VL2_P4+K(4)*VL3_P4+K(5)*VL4_P4];
% [S0,S1,S2,S3,S4]=jifen(nP0,nP1,nP2,nP3,nP4,X);
% S0=trapz(X(2:length(X)-1),P0(2:length(X)-1));
% S1=trapz((2:length(X)-1)X,(2:length(X)-1)P1);
% S2=trapz(X(2:length(X)-1),P2(2:length(X)-1));
% S3=trapz(X(2:length(X)-1),P3(2:length(X)-1));
% S4=trapz(X(2:length(X)-1),P4(2:length(X)-1));
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
S=-(S0+S1+S2+S3+S4);
[Pi]=DPi(dP0,dP1,dP2,dP3,dP4,dX);%?
[DP0,DP1,DP2,DP3,DP4]=DnewP(dP0,dP1,dP2,dP3,dP4,Pi);%?
D0=trapz(dX,DP0);
D1=trapz(dX,DP1);
D2=trapz(dX,DP2);
D3=trapz(dX,DP3);
D4=trapz(dX,DP4);
% [D0,D1,D2,D3,D4]=jifen(DP1,DP2,DP3,DP4,DP5,dX);
D=-(D0+D1+D2+D3+D4);
subplot(2,3,1);
plot(dX,nP0)
text(Ve,0.2,'\downarrowVe');text(VL_star,0.2,'\downarrowVL')
xlabel('V');ylabel('P(0,V)');
subplot(2,3,2);
plot(dX,nP1)
text(Ve,0.2,'\downarrowVe');text(VL_star,0.2,'\downarrowVL')
xlabel('V');ylabel('P(1,V)');
subplot(2,3,3);
plot(dX,nP2)
text(Ve,0.2,'\downarrowVe');text(VL_star,0.2,'\downarrowVL')
xlabel('V');ylabel('P(2,V)');
subplot(2,3,4);
plot(dX,nP3)
text(Ve,0.2,'\downarrowVe');text(VL_star,0.2,'\downarrowVL')
xlabel('V');ylabel('P(3,V)');
subplot(2,3,5);
plot(dX,nP4)
text(Ve,0.0002,'\downarrowVe');text(VL_star,0.0002,'\downarrowVL')
xlabel('V');ylabel('P(4,V)');
subplot(2,3,6);
plot(dX,nP0+nP1+nP2+nP3+nP4)
text(Ve,0.2,'\downarrowVe');text(VL_star,0.2,'\downarrowVL')
xlabel('V');ylabel('P(V)');