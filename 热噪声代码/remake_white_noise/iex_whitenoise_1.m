syms x m
M=0.02;%����Ĥ�������λ��um^2%��ʵ��Χ 0.15 
Iex=0 ;%�ⲿ��������λ fA ����   1 A = 10^15 fA ��%��ʵ��Χ 1
%M0.02��0.08|Iex-1��1fA 
%��H_20
ddai=210;
%iex=0,10,-10,ȷ��Ĥ���
%10^-3�ǹ���ʱ��Ļ���ϵ�������L��E����ϵ��һ��
%��������0��ʱ��M����ķ�Χ
gL=3*M;%��λ�ǵ絼�� pS
gK=20;%��λ�ǵ絼�� pS
VL1=-54.4;%��λ�� mV
VK=-77;%��λ�� mV
ge=gL+gK;%��λ�ǵ絼�� pS
Cm=0.01*M;%��λ�ǵ��� pF
E=0.001*ge/Cm;%��λ�� pS/pF
L=0.001*gL/Cm;%��λ�� pS/pF
VL=(gL*VL1+Iex)/gL;%��λ�� mV
Ve=(gL*VL1+gK*VK+Iex)/ge;%��λ�� mV
V_mid=(VL-Ve)/2+Ve;%չ����
Ve_0=Ve;%�߽��
VL_0=VL;%�߽��
K_b=1.380649*10^-23;%������������,��λΪJ/K��JΪ������K������ѧ�¶�
T=294.15;%�����¶ȣ���λK
delta_1=sqrt(2*K_b*T*gL*10^12*10^3)/Cm;%��ʾ��ɢ��ǿ��
delta_2=sqrt(2*K_b*T*ge*10^12*10^3)/Cm;%�������1
o=10;
fa=taylor(0.01*(x+55)/(1-exp(-(x+55)/10)),x,V_mid,'order',o);%�ٶȦ�(v)��V0������̩��չʽ
fb=taylor(0.125*exp(-(x+65)/80),x,V_mid,'order',o);%�ٶȦ�(v)��V0������̩��չʽ
newfa=subs(fa,(x-V_mid),m);%��(x-V0)�滻��m
newfb=subs(fb,(x-V_mid),m);%��(x-V0)�滻��m
Kaa=vpa(coeffs(newfa));%����������̩��չʽ��ϵ��
Lbb=vpa(coeffs(newfb));%����������̩��չʽ��ϵ��
Ka=zeros(1,ddai+1);
Lb=zeros(1,ddai+1);
Ka(1:o)=Kaa;
Lb(1:o)=Lbb;
data=cell(1,10);%����һ��ʮ�е�Ԫ��һ��������һ��������ݼ���ϵ��
intP=cell(1,10);%��һ���õ�ʮ���ұߵ����������Ļ���
P1_V0=cell(1,10);%Ԫ���Ԫ����[P0,P1,P2,P3,P4]�ұ߽߱���ֵ
P2_V0=cell(1,10);%Ԫ���Ԫ����[P0,P1,P2,P3,P4]��߽߱���ֵ
P1_V0_2=cell(1,10);%Ԫ���Ԫ����[P0,P1,P2,P3,P4]�ұ߽߱�㵼����ֵ
P2_V0_2=cell(1,10);%Ԫ���Ԫ����[P0,P1,P2,P3,P4]��߽߱�㵼����ֵ
ZY=eye(10);%����δ֪��
for i=1:10
data{1,i}=V_diedai_noise(ddai,Ve,VL,V_mid,E,Ka,Lb,L,ZY(1:5,i),ZY(6:10,i),delta_1,delta_2);
%ʮ��������              
intP{1,i}=V_int_noise(VL_0,Ve_0,V_mid,data{1,i},ddai);
%��һ���õ���
P1_V0{1,i} = P_zeros_noise(VL_0,V_mid,data{1,i},ddai);%�ұ߽߱���ֵ
P2_V0{1,i} = P_zeros_noise(Ve_0,V_mid,data{1,i},ddai);%��߽߱���ֵ
P1_V0_2{1,i} = P_zeros_noise_2(VL_0,V_mid,data{1,i},ddai);%�ұߵ��������
P2_V0_2{1,i} = P_zeros_noise_2(Ve_0,V_mid,data{1,i},ddai);%��ߵ��������
end
%ȡ-50mv��-70mv
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=zeros(11,10);
for i=1:11
    for j=1:10
        if(i==1)
            A(i,j)=P2_V0_2{1,j}(5,1);%��Ve��n4״̬�ĵ���Ϊ0
        elseif(i<6)
             A(i,j)=P2_V0_2{1,j}(i-1,1)-L*(VL-Ve)*P2_V0{1,j}(i-1,1)*2/delta_1^2;
             %��Ve��n0 1 2 3״̬��Pi��������L*(VL-Ve)*(2/delta_1^2)*Pi
         elseif(i<10)
            A(i,j)=P1_V0_2{1,j}(i-5,1);  %��VL��n0 1 2 3 ״̬�ĵ���Ϊ0
          elseif(i==10)
             A(i,j)=P1_V0_2{1,j}(5,1)-E*(Ve-VL)*(2/delta_2^2)*P1_V0{1,j}(5,1);
             %��VL��n4״̬��Pi��������E*(Ve-VL)*(2/delta_2^2)*P4
              else
               A(i,j)=intP{1,j}; %��һ��
        end
    end
end
% Aa=size(A,1);
% while rank(A)~=Aa %���������H����������ʱ���ҳ�����������ļ����޹������newH����H
% [R,jb]=rref(A'); %�����RΪrref������Ŀ�������н���������jb��һ��������ΪĿ�������м����޹������ڵ�������ע��������Ҫ�Ƚ�Hת�ã���Ϊ���ǵ�Ŀ����Ҫ��H����������ļ����޹���
% newH=A(jb,:);
% A=newH;
% Aa=size(A,1);
% end
B=zeros(11,1);
B(11,1)=1;
K=A\B;
X=linspace(Ve_0,VL_0,1000);
X_P=cell(1,10);
P=cell(1,5);
for j=1:10
    X_P{1,j}=leijia2_V_noise(X,V_mid,data{1,j},ddai);
end
for i=1:5    
P{1,i}=K(1)*X_P{1,1}{1,i}+K(2)*X_P{1,2}{1,i}+K(3)*X_P{1,3}{1,i}+K(4)*X_P{1,4}{1,i}...
        +K(5)*X_P{1,5}{1,i}+K(6)*X_P{1,6}{1,i}+K(7)*X_P{1,7}{1,i}+K(8)*X_P{1,8}{1,i}...
        +K(9)*X_P{1,9}{1,i}+K(10)*X_P{1,10}{1,i};
end
Pp=P{1,5}+P{1,4}+P{1,3}+P{1,2}+P{1,1};
%�򿪸���
P0=trapz(X,P{1,1});
P1=trapz(X,P{1,2});
P2=trapz(X,P{1,3});
P3=trapz(X,P{1,4});
P4=trapz(X,P{1,5});
Pall=P0+P1+P2+P3+P4;
P0=P0/Pall;
P1=P1/Pall;
P2=P2/Pall;
P3=P3/Pall;
P4=P4/Pall;
figure(1)
subplot(2,3,1);
plot(X,P{1,1})
xlabel('V');ylabel('P(0,V)');
set(gca,'xlim',[Ve_0,VL_0]);%x�����᷶Χ
subplot(2,3,2);
plot(X,P{1,2})
xlabel('V');ylabel('P(1,V)');
set(gca,'xlim',[Ve_0,VL_0]);%x�����᷶Χ
% set(gca,'ylim',[0,0.025]);%y�����᷶Χ
subplot(2,3,3);
plot(X,P{1,3})
xlabel('V');ylabel('P(2,V)');
set(gca,'xlim',[Ve_0,VL_0]);%x�����᷶Χ
% set(gca,'ylim',[0,0.035]);%y�����᷶Χ
subplot(2,3,4);
plot(X,P{1,4})
xlabel('V');ylabel('P(3,V)');
set(gca,'xlim',[Ve_0,VL_0]);%x�����᷶Χ
% set(gca,'ylim',[0,0.016]);%y�����᷶Χ
subplot(2,3,5);
plot(X,P{1,5})
xlabel('V');ylabel('P(4,V)');
set(gca,'xlim',[Ve_0,VL_0]);%x�����᷶Χ
subplot(2,3,6);
plot(X,P{1,5}+P{1,4}+P{1,3}+P{1,2}+P{1,1})
%text(Vt,0.2,'\downarrowVt');text(Ve,0.2,'\downarrowVe');text(VL,0.2,'\downarrowVL')
xlabel('V');ylabel('P(V)');
set(gca,'xlim',[Ve_0,VL_0]);%x�����᷶Χ
% set(gca,'ylim',[0,0.08]);%y�����᷶Χ