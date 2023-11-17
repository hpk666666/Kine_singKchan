tic;
clear;
Iex=0.6*10^-12;%�ⲿ��������λ A
M_area=6*10^-12;%����Ĥ�������λ��um^2
Cm=0.01*M_area;%��λ�ǵ��� pF
gL=3*M_area;%��λ�ǵ絼�� pS
gK=20*10^-12;%��λ�ǵ絼�� S
VL=-54.4*10^-3;%��λ�� mV
V_K=-77*10^-3;%��λ�� mV
ge=gL+gK;%��λ�ǵ絼�� pS
V_L=(gL*VL+Iex)/gL;%��λ�� V
V_e=(gL*VL+gK*V_K+Iex)/ge;%��λ�� mV
Cond_L=gL/Cm;
Cond_k=ge/Cm;
%Alpha_V=0.01*(V+55)/(1-exp(-(V+55)/10));
%Beta_V=0.125*exp(-(V+65)/80);
%K01=4*Alpha_V;
%K10=Beta_V;K12=3*Alpha_V;
%K21=2*Beta_V;K23=2*Alpha_V;
%K32=3*Beta_V;K34=Alpha_V;
%K43=4*Beta_V;
%ת������
%lambda_t=[K01,K10+K12,K21+K23,K32+K34,K43];
%lambda_tt=0;
%lambda_tt=[K10/(K10+K12),K21/(K21+K23),K32/(K32+K34)];

delta_t=0.001*10^-3;
T_step=2^18;%2^20
% T=linspace(0,0.001*2^19*10^-3,2^19);
V=zeros(1,T_step);
S_k=zeros(1,T_step);
%����Ĳ���ֵҪע��һ��.  ���ﶨ����:C=1
V(1)=-55*10^-3;
n_i=1;
while n_i<T_step
    
   switch S_k(n_i)
      %-------------------------  
              case 0
                  %calculate time that the channel stay at "0"   lambda=4*alfa_n; 
                   U1=rand;
                   F=1;
                   integ_lambda=0;        %the integrel of lambda(s)
                   while F>U1  && n_i<T_step
                       V(n_i+1)=V(n_i)+(-Cond_L*(V(n_i)-V_L))*delta_t;
                       S_k(n_i+1)=S_k(n_i);
                       alfa_n=0.01*10^6*(V(n_i)+55*10^-3)/(1-exp(-(V(n_i)+55*10^-3)/0.01));
                       lambda=4*alfa_n;
                       integ_lambda=integ_lambda+lambda*delta_t;
                       F=exp(-integ_lambda);
                       n_i=n_i+1;
                   end
            
           
                 %at the jump time , system jump to "1"; definitely
            
                 S_k(n_i)=1;
         
            
              %-------------------------  
              case 1
                  %calculate time that the channel stay at "1"   lambda=3*alfa_n+1*bata_n; 
                   U1=rand;
                   F=1;
                   integ_lambda=0;        %the integrel of lambda(s)
                   while F>U1   &&   n_i<T_step
                       V(n_i+1)=V(n_i)+(-Cond_L*(V(n_i)-V_L))*delta_t;
                       S_k(n_i+1)=S_k(n_i);
                       alfa_n=0.01*10^6*(V(n_i)+55*10^-3)/(1-exp(-(V(n_i)+55*10^-3)/0.01));
                       bata_n=0.125*10^3*exp(-(V(n_i)+65*10^-3)/0.08); 
                       lambda=3*alfa_n+bata_n;
                       integ_lambda=integ_lambda+lambda*delta_t;
                       F=exp(-integ_lambda);
                       n_i=n_i+1;
                   end
           
                 
                 %at the jump time , system jump to "0" or "2"; with
                 %probability ...rate            
            
                U2=rand;
                if U2<=3*alfa_n/lambda; 
                    S_k(n_i)=2;   %jump to state "2"
                else
                   S_k(n_i)=0;   %jump to state "0"
                 end 
            
            
          %-------------------------  
              case 2
                   %calculate time that the channel stay at "2"   lambda=2*alfa_n+2*bata_n; 
                   U1=rand;
                   F=1;
                   integ_lambda=0;        %the integrel of lambda(s)
                   while F>U1  && n_i<T_step
                       V(n_i+1)=V(n_i)+(-Cond_L*(V(n_i)-V_L))*delta_t;
                       S_k(n_i+1)=S_k(n_i);
                       alfa_n=0.01*10^6*(V(n_i)+55*10^-3)/(1-exp(-(V(n_i)+55*10^-3)/0.01));
                       bata_n=0.125*10^3*exp(-(V(n_i)+65*10^-3)/0.08); 
                       lambda=2*alfa_n+2*bata_n;
                       integ_lambda=integ_lambda+lambda*delta_t;
                       F=exp(-integ_lambda);
                       n_i=n_i+1;
                   end
           
                  %at the jump time , system jump to "1" or "3"; with
                 %probability ...rate 
            
            
            
                   U2=rand;
                   if U2<=2*alfa_n/lambda;
                       S_k(n_i)=3;   %jump to state "3"
                   else
                       S_k(n_i)=1;   %jump to state "1"
                   end 
            
            
            
             %-------------------------  3
              case 3
                   %calculate time that the channel stay at "3"   lambda=1*alfa_n+3*bata_n; 
                   U1=rand;
                   F=1;
                   integ_lambda=0;        %the integrel of lambda(s)
                   while F>U1  &&   n_i<T_step
                       V(n_i+1)=V(n_i)+(-Cond_L*(V(n_i)-V_L))*delta_t;
                       S_k(n_i+1)=S_k(n_i);
                       alfa_n=0.01*10^6*(V(n_i)+55*10^-3)/(1-exp(-(V(n_i)+55*10^-3)/0.01));
                       bata_n=0.125*10^3*exp(-(V(n_i)+65*10^-3)/0.08); 
                       lambda=alfa_n+3*bata_n;
                       integ_lambda=integ_lambda+lambda*delta_t;
                       F=exp(-integ_lambda);
                       n_i=n_i+1;
                   end
           
                  %at the jump time , system jump to "2" or "4"; with
                 %probability ...rate 
            
                   U2=rand;
                   if U2<=alfa_n/lambda;
                       S_k(n_i)=4;   %jump to state "4"
                   else
                       S_k(n_i)=2;   %jump to state "2"
                   end 
         
            
            
             %-------------------------  4
              case 4
                    %calculate time that the channel stay at "4"   lambda=4*bata_n; 
                   U1=rand;
                   F=1;
                   integ_lambda=0;        %the integrel of lambda(s)
                   while F>U1    &&   n_i<T_step
                       V(n_i+1)=V(n_i)+(-Cond_k*(V(n_i)-V_e))*delta_t;
                       S_k(n_i+1)=S_k(n_i); 
                       bata_n=0.125*10^3*exp(-(V(n_i)+65*10^-3)/0.08); 
                       lambda=4*bata_n;
                       integ_lambda=integ_lambda+lambda*delta_t;
                       F=exp(-integ_lambda);
                       n_i=n_i+1;
                   end
             %   %time that the system stay in state "4"
            %at jump time , system jump to  "3"; definitely
                  S_k(n_i)=3;
         
            
   end

end
subplot(2,1,1);
plot(S_k);   
xlabel('Time(ms)');
ylabel('States');
set(gca,'ytick',[0 1 2 3 4]) %x�������Ͽ̶ȵ����ݵ�λ��
set(gca,'yticklabel',{'N0','N1','N2','N3','N4'});
% set(gca,'xlim',[0,500000]);%x�����᷶Χ
subplot(2,1,2);                                                                                 
plot(V);   
xlabel('Time(ms)');
ylabel('Voltage(mV)');
% set(gca,'xlim',[0,500000]);%x�����᷶Χ
%������������̶�*0.001,�����ԭ��������T_step,������T_step*delta_t