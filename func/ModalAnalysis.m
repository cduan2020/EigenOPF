function [freq,damping_ratio,pfac,mode2gen]=ModalAnalysis(lambda,rv,lv)

T=find(imag(lambda)>=0);
lambda=lambda(T);
rv=rv(:,T);
lv=lv(:,T);

freq=abs(imag(lambda))/(2*pi);%in hz
damping_ratio=-real(lambda)./abs(lambda);
N_State=size(rv,1);
N_Machine=N_State/2;

for i=1:length(T)
    if(abs(lambda(i))<=1e-10)
    damping_ratio(i)=1;
    end
end
[Dr,Idx]=sort(damping_ratio*100);

pf1=zeros(length(T),N_State);
pf=zeros(length(T),N_State);
for i=1:length(T)
    pf1(i,:)=(rv(:,i).*(lv(:,i))).';
    pf(i,:)=abs(pf1(i,:)/max(abs(pf1(i,:))));
end
[PF,Indx]=sort(pf,2,'descend');
pfac=PF;
mode2gen=mod(Indx,N_State/2)+1;


ms3=rv(N_Machine+1:2*N_Machine,Idx(1:3)).';
eig2=[damping_ratio(Idx(1:3)) freq(Idx(1:3))].';

figure(2);

subplot(2,3,1);
compass(ms3(1,:).')
subplot(2,3,4);
feather(ms3(1,:).')
xlim([0 N_Machine+1]);
xlabel 'Machine No.'
ylabel 'Polar plot'
title(sprintf('Mode %d \n damping ratio = %3.2f percent \n frequency = %3.2f Hz',1,round(10000*eig2(1,1))/100,round(1000*eig2(2,1))/1000));


subplot(2,3,2);
compass(ms3(2,:).')
subplot(2,3,5);
feather(ms3(2,:).')
xlim([0 N_Machine+1]);
xlabel 'Machine No.'
ylabel 'Polar plot'
title(sprintf('Mode %d \n damping ratio = %3.2f percent \n frequency = %3.2f Hz',2,round(10000*eig2(1,2))/100,round(1000*eig2(2,2))/1000));

subplot(2,3,3);
compass(ms3(3,:).')
subplot(2,3,6);
feather(ms3(3,:).')
xlim([0 N_Machine+1]);
xlabel 'Machine No.'
ylabel 'Polar plot'
title(sprintf('Mode %d \n damping ratio = %3.2f percent \n frequency = %3.2f Hz',3,round(10000*eig2(1,3))/100,round(1000*eig2(2,3))/1000));

end