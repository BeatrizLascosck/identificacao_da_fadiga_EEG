%clear all
clc
close all

Fs=512;

%load('sinais_freq.mat')
%load('frequencias_estimulos.mat')



for v = 1:1 % 1 até 20 são os voluntários, 21 é a média
    % v => para cada voluntário
    %v = 21; %média dos voluntários e do canal
    
    Alpha_estimulos=zeros(1, 125);
    
    for f = 1:125
        % f => para cada frequência de estímulo
        
        sinal_EEG=sinais_freq{v,f};
   
        [sigma_power,teta_power,alpha_power,beta_power] = PSD_banda(sinal_EEG,Fs);
        Alpha_estimulos(1,f)=alpha_power;   
    end
   
   
    sinal=sinais_freq{1,1};
    %sinal(:,12)
    Vol_power_density(sinais_freq(v,:),Fs,v,frequencias_estimulos)
  
    %Vol_estimulos_ritmo('Alpha',Alpha_estimulos,frequencias_estimulos,v)
        
 end
