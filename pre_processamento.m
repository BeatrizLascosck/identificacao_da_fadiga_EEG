clear all, clc, close all
addpath('fases');



load('estimulos_BMA.mat');
[nv,nf] = size(estimulos_BMA);

Fs=512;
fftSize = 1024;

pre_proc_estimulos_BMA = cell(nv, nf);

for indice_f = 1:nf
    for v = 1:nv
        eeg_signal = estimulos_BMA{v, indice_f}.signal;
        [A,B]= size(eeg_signal);
        
        %Retirando componente DC dos sinais  
        valorDC = mean(eeg_signal);
        vetDC = repmat(valorDC,A,1);
        eeg_signal_semDC = eeg_signal - vetDC;

        % Fazendo uma filtragem espacial CAR
        media_car = mean(eeg_signal_semDC,2);
        vetcar = repmat(media_car,1,B);
        eeg_signal_car = eeg_signal_semDC - vetcar;

        [nSamples,canais] = size(eeg_signal_car);

        % Aplicando o filtro Notch
        wo = 50/(Fs/2);  % frequência central do filtro em 50 Hz
        bw = wo/10;  % largura de banda em 10 Hz
        [b,a] = iirnotch(wo,bw);
         for k = 1:canais
            notch_signal(:,k) = filter(b, a, eeg_signal_car(:,k)); 
         end
         pre_proc_estimulos_BMA{v, indice_f} = notch_signal;
    end
end





