clear all, clc, close all

load('estimulosBMA.mat');

%load('Simulacao_3.mat');
%X=Simulacao_3;

Fs=512;
fftSize = 1024;

canal = 10;
[nv,nf] = size(X);

% EIXO DE FREQUÊNCIA
vetor_eixo_freq = -Fs/2:Fs/fftSize:Fs/2-Fs/fftSize;

% Indices of the frequencies of interest
freq_values_of_interest = [9, 21, 37.5, 66];


for indice_f = 1:nf
    
    % VETOR COM VALORES DAS FREQUÊNCIAS DE ESTÍMULOS
    freq_estimulos(indice_f) = X{1,indice_f}.freq;
    
    for v = 1:nv
        sinal = X{v,indice_f}.signal;
        % MATRIZ DE ESPECTRO DE FREQUÊNCIA
        % Transformada Rápida de Fourier (FFT)
        % Apenas do canal 10
        S(:,v)=fft(sinal(:,canal),fftSize);         
    end
    % ESPECTRO DE POTÊNCIA MÉDIO 
    % Média de todos os voluntários para cada frequência de estímulo
    med(:,indice_f) = abs(sum(S.^2,2)./nv);
      
    
end



% Loop through the frequency values of interest
for freq_val = freq_values_of_interest
    % Find the index in the frequency values that is closest to the desired frequency
    [~, freq_index] = min(abs(freq_estimulos - freq_val));
    
    % Create a new figure for each frequency value
    figure;
    for v = 1:nv
        sinal = X{v, freq_index}.signal;
        
        % Plot the EEG signal for the current frequency and volunteer
        subplot(nv, 1, v);
        plot(sinal(:, canal));
        title(['Sinal EEG - Frequência: ' num2str(freq_val) ' Hz - Voluntário: ' num2str(v)]);
        xlabel('Amostras');
        ylabel('Amplitude');
    end
end
