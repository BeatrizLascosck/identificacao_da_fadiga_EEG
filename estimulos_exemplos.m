clear all, clc, close all

load('estimulosBMA.mat');

%load('Simulacao_3.mat');
%X=Simulacao_3;

Fs=512;
fftSize = 1024;

canal = 10;
[nv,nf] = size(X);

% EIXO DE FREQU�NCIA
vetor_eixo_freq = -Fs/2:Fs/fftSize:Fs/2-Fs/fftSize;

% Indices of the frequencies of interest
freq_values_of_interest = [9, 21, 37.5, 66];


for indice_f = 1:nf
    
    % VETOR COM VALORES DAS FREQU�NCIAS DE EST�MULOS
    freq_estimulos(indice_f) = X{1,indice_f}.freq;
    
    for v = 1:nv
        sinal = X{v,indice_f}.signal;
        % MATRIZ DE ESPECTRO DE FREQU�NCIA
        % Transformada R�pida de Fourier (FFT)
        % Apenas do canal 10
        S(:,v)=fft(sinal(:,canal),fftSize);         
    end
    % ESPECTRO DE POT�NCIA M�DIO 
    % M�dia de todos os volunt�rios para cada frequ�ncia de est�mulo
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
        title(['Sinal EEG - Frequ�ncia: ' num2str(freq_val) ' Hz - Volunt�rio: ' num2str(v)]);
        xlabel('Amostras');
        ylabel('Amplitude');
    end
end
