clear all, clc, close all

load('estimulos_BMA.mat');
X=estimulos_BMA;

Fs=512;
fftSize = 1024;

canal = 10;
[nv,nf] = size(X);

% EIXO DE FREQU�NCIA
vetor_eixo_freq = -Fs/2:Fs/fftSize:Fs/2-Fs/fftSize;

%GR�FICO 0 A 100HZ
freq_inicial = 0;
freq_final = 100;

indice_inicial = find(vetor_eixo_freq==freq_inicial);
indice_final = find(vetor_eixo_freq==freq_final);
fids_0_100 = indice_inicial:indice_final;

% LIMIARES DO RITMO ALFA
freq_inicial_alfa = 8;
freq_final_alfa = 13;

% MOSTRAR EIXO DE 8 A 13HZ 
indice_freq_inicial = find(vetor_eixo_freq==freq_inicial_alfa);
indice_freq_final = find(vetor_eixo_freq==freq_final_alfa);
fids_alfa = indice_freq_inicial:indice_freq_final;

% �NDICES DOS INTERVALOS DE FREQU�NCIA BAIXA, M�DIA E ALTA
z=[0 38 73 125];

hold all

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
      
    potencia_alfa_por_estimulo(indice_f) = mean(med(fids_alfa, indice_f));  
    
end

     
%%
%GR�FICO MOSTRANDO APENAS EXEMPLAR DE CADA GRUPO DE FREQU�NCIA
% Valores das frequ�ncias de est�mulo correspondentes aos �ndices 32, 56, 89 e 115
freq_legenda = [9, 21, 37.5, 66];
% Plot the power spectrum
figure;
%semilogy(freq_inicial:Fs/fftSize:freq_final, med(fids_0_100, 32));
plot(freq_inicial:Fs/fftSize:freq_final, med(fids_0_100, 32));
hold on
%semilogy(freq_inicial:Fs/fftSize:freq_final, med(fids_0_100, 56));
plot(freq_inicial:Fs/fftSize:freq_final, med(fids_0_100, 56));

hold on
%semilogy(freq_inicial:Fs/fftSize:freq_final, med(fids_0_100, 89));
plot(freq_inicial:Fs/fftSize:freq_final, med(fids_0_100, 89));
hold on
%semilogy(freq_inicial:Fs/fftSize:freq_final, med(fids_0_100, 115));
plot(freq_inicial:Fs/fftSize:freq_final, med(fids_0_100, 115));

% Adicione a legenda usando os valores das frequ�ncias de est�mulo
legend(num2str(freq_legenda(1)), num2str(freq_legenda(2)), num2str(freq_legenda(3)), num2str(freq_legenda(4)));

grid on;
title('Espectro de pot�ncia de quatro exemplares de est�mulos');
xlabel('Frequ�ncia (Hz)');
%ylabel('Magnitude (�V^2)');

hold off
%%

% Loop para gerar os gr�ficos separadamente para cada intervalo de frequ�ncia
for k = 1:length(z) - 1
    % Criar o gr�fico para cada intervalo de frequ�ncia
    figure;
    semilogy(freq_inicial:Fs/fftSize:freq_final, med(fids_0_100, z(k) + 1:z(k + 1)));
    
    xlabel('Frequ�ncia (Hz)');
    ylabel('Magnitude (�V^2)');
    title(['Espectros de Pot�ncia M�dios (Intervalo ', num2str(z(k) + 1), ' a ', num2str(z(k + 1)), ')']);
    grid on;
    xlim([0 100]);

    % Selecionar as frequ�ncias de est�mulo para o intervalo atual
    freq_estimulos_intervalo = freq_estimulos(z(k) + 1:z(k + 1));

    % Adicionar as legendas ao gr�fico
    legend(arrayfun(@(f) num2str(f), freq_estimulos_intervalo, 'UniformOutput', false));

    % Exibir o gr�fico
    hold off;
    
end



% M�DIA DOS ESPECTROS DE POT�NCIA DO INTERVALO
for k = 1:length(z)-1
    %Mmed(:,1) -> espectro de pot�ncia das baixas frequ�ncias de est�mulo
    %Mmed(:,2) -> espectro de pot�ncia das m�dias frequ�ncias de est�mulo
    %Mmed(:,3) -> espectro de pot�ncia das altas frequ�ncias de est�mulo
    Mmed(:,k) = sum(med(:,z(k)+1:z(k+1)),2)./(z(k+1)-z(k));  
end

%%
% Criar o gr�fico dos n�veis de frequ�ncia
figure;
semilogy(freq_inicial:Fs/fftSize:freq_final,Mmed(fids_0_100,:));

xlabel('Frequ�ncia (Hz)');
ylabel('Magnitude (�V^2)');
title('Espectros de Pot�ncia M�dios');
legend('Baixa Frequ�ncia', 'M�dia Frequ�ncia', 'Alta Frequ�ncia');
grid on;
xlim([0 100]);

% Exibir o gr�fico
hold off;





%%
% Criar o gr�fico dos n�veis de frequ�ncia
figure;
semilogy(freq_inicial:Fs/fftSize:freq_final,Mmed(fids_0_100,:));

xlabel('Frequ�ncia (Hz)');
ylabel('Magnitude (�V^2)');
title('Espectros de Pot�ncia M�dios');
legend('Baixa Frequ�ncia', 'M�dia Frequ�ncia', 'Alta Frequ�ncia');
grid on;
xlim([0 100]);

% Exibir o gr�fico
hold off;

%% GR�FICO DE BARRAS PARA CADA RITMO

% �ndices das frequ�ncias das bandas de interesse
freq_delta = [1 4];  % Delta (1-4 Hz)
freq_teta = [4 8];   % Teta (4-8 Hz)
freq_alfa = [8 13];  % Alfa (8-13 Hz)
freq_beta = [13 30]; % Beta (13-30 Hz)

% Calcular os �ndices correspondentes a cada frequ�ncia de interesse
indice_freq_delta = find(vetor_eixo_freq >= freq_delta(1) & vetor_eixo_freq <= freq_delta(2));
indice_freq_teta = find(vetor_eixo_freq >= freq_teta(1) & vetor_eixo_freq <= freq_teta(2));
indice_freq_alfa = find(vetor_eixo_freq >= freq_alfa(1) & vetor_eixo_freq <= freq_alfa(2));
indice_freq_beta = find(vetor_eixo_freq >= freq_beta(1) & vetor_eixo_freq <= freq_beta(2));

ritmos = {'Delta', 'Teta', 'Alfa', 'Beta'};
freq_indices = {indice_freq_delta, indice_freq_teta, indice_freq_alfa, indice_freq_beta};


% Plotar gr�ficos de barras por ritmo
for ritmo_indice = 1:length(ritmos)
    ritmo = ritmos{ritmo_indice};
    indice_freq = freq_indices{ritmo_indice};
    
    % Calcular a m�dia do PSD para cada banda em cada n�vel de frequ�ncia
    med_PSD_baixa = mean(Mmed(indice_freq(1):indice_freq(2), 1));
    med_PSD_media = mean(Mmed(indice_freq(1):indice_freq(2), 2));
    med_PSD_alta = mean(Mmed(indice_freq(1):indice_freq(2), 3));
    
    % Plotar o gr�fico de barras
    figure;
    bar([med_PSD_baixa, med_PSD_media, med_PSD_alta])
    xticklabels({'Baixa', 'Media', 'Alta'})
    %bar([med_PSD_media, med_PSD_alta])
    %xticklabels({'Media', 'Alta'})
    title(['M�dia do PSD na banda ', ritmo, ' para os tr�s n�veis de frequ�ncia'])
    ylabel('PSD')
    
end






