clear all, clc, close all

load('estimulos_BMA.mat');
X=estimulos_BMA;

Fs=512;
fftSize = 1024;

canal = 10;
[nv,nf] = size(X);

% EIXO DE FREQUÊNCIA
vetor_eixo_freq = -Fs/2:Fs/fftSize:Fs/2-Fs/fftSize;

%GRÁFICO 0 A 100HZ
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

% ÍNDICES DOS INTERVALOS DE FREQUÊNCIA BAIXA, MÉDIA E ALTA
z=[0 38 73 125];

hold all

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
      
    potencia_alfa_por_estimulo(indice_f) = mean(med(fids_alfa, indice_f));  
    
end

     
%%
%GRÁFICO MOSTRANDO APENAS EXEMPLAR DE CADA GRUPO DE FREQUÊNCIA
% Valores das frequências de estímulo correspondentes aos índices 32, 56, 89 e 115
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

% Adicione a legenda usando os valores das frequências de estímulo
legend(num2str(freq_legenda(1)), num2str(freq_legenda(2)), num2str(freq_legenda(3)), num2str(freq_legenda(4)));

grid on;
title('Espectro de potência de quatro exemplares de estímulos');
xlabel('Frequência (Hz)');
%ylabel('Magnitude (µV^2)');

hold off
%%

% Loop para gerar os gráficos separadamente para cada intervalo de frequência
for k = 1:length(z) - 1
    % Criar o gráfico para cada intervalo de frequência
    figure;
    semilogy(freq_inicial:Fs/fftSize:freq_final, med(fids_0_100, z(k) + 1:z(k + 1)));
    
    xlabel('Frequência (Hz)');
    ylabel('Magnitude (µV^2)');
    title(['Espectros de Potência Médios (Intervalo ', num2str(z(k) + 1), ' a ', num2str(z(k + 1)), ')']);
    grid on;
    xlim([0 100]);

    % Selecionar as frequências de estímulo para o intervalo atual
    freq_estimulos_intervalo = freq_estimulos(z(k) + 1:z(k + 1));

    % Adicionar as legendas ao gráfico
    legend(arrayfun(@(f) num2str(f), freq_estimulos_intervalo, 'UniformOutput', false));

    % Exibir o gráfico
    hold off;
    
end



% MÉDIA DOS ESPECTROS DE POTÊNCIA DO INTERVALO
for k = 1:length(z)-1
    %Mmed(:,1) -> espectro de potência das baixas frequências de estímulo
    %Mmed(:,2) -> espectro de potência das médias frequências de estímulo
    %Mmed(:,3) -> espectro de potência das altas frequências de estímulo
    Mmed(:,k) = sum(med(:,z(k)+1:z(k+1)),2)./(z(k+1)-z(k));  
end

%%
% Criar o gráfico dos níveis de frequência
figure;
semilogy(freq_inicial:Fs/fftSize:freq_final,Mmed(fids_0_100,:));

xlabel('Frequência (Hz)');
ylabel('Magnitude (µV^2)');
title('Espectros de Potência Médios');
legend('Baixa Frequência', 'Média Frequência', 'Alta Frequência');
grid on;
xlim([0 100]);

% Exibir o gráfico
hold off;





%%
% Criar o gráfico dos níveis de frequência
figure;
semilogy(freq_inicial:Fs/fftSize:freq_final,Mmed(fids_0_100,:));

xlabel('Frequência (Hz)');
ylabel('Magnitude (µV^2)');
title('Espectros de Potência Médios');
legend('Baixa Frequência', 'Média Frequência', 'Alta Frequência');
grid on;
xlim([0 100]);

% Exibir o gráfico
hold off;

%% GRÁFICO DE BARRAS PARA CADA RITMO

% Índices das frequências das bandas de interesse
freq_delta = [1 4];  % Delta (1-4 Hz)
freq_teta = [4 8];   % Teta (4-8 Hz)
freq_alfa = [8 13];  % Alfa (8-13 Hz)
freq_beta = [13 30]; % Beta (13-30 Hz)

% Calcular os índices correspondentes a cada frequência de interesse
indice_freq_delta = find(vetor_eixo_freq >= freq_delta(1) & vetor_eixo_freq <= freq_delta(2));
indice_freq_teta = find(vetor_eixo_freq >= freq_teta(1) & vetor_eixo_freq <= freq_teta(2));
indice_freq_alfa = find(vetor_eixo_freq >= freq_alfa(1) & vetor_eixo_freq <= freq_alfa(2));
indice_freq_beta = find(vetor_eixo_freq >= freq_beta(1) & vetor_eixo_freq <= freq_beta(2));

ritmos = {'Delta', 'Teta', 'Alfa', 'Beta'};
freq_indices = {indice_freq_delta, indice_freq_teta, indice_freq_alfa, indice_freq_beta};


% Plotar gráficos de barras por ritmo
for ritmo_indice = 1:length(ritmos)
    ritmo = ritmos{ritmo_indice};
    indice_freq = freq_indices{ritmo_indice};
    
    % Calcular a média do PSD para cada banda em cada nível de frequência
    med_PSD_baixa = mean(Mmed(indice_freq(1):indice_freq(2), 1));
    med_PSD_media = mean(Mmed(indice_freq(1):indice_freq(2), 2));
    med_PSD_alta = mean(Mmed(indice_freq(1):indice_freq(2), 3));
    
    % Plotar o gráfico de barras
    figure;
    bar([med_PSD_baixa, med_PSD_media, med_PSD_alta])
    xticklabels({'Baixa', 'Media', 'Alta'})
    %bar([med_PSD_media, med_PSD_alta])
    %xticklabels({'Media', 'Alta'})
    title(['Média do PSD na banda ', ritmo, ' para os três níveis de frequência'])
    ylabel('PSD')
    
end






