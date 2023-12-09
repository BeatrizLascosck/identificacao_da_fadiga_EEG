clear all, clc, close all
% 
load('pre_proc_estimulos_BMA.mat');
X=pre_proc_estimulos_BMA;
[nv,nf] = size(X);

grupos_freq=[0 38 73 125];
fases_B = [0 5 9 13];
fases_M = [0 4 9 13];
fases_A = [0 4 8 13];
canal=10;
Fs=512;
fftSize = 1024;

Simulacoes = cell(1, 4);

% EIXO DE FREQUÊNCIA
vetor_eixo_freq = -Fs/2:Fs/fftSize:Fs/2-Fs/fftSize;

%% SIMULAÇÕES

for fase = 1:3
    S = [];
    %Baixa
    for indice_f = grupos_freq(fase)+1:grupos_freq(fase+1)
        for v_B = fases_B(fase)+1:fases_B(fase+1)
            sinal = X{v_B,indice_f};
            S(:,v_B)=fft(sinal(:,canal),fftSize);         
        end     
    end
    
     %Média
    for indice_f = grupos_freq(fase)+1:grupos_freq(fase+1)
        for v_M = fases_M(fase)+1:fases_M(fase+1)
            sinal = X{v_M,indice_f};
            S(:,(size(S, 2)+1))=fft(sinal(:,canal),fftSize);         
        end      
    end
    
    %Alta
    for indice_f = grupos_freq(fase)+1:grupos_freq(fase+1)
        for v_A = fases_A(fase)+1:fases_A(fase+1)
            sinal = X{v_A,indice_f};
            S(:,(size(S, 2)+1))=fft(sinal(:,canal),fftSize);         
        end            
    end
    
   for indice_f = 1:nf
    % ESPECTRO DE POTÊNCIA MÉDIO PARA CADA SIMULAÇÃO
    % Média de todos os voluntários para cada frequência de estímulo
    simulacao_med(:,indice_f) = abs(sum(S.^2,2)./nv);
   end
 
   Simulacoes{1,fase} = simulacao_med;

end


for indice_f = 1:nf
    for v = 1:nv
        sinal = X{v,indice_f};
        S(:,v)=fft(sinal(:,canal),fftSize);         
    end
    todos_med(:,indice_f) = abs(sum(S.^2,2)./nv);
end

Simulacoes{1,4} = todos_med;

%% GRÁFICO 0 A 100HZ
freq_inicial = 0;
freq_final = 100;

indice_inicial = find(vetor_eixo_freq==freq_inicial);
indice_final = find(vetor_eixo_freq==freq_final);
fids_0_100 = indice_inicial:indice_final;


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


for fase = 1:4
    med = [];
    Mmed=[];
    
    med = Simulacoes{1,fase};
    
    % MÉDIA DOS ESPECTROS DE POTÊNCIA DO INTERVALO
    for k = 1:length(grupos_freq)-1
        Mmed(:,k) = sum(med(:,grupos_freq(k)+1:grupos_freq(k+1)),2)./(grupos_freq(k+1)-grupos_freq(k));  
    end

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
        title(['Média do ritmo ', ritmo, ' para os três níveis de frequência na Simulação ',num2str(fase)])

    end
   
end


