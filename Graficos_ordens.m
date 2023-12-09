close all, clc 
%clear all

%%
% simulacoes =cell(1, 5);
% 
% load('G:\Funcoes_TCC\fases\simulacaoG2_1.mat')
% simulacoes{1,1}=simulacaoG2_1;
% 
% load('G:\Funcoes_TCC\fases\simulacaoG2_2.mat')
% simulacoes{1,2}=simulacaoG2_2;
% 
% load('G:\Funcoes_TCC\fases\simulacaoG2_3.mat')
% simulacoes{1,3}=simulacaoG2_3;
% 
% load('G:\Funcoes_TCC\fases\simulacaoG2_4.mat')
% simulacoes{1,4}=simulacaoG2_4;
% 
% 
% load('G:\Funcoes_TCC\fases\dados_baseG2.mat')
% simulacoes{1,5}=dados_baseG2;
%%

ordens =cell(1, 7);

load('G:\Funcoes_TCC\fases\ABM\ABM.mat')
ordens{1,1}=ordem;

load('G:\Funcoes_TCC\fases\AMB\AMB.mat')
ordens{1,2}=ordem;

load('G:\Funcoes_TCC\fases\BAM\BAM.mat')
ordens{1,3}=ordem;

load('G:\Funcoes_TCC\fases\BMA\BMA.mat')
ordens{1,4}=ordem;

load('G:\Funcoes_TCC\fases\MAB\MAB.mat')
ordens{1,5}=ordem;

load('G:\Funcoes_TCC\fases\MBA\MBA.mat')
ordens{1,6}=ordem;

load('G:\Funcoes_TCC\dadosArg\Sinal Base\dados_base.mat')
ordens{1,7}=dados_base;

%%

Fs=512;
fftSize = 1024;
canal = 12;

% EIXO DE FREQUÊNCIA
vetor_eixo_freq = fftshift(-Fs/2:Fs/fftSize:Fs/2-Fs/fftSize);

%GRÁFICO 0 A 100HZ
freq_inicial = 0;
freq_final = 100;

indice_inicial = find(vetor_eixo_freq==freq_inicial);
indice_final = find(vetor_eixo_freq==freq_final);
fids_0_100 = indice_inicial:indice_final;


% ÍNDICES DOS INTERVALOS DE FREQUÊNCIA BAIXA, MÉDIA E ALTA
z=[0 38 73 125];

load('frequencias_estimulos.mat');
freq_estimulos = frequencias_estimulos;

indices_niveis_rodadas = cell(1, 4);
delta_niveis_rodadas = cell(1, 4);
teta_niveis_rodadas = cell(1, 4);
alfa_niveis_rodadas = cell(1, 4);
beta_niveis_rodadas = cell(1, 4);

%%
%SINAL BASE
%X=simulacoes(5);
X=ordens(7);
[nv,c] = size(X{1,1});

%%
% Inicialize a variável pre_proc_X_base
pre_proc_X_base = cell(nv, c);

%%
for v = 1:nv
    
    eeg_signal = X{1,1}{v,1}.signal;
    [A,B]= size(eeg_signal);

    %Retirando componente DC dos sinais  
    valorDC = mean(eeg_signal);
    vetDC = repmat(valorDC,A,1);
    eeg_signal_semDC = eeg_signal - vetDC;

    %Fazendo uma filtragem espacial CAR
    media_car = mean(eeg_signal_semDC,2);
    vetcar = repmat(media_car,1,B);
    eeg_signal_car = eeg_signal_semDC - vetcar;

     pre_proc_X_base{v, 1} = eeg_signal_car;
end

S=[];
for v = 1:nv
    
    sinal = pre_proc_X_base{v,1};
    % MATRIZ DE ESPECTRO DE FREQUÊNCIA
    % Transformada Rápida de Fourier (FFT)
    % Apenas do canal 12
    S(:,v)=fft(sinal(:,canal),fftSize);
   
end
% ESPECTRO DE POTÊNCIA MÉDIO 
% Média de todos os voluntários para cada frequência de estímulo
base_med = abs(sum(S.^2,2)./nv);

figure

semilogy(freq_inicial:Fs/fftSize:freq_final,base_med(fids_0_100,1));
title(['Média do sinal base']);
grid on;
xlim([0 100]);


%%
for fase = 1:6
    X=ordens(fase);
    [nv,nf] = size(X{1,1});

  % PRÉ-PROCESSAMENTO
    pre_proc_X = cell(nv, nf);

    for indice_f = 1:nf
        for v = 1:nv
            
            eeg_signal = X{1,1}{v, indice_f}.signal;
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
%             wo = 50/(Fs/2);  % frequência central do filtro em 50 Hz
%             bw = wo/10;  % largura de banda em 10 Hz
%             [b,a] = iirnotch(wo,bw);
%              for k = 1:canais
%                 notch_signal(:,k) = filter(b, a, eeg_signal_car(:,k)); 
%              end
             pre_proc_X{v, indice_f} = eeg_signal_car;
            
        end
    end
  %%  
    for indice_f = 1:nf
        S=[];
        for v = 1:nv
            sinal = pre_proc_X{v,indice_f};
            % MATRIZ DE ESPECTRO DE FREQUÊNCIA
            % Transformada Rápida de Fourier (FFT)
            % Apenas do canal 12
            S(:,v)=fft(sinal(:,canal),fftSize);          
        end
        % ESPECTRO DE POTÊNCIA MÉDIO 
        % Média de todos os voluntários para cada frequência de estímulo
       
        med(:,indice_f) = abs(sum(S.^2,2)./nv);
    end
    
    %% MÉDIA DOS ESPECTROS DE POTÊNCIA DO INTERVALO

    for k = 1:length(z)-1
        Mmed(:,k) = sum(med(:,z(k)+1:z(k+1)),2)./(z(k+1)-z(k));  
    end
    
    %% Criar o gráfico dos níveis de frequência

    nome_figura = sprintf('Figura_1_%d', fase);
    figure('Name', nome_figura);

    semilogy(freq_inicial:Fs/fftSize:freq_final,Mmed(fids_0_100,:));
    hold on;
    semilogy(freq_inicial:Fs/fftSize:freq_final,base_med(fids_0_100,1));
 
   
    xlabel('Frequência (Hz)');
    ylabel('Magnitude (µV^2)');
    %title(['Espectros de potência médios para os três níveis de frequência na Rodada ',num2str(fase)']);
    %title(['Espectros de Potência na Rodada ',num2str(fase)']);
    title(['Espectros de Potência na Ordem ',num2str(fase)']);
    legend('Baixa Frequência', 'Média Frequência', 'Alta Frequência','Sem Estímulo');
    grid on;
    xlim([0 100]);


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
    ritmos_fase = [];
    ritmos_base = [];
    for ritmo_indice = 1:length(ritmos)
        ritmo = ritmos{ritmo_indice};
        indice_freq = freq_indices{ritmo_indice};

        % Calcular a média do PSD para cada banda em cada nível de frequência
        med_PSD_baixa = mean(Mmed(indice_freq(1):indice_freq(2), 1));
        med_PSD_media = mean(Mmed(indice_freq(1):indice_freq(2), 2));
        med_PSD_alta = mean(Mmed(indice_freq(1):indice_freq(2), 3));
        
        ritmo_base = mean(base_med(indice_freq(1):indice_freq(2)));
        
        y = [med_PSD_baixa med_PSD_media med_PSD_alta];
        
        % Plotar o gráfico de barras
        nome_figura = sprintf('Figura_2_%d', fase);
        figure('Name', nome_figura);


        b = bar(y,'FaceColor','flat');
        b.CData(1,:) = [0 0.4470 0.7410];
        b.CData(2,:) = [0.8500 0.3250 0.0980];
        b.CData(3,:) = [0.9290 0.6940 0.1250];
        
        xticklabels({'Baixa', 'Media', 'Alta'})
        %title(['Média do ritmo ', ritmo, ' para os três níveis de frequência na Rodada ',num2str(fase)'])
        %title(['Ritmo ', ritmo, ' na Rodada ',num2str(fase)'])
        title(['Ritmo ', ritmo, ' na Ordem ',num2str(fase)'])
        
        hold on;
        base = plot([0, 4], [ritmo_base, ritmo_base], 'r-', 'Color', [0.5, 0, 0.5], 'LineWidth', 2);
        legend([base], {'Sem Estímulo'})
        
        ritmo_BMA = [med_PSD_baixa,med_PSD_media,med_PSD_alta];
        ritmos_fase = vertcat(ritmos_fase, ritmo_BMA);
        ritmos_base = vertcat(ritmos_base, ritmo_base);
        
        if ritmo_indice ==1
            Delta_niveis_rodadas{1, fase} = ritmo_BMA;
        end

        if ritmo_indice ==2
            Teta_niveis_rodadas{1, fase} = ritmo_BMA;
        end        
        
        if ritmo_indice ==3
            Alfa_niveis_rodadas{1, fase} = ritmo_BMA;
        end
        
        if ritmo_indice ==4
            Beta_niveis_rodadas{1, fase} = ritmo_BMA;
        end
        
        ritmo_BMA = [];
              
    end
    
    %disp(['RODADA ' num2str(fase)])
    %disp(ritmos_fase)
    
    % indice_1 => teta/alfa
    % indice_2 => beta/alfa
    % indice_3 => (teta + alfa)/beta
   %% 

    indices_niveis_f = [];
    for nivel_f = 1:3
        indice_1 = ritmos_fase(2,nivel_f)/ritmos_fase(3,nivel_f);
        indice_2 = ritmos_fase(4,nivel_f)/ritmos_fase(3,nivel_f);
        indice_3 = (ritmos_fase(2,nivel_f) + ritmos_fase(3,nivel_f))/ritmos_fase(4,nivel_f);
        indices_nivel_f = [indice_1;indice_2;indice_3];
        indices_niveis_f = [indices_niveis_f,indices_nivel_f];
    end
    
    indices_niveis_rodadas{1, fase} = indices_niveis_f;
    
    indice_1_base = ritmos_base(2)/ritmos_base(3);
    indice_2_base = ritmos_base(4)/ritmos_base(3);
    indice_3_base = (ritmos_base(2) + ritmos_base(3))/ritmos_base(4);
    
    indices_base = [indice_1_base, indice_2_base, indice_3_base];
    
    
    for indice = 1:3
        % Plotar o gráfico de barras

        nome_figura = sprintf('Figura_3_%d', fase);
        figure('Name', nome_figura);
        
        y = indices_niveis_f(indice,:);
        
        b = bar(y,'FaceColor','flat');
        b.CData(1,:) = [0 0.4470 0.7410];
        b.CData(2,:) = [0.8500 0.3250 0.0980];
        b.CData(3,:) = [0.9290 0.6940 0.1250];
        
        hold on;
        
        ind_numeric = indices_base(indice);
        
        base = plot([0, 4], [ind_numeric, ind_numeric], 'r-', 'Color', [0.5, 0, 0.5], 'LineWidth', 2);
        legend([base], {'Sem Estímulo'})
        
        xticklabels({'Baixa', 'Media', 'Alta'})
        %title(['Média do ritmo ', ritmo, ' para os três níveis de frequência na Rodada ',num2str(fase)'])
        %title(['Indice ', num2str(indice), ' na Rodada ',num2str(fase)'])
        title(['Indice ', num2str(indice), ' na Ordem ',num2str(fase)'])
        
    end
  
end

%%
for indice = 1:3
    indices_B = [indices_niveis_rodadas{1,1}(indice,1),indices_niveis_rodadas{1,2}(indice,1),indices_niveis_rodadas{1,3}(indice,1)];
    indices_M = [indices_niveis_rodadas{1,1}(indice,2),indices_niveis_rodadas{1,2}(indice,2),indices_niveis_rodadas{1,3}(indice,2)];
    indices_A = [indices_niveis_rodadas{1,1}(indice,3),indices_niveis_rodadas{1,2}(indice,3),indices_niveis_rodadas{1,3}(indice,3)];

    y_Baixa = indices_B; 
    y_Media = indices_M;
    y_Alta = indices_A;


    nome_figura = sprintf('Figura_4_%d', indice);
    figure('Name', nome_figura);
    
    %Gráfico de dispersão
    plot(y_Baixa, '-o');
    hold on; % Mantém o gráfico atual ativo

    plot(y_Media, '-s');
    plot(y_Alta, '-d');  
    
    
    hold on;    
    ind_numeric = indices_base(indice);
    plot([0.5, 3.5], [ind_numeric, ind_numeric], 'r-', 'Color', [0.5, 0, 0.5], 'LineWidth', 2)
    
    
    xticklabels({' ', 'R1',' ', 'R2', ' ', 'R3'});
    title(['Indice ', num2str(indice)]);
    
    legend('Baixa Frequência', 'Média Frequência', 'Alta Frequência', 'Sem Estímulo');
    

end


%%

ritmos = {'Delta', 'Teta', 'Alfa', 'Beta'};

for ritmo_idx = 1:length(ritmos)
    ritmo = ritmos{ritmo_idx};
    
    % Use eval para criar variáveis com nomes dinâmicos
    eval([ritmo '_B = [', ritmo '_niveis_rodadas{1, 1}(1), ', ritmo '_niveis_rodadas{1, 2}(1), ', ritmo '_niveis_rodadas{1, 3}(1)];']);
    eval([ritmo '_M = [', ritmo '_niveis_rodadas{1, 1}(2), ', ritmo '_niveis_rodadas{1, 2}(2), ', ritmo '_niveis_rodadas{1, 3}(2)];']);
    eval([ritmo '_A = [', ritmo '_niveis_rodadas{1, 1}(3), ', ritmo '_niveis_rodadas{1, 2}(3), ', ritmo '_niveis_rodadas{1, 3}(3)];']);

    y_Baixa = eval([ritmo '_B']); 
    y_Media = eval([ritmo '_M']);
    y_Alta = eval([ritmo '_A']);

    nome_figura = sprintf('Figura_5_%d', ritmo_idx);
    figure('Name', nome_figura);
    % Gráfico de dispersão
    plot(y_Baixa, '-o');
    hold on;
    plot(y_Media, '-s');
    plot(y_Alta, '-d');  
    
    ritmo_base = ritmos_base(ritmo_idx);
    plot([0.5, 3.5], [ritmo_base, ritmo_base], 'r-', 'Color', [0.5, 0, 0.5], 'LineWidth', 2)

    xticks(1:3);
    xticklabels({'R1', 'R2', 'R3'});
    title([ritmo]);

    legend('Baixa Frequência', 'Média Frequência', 'Alta Frequência', 'Sem Estímulo');
  
end


