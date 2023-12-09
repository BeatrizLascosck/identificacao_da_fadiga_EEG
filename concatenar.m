clear all, clc, close all
% 
% load('G:\Funcoes_TCC\fases\B_all.mat')
% B=all_B;
% load('G:\Funcoes_TCC\fases\M_all.mat')
% M=all_M;
% load('G:\Funcoes_TCC\fases\A_all.mat')
% A=all_A;
% 
% estimulos_BMA=[B,M,A];


%%
dados_base =cell(20, 1);


diretorio = 'G:\Funcoes_TCC\dadosArg\Sinal Base\'; 

dados_concatenados = [];

for v = 1:20
    
    nome_arquivo = sprintf('vol_%d_A_BASE.mat', v);
    caminho_completo = fullfile(diretorio, nome_arquivo);
    
    load(caminho_completo);
    dados_base{v,1}.signal = data;
     
end
