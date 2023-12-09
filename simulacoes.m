clear all, clc, close all

load('G:\Funcoes_TCC\fases\F1B\B1RG1.mat')
F1B=X;

load('G:\Funcoes_TCC\fases\F1M\M1RG1.mat')
F1M=X;

load('G:\Funcoes_TCC\fases\F1A\A1RG1.mat')
F1A=X;

simulacao_1=[F1B,F1M,F1A];

%%

load('G:\Funcoes_TCC\fases\F2B\B2RG1.mat')
F2B=X;

load('G:\Funcoes_TCC\fases\F2M\M2RG1.mat')
F2M=X;

load('G:\Funcoes_TCC\fases\F2A\A2RG1.mat')
F2A=X;

simulacao_2=[F2B,F2M,F2A];


%%

load('G:\Funcoes_TCC\fases\F3B\B3RG1.mat')
F3B=X;

load('G:\Funcoes_TCC\fases\F3M\M3RG1.mat')
F3M=X;

load('G:\Funcoes_TCC\fases\F3A\A3RG1.mat')
F3A=X;

simulacao_3=[F3B,F3M,F3A];

%
simulacao_4= [simulacao_1;simulacao_2;simulacao_3];

%%
% 
% simulacoes =cell(1, 4);
% simulacoes{1,1}=simulacao_1;
% simulacoes{1,2}=simulacao_2;
% simulacoes{1,3}=simulacao_3;
% simulacoes{1,4}=simulacao_4;

