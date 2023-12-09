clear all, clc, close all
addpath('fases');
pasta='G:\Funcoes_TCC\dadosArg\Sinais\';

vols = [4 6 10 12];

% B3RG1
k=1;
for k=1:length(vols)
    load([pasta 'Vol_' num2str(vols(k)) '_LOW.mat']);
    for q = 1:length(ans)
        X{k,q} = ans{1,q};
    end
    k=k+1;
end


%%

%  for k=1:38
%      X{3,k}.signal = [];
%  end
