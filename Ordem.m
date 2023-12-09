clear all, clc, close all
addpath('fases');
pasta='G:\Funcoes_TCC\fases\MBA\';

vols = [3 9 15 16];

B=[];
M=[];
A=[];

for i=1:3
    k=1;
    if i==1
        for k=1:length(vols)
            load([pasta 'Vol_' num2str(vols(k)) '_LOW.mat']);
            for q = 1:length(ans)
                B{k,q} = ans{1,q};
            end
            k=k+1;
        end
    end
      
    k=1;
    if i==2
        for k=1:length(vols)
            load([pasta 'Vol_' num2str(vols(k)) '_MEDIUM_A.mat']);
            for q = 1:length(ans)
                M{k,q} = ans{1,q};
            end
            k=k+1;
        end
    end
    
    k=1;
    if i==3
        for k=1:length(vols)
            load([pasta 'Vol_' num2str(vols(k)) '_HIGH_A.mat']);
            for q = 1:length(ans)
                A{k,q} = ans{1,q};
            end
            k=k+1;
        end
    end
    
    ordem = [B,M,A];
    i=i+1;
end


%%
