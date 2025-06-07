function [ drug_similarity,virus_similarity ] = integration_drug_virus_similarity(drug_gauss, virus_gauss)
%INT 此处显示有关此函数的摘要
%   此处显示详细说明

% wd = xlsread('weight_protein.xlsx');  %权重矩阵
% wv = xlsread('weight_disease.xlsx');
% ss = xlsread('Protein_Sequence_Similarity_Matrix.xlsx'); %药物化学结构相似性矩阵
% 
% vs = xlsread('Disease_Semantic_Similarity_Matrix.xlsx');          % 病毒序列相似性
% 
% drug_similarity = wd.*ss+(~wd).*drug_gauss;     %药物最终相似性矩阵
% 
% virus_similarity= wv.*vs+(~wv).*virus_gauss;  %病毒最终相似性矩阵




% wd = xlsread('weight_protein.xlsx');  %权重矩阵
% wv = xlsread('weight_disease.xlsx');
% ss = xlsread('pp_decline.xlsx'); %药物化学结构相似性矩阵
% 
% vs = xlsread('disease_decline.xlsx');          % 病毒序列相似性
% 
% drug_similarity = wd.*ss+(~wd).*drug_gauss;     %药物最终相似性矩阵
% 
% virus_similarity= wv.*vs+(~wv).*virus_gauss;  %病毒最终相似性矩阵




% wd = xlsread('weight_protein.xlsx');  %权重矩阵
% wv = xlsread('weight_disease.xlsx');
ss = xlsread('Protein_Sequence_Similarity_Matrix.xlsx'); %药物化学结构相似性矩阵

vs = xlsread('Disease_Semantic_Similarity_Matrix.xlsx');          % 病毒序列相似性

drug_similarity = ss+drug_gauss;     %药物最终相似性矩阵

virus_similarity= vs+virus_gauss;  %病毒最终相似性矩阵
end

