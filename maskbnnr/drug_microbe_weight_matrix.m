m = xlsread('Disease_Semantic_Similarity_Matrix.xlsx');
[m1,m2]=size(m);
weight_disease = zeros(m1,m2);      %weight_microbe 为病毒对应加权矩阵
for i=1:m1
    for j=1:m1
        if m(i,j)>0
            weight_disease(i,j)=1;
        end
    end
end
xlswrite('E:\L_B\code for MDHGI and BNNR\749-277\five_fold\代码\weight_disease.xlsx',weight_disease);


d1= xlsread('Protein_Sequence_Similarity_Matrix.xlsx');
[c1,c2]=size(d1);
weight_protein=zeros(c1,c2);
for i=1:c1
    for j=1:c1
        if d1(i,j)>0 
            weight_protein(i,j)=1;
        end
    end
end
xlswrite('E:\L_B\code for MDHGI and BNNR\749-277\five_fold\代码\weight_protein.xlsx',weight_protein); 
disp(99999);
    