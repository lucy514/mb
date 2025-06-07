function predict_matrix=ensemble_average(predict_matrix1,predict_matrix2,train_index_matrix,w) %集成方式的平均值
%addpath('D:\matlab2016\side effects\source code and datasets\RBMLIB\RBM');
temp1=(predict_matrix1(:));
mean_value1=mean(temp1(train_index_matrix(:)));
std_value1=std(temp1(train_index_matrix(:)));
temp2=(predict_matrix2(:));
mean_value2=mean(temp2(train_index_matrix(:)));
std_value2=std(temp2(train_index_matrix(:)));
% temp3=(predict_matrix3(:));
% mean_value3=mean(temp3(train_index_matrix(:)));
% std_value3=std(temp3(train_index_matrix(:)));
predict_matrix=w*normalize_01(predict_matrix1,mean_value1, std_value1)+(1-w)*normalize_01(predict_matrix2,mean_value2, std_value2);
end

function normalize_score=normalize_01(score,mean_value,std_value)
normalize_score=0.5*(tanh(0.1*(score-mean_value)./std_value)+1);
end
