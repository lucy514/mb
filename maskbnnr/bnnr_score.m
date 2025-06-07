A=load('Protein_Disease_adj.txt');                                                       %注意  RBM 代表 MDHGI, INMB 代表 BNNR
nd=max(A(:,1));                                                           % nd:the number of drug
nv=max(A(:,2));                                                           % nv:the number of virus
[pp,qq]=size(A);                                                          % pp:the number of known drug-virus associations,pp=933,qq=3

interaction_matrix = readmatrix('F:\BNNR\run_maskBNNR\749-276\five_fold\code\Protein_Disease_Associations.xlsx');
[row,col]=size(interaction_matrix);
number=100;                                                               %define the parameters

C1_INBM=[];


five_k_INBM=zeros(number,pp);


% BNNR 模型的一些参数
maxiter = 300;
alpha = 1;
beta = 10;
tol1 = 2*1e-3;
tol2 = 1*1e-5;

%gNN is used to predict the association between protein and disease.
array1= readmatrix('array_P-D.xlsx');
tic;
for k= 91:100
    fprintf('第%d轮\n',k);
    current_time = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss');
    disp(current_time);
    array=array1(k,:);                                      %array是100*933的矩阵，每一行是1-933的随机序列
    for circle=0:4
        if circle<4                                           %因为933除以5除不尽，所以前四折是187，最后一折为185
            nu=181;
            interaction_matrix1=interaction_matrix;
            interaction_matrix1_index=logical(zeros(row,col)~=0);    %interaction_matrix1_index为训练样本对应位置的矩阵
            for i=1:pp
                interaction_matrix1_index(A(i,1),A(i,2))=1;
            end
            new_array = array(1,1+circle*181:(circle+1)*181);      % 933个药物-病毒关联对，前面四个nu=187，最后一个n=185
            
            for j=1:181                                            %j=1:nu
                o=A(new_array(1,j),1);
                l=A(new_array(1,j),2);
                interaction_matrix1(o,l)=0;                       % 依次置将1变为0
                interaction_matrix1_index(o,l)=0;
            end
            
            drug_gauss = readmatrix('Protein_Gaussian_Similarity_Matrix.xlsx');
            virus_gauss = readmatrix('Disease_Gaussian_Similarity_Matrix.xlsx');
            [drug_integration_similarity,virus_integration_similarity] = integration_drug_virus_similarity(drug_gauss, virus_gauss ); %调用函数获得整合药物副作用，药物化学结构相似性，高斯相似性的药物相似性矩阵，病毒相似性也是如此
            
            
            %BNNR 预测的关联得分
            Wdd = drug_integration_similarity;
            Wvv= virus_integration_similarity;
            Wvd = interaction_matrix1';
            [dn,dr] = size(Wvd);

            T = [Wdd, Wvd'; Wvd, Wvv];
            [t1, t2] = size(T);
            trIndex = double(T ~= 0);
            [WW,iter] = BNNR(alpha, beta, T, trIndex, tol1, tol2, maxiter, 0, 1);
            M_recovery = (WW((t1-dn+1) : t1, 1 : dr))';
            A1 = M_recovery;
            A_Opposite = (1-interaction_matrix);
            A11 = A1.*A_Opposite;
            col_sum = sum(A11); % 计算每列的总和
            A11_norm = A11 ./ col_sum; % 归一化，使每列的和为1
            A_M = interaction_matrix.*A1;
            A_MASK = A11_norm+A_M;
            predict_score_matrix2 = A_MASK;
            
            
            
            %预测排名：INBM方法
            Sco2= predict_score_matrix2;
            final_score=Sco2( interaction_matrix==0);
            for i=1:181
                q=A(new_array(1,i),1);
                w=A(new_array(1,i),2);
                s_score = Sco2(q, w);
                five_score(i) = s_score;

            end
            
           
          
            
           C1_INBM = [C1_INBM, five_score];
            
            
        end
        
        if circle==4
            nu=181;
            interaction_matrix1=interaction_matrix;
            interaction_matrix1_index=logical(zeros(row,col)~=0);
            for i=1:pp
                interaction_matrix1_index(A(i,1),A(i,2))=1;
            end
            new_array = array(1,1+circle*181:pp);
            for j=1:181
                o=A(new_array(1,j),1);
                l=A(new_array(1,j),2);
                interaction_matrix1(o,l)=0;
                interaction_matrix1_index(o,l)=0;
            end
            drug_gauss = similarity_drug(interaction_matrix1);
            virus_gauss = similarity_microbe(interaction_matrix1);
            [drug_integration_similarity,virus_integration_similarity] = integration_drug_virus_similarity(drug_gauss, virus_gauss ); %调用函数获得整合药物副作用，药物化学结构相似性，高斯相似性的药物相似性矩阵，病毒相似性也是如此

            
            %BNNR 预测的关联得分
            Wdd = drug_integration_similarity;
            Wvv= virus_integration_similarity;
            Wvd = interaction_matrix1';
            [dn,dr] = size(Wvd);

            T = [Wdd, Wvd'; Wvd, Wvv];
            [t1, t2] = size(T);
            trIndex = double(T ~= 0);
            [WW,iter] = BNNR(alpha, beta, T, trIndex, tol1, tol2, maxiter, 0, 1);
            M_recovery = (WW((t1-dn+1) : t1, 1 : dr))';
            A1 = M_recovery;
            A_Opposite = (1-interaction_matrix);
            A11 = A1.*A_Opposite;
            col_sum = sum(A11); % 计算每列的总和
            A11_norm = A11 ./ col_sum; % 归一化，使每列的和为1
            A_M = interaction_matrix.*A1;
            A_MASK = A11_norm+A_M;
            predict_score_matrix2 = A_MASK;
            
           
            
            %第二种：INBM方法
            Sco2= predict_score_matrix2;
            final_score=Sco2( interaction_matrix==0);
            for i=1:181
                q=A(new_array(1,i),1);
                w=A(new_array(1,i),2);
                s_score = Sco2(q, w);
                five_score_1(i) = s_score;

            end
            
          
         C1_INBM = [C1_INBM, five_score_1];

           
        end
    end
    
    
   
 
   five_k_INBM(k,:) = C1_INBM;
   
    C1_INBM=[];
    
    
end

xlswrite('F:\BNNR\run_maskBNNR\749-276\five_fold\result_maskbnnr\result_91-100.xlsx',five_k_INBM);



toc;
fprintf('five_fold结束\n');


