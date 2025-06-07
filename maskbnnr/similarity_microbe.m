function KV=similarity_microbe(interaction)
  [~,nm]=size(interaction);
   for i=1:nm
sv(i)=norm(interaction(:,i))^2;
   end
    gamam=nm/sum(sv');
    for i=1:nm
        for j=1:nm
    kv(i,j)=exp(-gamam*(norm(interaction(:,i)-interaction(:,j)))^2);        %calculate Gaussian kernel for the similarity between miRNA: km
       end
    end 
     KV=kv;                                                    %the integrated similarity between miRNA 
