function KD=similarity_drug(interaction)
[nd,~]=size(interaction);
for i=1:nd
sh(i)=norm(interaction(i,:))^2;                                              %calculate gamad for Gaussian kernel calculation
end
    gamad=nd/sum(sh');
    for i=1:nd
        for j=1:nd
    kh(i,j)=exp(-gamad*(norm(interaction(i,:)-interaction(j,:)))^2);        %calculate gamad for Gaussian kernel calculation
        end
    end
    KD=kh;

% KD=kh;                                              %the integrated similarity between disease



%      for i=1:nd
%        for j=1:nd
%           KD(i,j)=V1(i,j)./sqrt(sum(V1(i,:))*sum(V1(j,:)));
% end
%      end
%       
   
