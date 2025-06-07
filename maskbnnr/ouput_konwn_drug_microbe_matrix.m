A=load('adj.txt'); 
nd=max(A(:,1));                                                           % nd:the number of drug
nm=max(A(:,2));                                                           % nv:the number of microbe
[pp,qq]=size(A);                                                          % pp:the number of known drug-microbe associations,pp=2884,qq=3
interaction=zeros(nd,nm);   
for i=1:pp
        interaction(A(i,1),A(i,2))=1;                                     %interaction: adajency matrix for the drug-microbe association network,已知联系的置1
end
xlswrite('D:\matlab2016\173mic-drug\interaction.xlsx',interaction); 
disp(999);