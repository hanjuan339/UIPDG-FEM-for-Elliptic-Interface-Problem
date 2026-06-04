a=[1 1000];
theta=[0 0];
P=1;
N=16;   
beta=1;
gamma0=100;
gamma1=100;
fname=['ex1data_' num2str(a(1)) '_' num2str(a(2)) '_'];
tic

for p=P
    for n=N
        fprintf('\n p=%g,  n=%g.\n',p,n);
        [err1,err0,errinf,errflux,u,el1,el2,eli,ltg,cn,en1_2,en1_4,en2_2,en2_4]=ex1fun(a,theta,n,p,beta,gamma0,gamma1);
        save([fname 'n' num2str(n) '_p' num2str(p) '.mat'],...
             'a','n','p','beta','gamma0',...
             'err1','err0','errinf','errflux','u','el1','el2','eli','ltg','cn','en1_2','en1_4','en2_2','en2_4');
    end
end

% a=[1 1000];
% theta=[0 0];
% P=1:3;
% N=[16 32 64 128 256];  
% beta=1;
% gamma0=100;
% gamma1=100;
% fname=['ex1data_' num2str(a(1)) '_' num2str(a(2)) '_'];
% tic
% 
% for p=P
%     for n=N
%         fprintf('\n p=%g,  n=%g.\n',p,n);
%         [err1,err0,errinf,errflux,u,el1,el2,eli,ltg,cn,en1_2,en1_4,en2_2,en2_4]=ex1fun(a,theta,n,p,beta,gamma0,gamma1);
%         save([fname 'n' num2str(n) '_p' num2str(p) '.mat'],...
%              'a','n','p','beta','gamma0',...
%              'err1','err0','errinf','errflux','u','el1','el2','eli','ltg','cn','en1_2','en1_4','en2_2','en2_4');
%     end
% end
% 
% 
% toc
% disp(['运行时间: ',num2str(toc)]);