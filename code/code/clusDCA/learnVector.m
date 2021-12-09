function [ USA,USB] = learnVector(specie, d, bp, rspx, rspy,nnode)
% load network data
%
% [Input]
% weight: 0 represents the all the edges have same weight. 1 represents all
% the deges have different weight
% nnode: number of nodes
% bp: back propogation probability
%
% [Output]
% A : the sparse matrix of network.

go_network = ['Data/',specie,'Graph/go.txt'];

% get six ppi networks
for i=1:6
    file_name = ['Data/',specie,'Graph/ppi',num2str(i),'.txt'];
    nt = load_network(file_name,true,nnode);
    nnode = size(nt,1);
    tA = run_diffusion(nt, 'personalized-pagerank', struct('maxiter', 20, 'reset_prob', rspx));
    if i==1
        QA = tA;         
        continue
    end
    QA = [QA,tA];
 
end

alpha = 1/nnode;
QA = log(QA+alpha)-log(alpha);
QA=QA*QA';


% get label network
% label network loaded with back propagation probability built-in
B = load_network(go_network,false,0,bp);
nlabel=size(B,1);
QB = run_diffusion(B, 'personalized-pagerank', struct('maxiter', 20, 'reset_prob', rspy));
alpha = 1/nlabel;
QB = log(QB+alpha)-log(alpha);


% get low dimensional embedding for x and y
fprintf('run X SVD d=%d\n',d); 
QA = sparse(QA);
[U,S] = svds(QA,d);
LA = U;
USA = LA*sqrt(sqrt(S)); 


fprintf('run Y SVD d=%d\n',d); 
[U,S] = svds(QB,d);
LB = U;
USB = LB*sqrt(S); 

end

