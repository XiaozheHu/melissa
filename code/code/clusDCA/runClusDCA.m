% type 1: MF
% type 2: BP
% type 3: CC
% cat 1: 3-10
% cat 2: 11-30
% cat 3: 31-100
% cat 4 101-300

specie='Yeast';
path = '/Users/kaiyiwu/Desktop/clusDCA-master';
addpath(genpath('Code'));
addpath(genpath('Data'));
 

dim=2500;    % embedding dimension
rspx = 0.5;  % RWR restart probability for x
rspy = 0.8;  % RWR restart probability for y
bp = 0.8;  % back propagation parameter alpha
nnode = 6311;%16662;
nlabel = 4240;%13807;

%16662 node for human, 6311 node for yeast

% compute the low dimensional embedding for x(gene) and y(labels)
[lx, ly] = learnVector(specie, dim, bp, rspx, rspy,nnode);
% store this

our_score = clusDCA(specie,lx,ly);

%

% dca_score = DCA(specie,lx,ly);

