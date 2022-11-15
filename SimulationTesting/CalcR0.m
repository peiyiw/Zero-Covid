function [R0,seedStructure] = CalcR0(pop_stru_size,suscep,asym_prob,contact_matrix,beta,r1,r2,sigma,sigma_pre,gamma_a,gamma_s)

GroupSize = length(asym_prob);

suscep = reshape(suscep,GroupSize,1);
Beta_S = beta*contact_matrix.*suscep;
pop_stru_size = reshape(pop_stru_size,GroupSize,1);
% Beta_S = Beta_S.*pop_stru_size;
% Beta_S = Beta_S./pop_stru_size';

Beta_A = r1*Beta_S;
Beta_Lp = r2*Beta_S;

%%%%---E1,...,E7, IA1,....,IA7, Lp1,...,Lp7, IS1,...,IS7.
T = zeros(GroupSize*4);
T(1:GroupSize,(GroupSize+1):(2*GroupSize)) = Beta_A;
T(1:GroupSize,(2*GroupSize+1):(3*GroupSize)) = Beta_Lp;
T(1:GroupSize,(3*GroupSize+1):(4*GroupSize)) = Beta_S;

S = zeros(GroupSize*4);
S(1:GroupSize,1:GroupSize) = -sigma*eye(GroupSize);

S((GroupSize+1):(2*GroupSize),1:GroupSize) = sigma*diag(asym_prob);
S((GroupSize+1):(2*GroupSize),(GroupSize+1):(2*GroupSize)) = -gamma_a*eye(GroupSize);


S((2*GroupSize+1):(3*GroupSize),1:GroupSize) = sigma*diag(1-asym_prob);
S((2*GroupSize+1):(3*GroupSize),(2*GroupSize+1):(3*GroupSize)) = -sigma_pre*eye(GroupSize);

S((3*GroupSize+1):(4*GroupSize),(2*GroupSize+1):(3*GroupSize)) = sigma_pre*eye(GroupSize);
S((3*GroupSize+1):(4*GroupSize),(3*GroupSize+1):(4*GroupSize)) = -gamma_s*eye(GroupSize);

K = -T*inv(S);

[V,D]=eig(K);
[R0 i]=max(abs(diag(D))); 
seedStructure=abs(V(1:GroupSize,i));
seedStructure=seedStructure/sum(seedStructure);




