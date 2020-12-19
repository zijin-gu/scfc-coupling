function gamma = generate_gamma(subject_num)
% compute the covariance matrix: gamma (TKT', TlambdaT', TT')
% compute T: subject_num*4 * subject_num
T = zeros(subject_num*4, subject_num);
n = 1;
for k = 1:4:4*subject_num
    T(k:k+3,n) = 1;
    n = n + 1;
end

% generate sample K
% MZ twin pairs - 1/2; DZ twin pairs - 1/4; siblings - 1/4;
kinship = zeros(subject_num,subject_num);
% randomly sample 1 MZ twin pair, 1 DZ twin pair
twin_idx = randi(subject_num, 2);
kinship(twin_idx(1,1),twin_idx(1,2)) = 1/2;
kinship(twin_idx(1,2),twin_idx(1,1)) = 1/2;
kinship(twin_idx(2,1),twin_idx(2,2)) = 1/4;
kinship(twin_idx(2,2),twin_idx(2,1)) = 1/4;

K = 2*kinship;
K(1:1+size(kinship,1):end) = 1;
% compute TKT'
tkt = T * K * T.';

% compute lambda
lambda = zeros(subject_num,subject_num);
lambda(twin_idx(1,1),twin_idx(1,2)) = 1;
lambda(twin_idx(1,2),twin_idx(1,1)) = 1;
lambda(twin_idx(2,1),twin_idx(2,2)) = 1;
lambda(twin_idx(2,2),twin_idx(2,1)) = 1;

lambda(1:1+size(lambda,1):end) = 1;
% compute TlambdaT'
tlt = T * lambda * T.';

% % compute gamma
gamma = zeros(4*subject_num, 4*subject_num, 3);
gamma(:,:,1) = tkt; 
gamma(:,:,2) = tlt;
gamma(:,:,3) = T * T.';