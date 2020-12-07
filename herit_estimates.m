%function herit_estimates(phenotype, fctype)
data_dir = './data/heritability/';
subj_id = csvread(strcat(data_dir, 'herit_subject941.txt'));
covariates = csvread(strcat(data_dir, 'covariates.csv'));
subj_num = length(covariates)/4;

roi_number = 392;

data_prename = strcat(data_dir, 'herit2model_', phenotype, '_', fctype, '_session');
% load the sc-fc coupling
coupling = zeros(roi_number,subj_num*4);

cp1 = load(strcat(data_prename, '1.mat'));
cp2 = load(strcat(data_prename, '2.mat'));
cp3 = load(strcat(data_prename, '3.mat'));
cp4 = load(strcat(data_prename, '4.mat'));
cp1 = cp1.type;
cp2 = cp2.type;
cp3 = cp3.type;
cp4 = cp4.type;

% subj993 = csvread(strcat(data_dir, 'herit_subject993.txt'));
% 
% [sharedvals,idx] = intersect(subj993, subj_id, 'stable');
% cp1 = cp1(idx,:);
% cp2 = cp2(idx,:);
% cp3 = cp3(idx,:);
% cp4 = cp4(idx,:);
for k = 1:392
    n = 1;
    for i = 1:4:length(subj_id)*4
        coupling(k,i:i+3) = [cp1(n,k),cp2(n,k),cp3(n,k),cp4(n,k)];
        n = n + 1;
    end
end

cov_mat = load(strcat(data_dir,'herit_covar.mat'));
gamma = cov_mat.gamma;
X = covariates;
res = [];
for i = 1:392
    disp(strcat('---------- Region', int2str(i),' ----------'))
    y = coupling(i,:).';
    [res(i).flag, res(i).m2_tot, res(i).se_tot, res(i).M2, res(i).SE, res(i).Vc, res(i).Lnew] = Morphometricity_NRandEffects(y, X, gamma, 0, 1e-4, 100, 'True');
end
save(strcat('./results/herit_', phenotype, '_', fctype,'.mat'),'res');
