clear;
cwd = pwd;
data_dir = strcat(cwd, '/data/');
herit_dir = strcat(cwd, '/herit_res/');
if ~exist(herit_dir, 'dir')
       mkdir(herit_dir)
end

subject_num = 1000; % subject to change
region_num = 392;

% create your phenotype data with 4 sessions
cp1 = randn(subject_num, region_num);
cp2 = randn(subject_num, region_num);
cp3 = randn(subject_num, region_num);
cp4 = randn(subject_num, region_num);

% create your covariates: age, sex and handedness
age = randi([22,37], subject_num);

sex = ones(subject_num);
male_num = randi(subject_num);
female_num = subject_num - male_num;
sex(randperm(subject_num, female_num)) = -1;

handedness = randi([-100,100], subject_num);

% prepare the covariates for the model
covariates = zeros(subject_num*4, 3);
n = 1;
for i = 1:4:subject_num*4
    covariates(i:i+3,1) = age(n);
    covariates(i:i+3,2) = sex(n);
    covariates(i:i+3,3) = handedness(n);
    n = n + 1;
end

% prepare the phenotype data for the model
coupling = zeros(region_num,subject_num*4);
for k = 1:region_num
    n = 1;
    for i = 1:4:subject_num*4
        coupling(k,i:i+3) = [cp1(n,k),cp2(n,k),cp3(n,k),cp4(n,k)];
        n = n + 1;
    end
end

% generate your covariance matrix: gamma (TKT', TlambdaT', TT')
gamma = generate_gamma(subject_num);

X = covariates;
res = [];
for i = 1:region_num
    disp(strcat('---------- Region', int2str(i),' ----------'))
    y = coupling(i,:).';
    [res(i).flag, res(i).m2_tot, res(i).se_tot, res(i).M2, res(i).SE, res(i).Vc, res(i).Lnew] = Morphometricity_NRandEffects(y, X, gamma, 0, 1e-4, 100, 'True');
end
save(strcat(herit_dir, 'herit_cp.mat'),'res');


