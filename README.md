# SC-FC coupling

This is the analysis scripts of SC-FC coupling for paper [Heritability and interindividual variability of regional structure-function coupling](https://www.nature.com/articles/s41467-021-25184-4).
All data is created randomly for demonstration purpose.

## Requirements
1) Install Python by following instructions here: https://www.python.org/downloads/
2) Install Matlab by following instructions here: https://www.mathworks.com/help/install/install-products.html
3) `pip -r install requirements.txt`

## Steps
1) At your terminal, run `glm_demo.py` to get the SC-FC coupling and GLM results saved in data directory and glm_res directory.
2) In your MATLAB, run `herit_demo.m `to get the heritability estimates of some random phenotypes you created.

## Expected output
1) Under the data directory, you will have your self-generated SC, FC and SC-FC coupling data for the number of subjects you specify. Also there will be covariates data (age, sex, cognition, ICV and motion) if you do GLM analysis.
2) Under the glm_res directory, you will have a series of txt files named  `GLMvar?prec/corr_all/sigroi.txt` which contains the -log(p) values for different covariates (index by number ?). `allroi` files contains -log(p) value for all the regions while `sigroi` files only contains regions that are significant with other regions masked as 0. There is another mat file named GLM_'corr/prec'_confint.mat which contains the confident interval for the GLM results.
Under the herit_res directory, you will have a file named `herit_cp.mat` which is the heritability estimates using the LME model.

### References
Gu, Z., Jamison, K. W., Sabuncu, M. R., & Kuceyeski, A. (2021). Heritability and interindividual variability of regional structure-function coupling. Nature Communications, 12(1), 1-12.

Ge, T., Holmes, A. J., Buckner, R. L., Smoller, J. W., & Sabuncu, M. R. (2017). Heritability analysis with repeat measurements and its application to resting-state functional connectivity. Proceedings of the National Academy of Sciences, 114(21), 5521-5526.


