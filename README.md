To consider:
* To stay consistent with GTEx analyses, I have decided to not merge the Omni samples
* If/when we want to merge the Omni samples, consider that the bplink files for autosomes and chrX are in two separate bplink filesets; you may want to merge these first
* The covariates can be left as is as having the omni covariates won't affect the regression
* All covariates will be added to covariates file. If you'd like to remove one, simply do so manually as each covariate is a row; this is easily done with vim.
