# synco: Differences between synonymous codon-specific Ramachandran plots

The code presented in synco implements the two-sample goodness-of-fit test introduced in [1] to assess the effect of the translated codon on amino-acid backbone conformation. This work, described in [2], reproduces the analysis presented in [3] using the methodology defined in [1]. 

The file [codon_test.R](https://github.com/gonzalez-delgado/synco/blob/main/codon_test.R) performs the goodness-of-fit test assessing the null hypothesis "the codon-specific (phi, psi) distribution is the same for synonymous codons". Codon-specific Ramachandran plots are defined for a given amino-acid without taking into account the nearest neighbors effect. The same null hypothesis, but taking neighboring residues into account, is tested in [codon_test_tripeptides.R](https://github.com/gonzalez-delgado/synco/blob/main/codon_test_tripeptides.R). Here, the codon-specific Ramachandran plots contain conformations coming from a fixed tripeptide (fragment of three consecutive amino-acids). This analysis uses the function [get_tripeptides](https://github.com/gonzalez-delgado/synco/blob/main/get_tripeptides.R), which uses the sequence information of the [database](https://doi.org/10.7910/DVN/5P81D4) to extract neihgbors' identities.

#### References

[1] González-Delgado, J., González-Sanz, A., Cortés, J., & Neuvial, P. (2021). Two-sample goodness-of-fit tests on the flat torus based on Wasserstein distance and their relevance to structural biology. <i>Electron. J. Statist</i>., 17(1): 1547–1586, 2023. [[url]](https://doi.org/10.1214/23-EJS2135) [[code]](https://github.com/gonzalez-delgado/torustest).

[2] J. González-Delgado, Pau Bernadó, Pablo Mier, Pierre Neuvial, Juan Cortés. Statistical tests to detect differences between codon-specific Ramachandran plots. Submitted. [[bioRxiv]](https://www.biorxiv.org/content/10.1101/2022.11.29.518303v1).

[3] Rosenberg, A.A., Marx, A. & Bronstein, A.M. Codon-specific Ramachandran plots show amino acid backbone conformation depends on identity of the translated codon. Nat Commun 13, 2815 (2022). https://doi.org/10.1038/s41467-022-30390-9.

