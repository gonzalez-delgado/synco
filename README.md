# synco: Differences between synonymous codon-specific Ramachandran plots

This code reproduces the analyses presented in [González-Delgado _et al._ 2024](https://www.biorxiv.org/content/10.1101/2022.11.29.518303v2), showing that differences between the codon-specific Ramachandran plots provided in [Rosenberg _et al._ 2022](https://doi.org/10.1038/s41467-022-30390-9) are statistically non-significant. The two-sample goodness-of-fit tests that have been used for the analyses were introduced in [González-Delgado _et al._ 2023](https://doi.org/10.1214/23-EJS2135), whose implementation is available in the R package [torustest](https://github.com/gonzalez-delgado/torustest).

* The file [codon_test.R](https://github.com/gonzalez-delgado/synco/blob/main/codon_test.R) performs the goodness-of-fit test assessing the null hypothesis "Amino-acid backbone distribution does not depend on the identity of the translated codon". Codon-specific Ramachandran plots are defined for a given amino-acid without taking into account nearest neighbors effects. The analysis is carried out on the [experimental database](https://doi.org/10.7910/DVN/5P81D4) provided in [Rosenberg _et al._ 2022](https://doi.org/10.1038/s41467-022-30390-9) and on its [AlphaFold2 counterpart](https://github.com/gonzalez-delgado/synco/blob/main/nt_structure_2024.txt). Instructions to perform the analyses are provided in the file, that produces Figures 2a, 2b and S4 in [González-Delgado _et al._ 2024](https://www.biorxiv.org/content/10.1101/2022.11.29.518303v2).

* The same null hypothesis, but taking neighboring residues into account, is tested in [codon_test_tripeptides.R](https://github.com/gonzalez-delgado/synco/blob/main/codon_test_tripeptides.R). Here, the codon-specific Ramachandran plots contain conformations coming from a fixed tripeptide (fragment of three consecutive amino-acids). This analysis uses the function [get_tripeptides](https://github.com/gonzalez-delgado/synco/blob/main/get_tripeptides.R), which uses the sequence information of the [database](https://doi.org/10.7910/DVN/5P81D4) provided in [Rosenberg _et al._ 2022](https://doi.org/10.1038/s41467-022-30390-9) to extract neihgbors' identities.

* The simulation study described [González-Delgado _et al._ 2024](https://www.biorxiv.org/content/10.1101/2022.11.29.518303v2) (Section B of the SI) is reproduced in [simulation_study.R](https://github.com/gonzalez-delgado/synco/blob/main/simulation_study.R). This files produces Figures 1 and S1 in [González-Delgado _et al._ 2024](https://www.biorxiv.org/content/10.1101/2022.11.29.518303v2).

If you encounter any trouble to run the scripts or reproduce the results presented in [1], please file an [issue](https://github.com/gonzalez-delgado/synco/-/issues) or [contact us](mailto:javier.gonzalezdelgado@mcgill.ca).

#### References

[1] González-Delgado, J., González-Sanz, A., Cortés, J., & Neuvial, P. (2021). Two-sample goodness-of-fit tests on the flat torus based on Wasserstein distance and their relevance to structural biology. <i>Electron. J. Statist</i>., 17(1): 1547–1586, 2023. [https://doi.org/10.1214/23-EJS2135](https://doi.org/10.1214/23-EJS2135).

[2] J. González-Delgado, Pau Bernadó, Pablo Mier, Pierre Neuvial, Juan Cortés. Statistical tests to detect differences between codon-specific Ramachandran plots. Submitted. [bioRxiv](https://www.biorxiv.org/content/10.1101/2022.11.29.518303v2).

[3] Rosenberg, A.A., Marx, A. & Bronstein, A.M. Codon-specific Ramachandran plots show amino acid backbone conformation depends on identity of the translated codon. Nat Commun 13, 2815 (2022). [https://doi.org/10.1038/s41467-022-30390-9](https://doi.org/10.1038/s41467-022-30390-9).


