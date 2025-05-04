# synco: Differences between synonymous codon-specific Ramachandran plots

### Code description

This code reproduces the analyses presented in [[2]](https://www.biorxiv.org/content/10.1101/2022.11.29.518303v3), showing that differences between the codon-specific Ramachandran plots provided in [[3]](https://doi.org/10.1038/s41467-022-30390-9) are not statistically significant. The two-sample goodness-of-fit tests that have been used for the analyses were introduced in [[1]](https://doi.org/10.1214/23-EJS2135), whose implementation is available in the R package [torustest](https://github.com/gonzalez-delgado/torustest).

* The file [codon_test.R](https://github.com/gonzalez-delgado/synco/blob/main/codon_test.R) performs the goodness-of-fit test assessing the null hypothesis "Amino-acid backbone distribution does not depend on the identity of the translated codon". Codon-specific Ramachandran plots are defined for a given amino-acid without taking into account nearest neighbors effects. The analysis is carried out on the [experimental dataset](https://doi.org/10.7910/DVN/5P81D4) provided in [[3]](https://doi.org/10.1038/s41467-022-30390-9) and on the [same structures extracted from the AlphaFold database](https://zenodo.org/doi/10.5281/zenodo.11110092). Instructions to perform the analyses are provided in the file, that produces Figures 2a, 2b and results described in Section D of the SI in [[2]](https://www.biorxiv.org/content/10.1101/2022.11.29.518303v3).

* The same null hypothesis, but taking neighboring residues into account, is tested in [codon_test_tripeptides.R](https://github.com/gonzalez-delgado/synco/blob/main/codon_test_tripeptides.R). Here, the codon-specific Ramachandran plots contain conformations coming from a fixed tripeptide (fragment of three consecutive amino-acids). This analysis uses the function [get_tripeptides](https://github.com/gonzalez-delgado/synco/blob/main/get_tripeptides.R), which uses the sequence information of the [dataset](https://doi.org/10.7910/DVN/5P81D4) provided in [[3]](https://doi.org/10.1038/s41467-022-30390-9) to extract neighbors' identities. The code reproduces the analysis presented in Section E of the SI in [[2]](https://www.biorxiv.org/content/10.1101/2022.11.29.518303v3).

* The simulation study described in Section B of the SI in [[2]](https://www.biorxiv.org/content/10.1101/2022.11.29.518303v3) is reproduced in [simulation_study.R](https://github.com/gonzalez-delgado/synco/blob/main/simulation_study.R). This files produces also Figure 1 in [[2]](https://www.biorxiv.org/content/10.1101/2022.11.29.518303v3).

For any inquires, please file an [issue](https://github.com/gonzalez-delgado/synco/issues) or [contact us](mailto:javier.gonzalezdelgado@mcgill.ca).

### Related work

Check [Akeju & Cope paper](https://doi.org/10.1093/gbe/evae080) for a complementary analysis of the codon and bond angles correlation described in [Rosenberg _et al._ 2022](https://doi.org/10.1038/s41467-022-30390-9).

### References

[1] González-Delgado, J., González-Sanz, A., Cortés, J., & Neuvial, P. (2023). Two-sample goodness-of-fit tests on the flat torus based on Wasserstein distance and their relevance to structural biology. <i>Electron. J. Statist</i>., 17(1): 1547–1586. [https://doi.org/10.1214/23-EJS2135](https://doi.org/10.1214/23-EJS2135).

[2] González-Delgado, J., Bernadó, B., Mier, P., Neuvial, P. & Cortés, J. (2024). The dependence of the amino acid backbone conformation on the translated synonymous codon is not statistically significant. Submitted. [bioRxiv](https://www.biorxiv.org/content/10.1101/2022.11.29.518303v3).

[3] Rosenberg, A.A., Marx, A. & Bronstein, A.M. (2022). Codon-specific Ramachandran plots show amino acid backbone conformation depends on identity of the translated codon. <i>Nat Commun</i>, 13, 2815. [https://doi.org/10.1038/s41467-022-30390-9](https://doi.org/10.1038/s41467-022-30390-9).

[4] Akeju, O.J., Cope, A. L. (2024). Re-examining Correlations Between Synonymous Codon Usage and Protein Bond Angles in Escherichia coli. <i>Genome Biology and Evolution</i>, 16(5): evae080. [https://doi.org/10.1093/gbe/evae080](https://doi.org/10.1093/gbe/evae080).





