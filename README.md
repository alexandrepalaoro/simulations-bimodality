# Sigmoid allometries generate male dimorphism in secondary sexual traits: a comment on Packard (2023).

Authors: Bruno Buzatto, Glauco Machado, Alexandre V. Palaoro <br>
Journal: Evolutionary Ecology, DOI: 10.1007/s10682-024-10303-6 (<a href="https://link.springer.com/article/10.1007/s10682-024-10303-6" target="_blank">link</a>) <br>
Contact about code, and analyses: alexandre.palaoro@gmail.com

[![Zenodo](https://img.shields.io/badge/Zenodo-archive%20pending-lightgrey)](#archiving-on-zenodo)
[![paper](https://img.shields.io/badge/paper-10.1007%2Fs10682--024--10303--6-blue)](https://doi.org/10.1007/s10682-024-10303-6)
[![code license](https://img.shields.io/badge/code%20license-MIT-green)](LICENSE)
[![data license](https://img.shields.io/badge/data%20license-CC0%201.0-brightgreen)](LICENSE-DATA)

> [!IMPORTANT]
> **This repository is not archived on Zenodo yet.** Follow the three steps under
> *[Archiving on Zenodo](#archiving-on-zenodo)* at the bottom of this file. Then replace
> `ZENODO_DOI_HERE` below with the concept DOI Zenodo gives you, uncomment the `doi:` line
> in `CITATION.cff`, and delete this notice.

---

### This readme has been divided in three parts. First, we will talk about file structure, then the code, and then the dataset.

##### File structure:

We uploaded four different folders: "code", "data", "figures", and "simul_data". Each folder contain the files we used to run all analyses and the output files that came from the codes (the figures and simulation data). 

##### Code:

The folder "code" contains all code used in this paper. The code is divided in files used to run the analyses of the real data, and files used for the simulations. Codes that start with "Function" are required to run the analyses used in the "Simulating sigmoid allometries and detecting dimorphisms.R". These are just some functions to make it easier to analyze everything. The two code files that start with "simulation_" contain the raw code used for the simulations themselves. If you have the folder "simul_data", you do not need to use these codes and can use the resulting simulations directly with the main code: "Simulating sigmoid allometries and detecting dimorphisms.R"

##### Data:

Our data is divided in the folder "data" and the folder "simul_data".
In the "data" folder, we have the real harvestmen data used in the original papers. The name of the file relates to the species, and inside the folder you will have two columns, one related to the size of the body and another to the size of the weapon. 

In the "simul_data" you will find .RDS files that contains the simulations. You can run your own simulations, but these are the files we used to generate the tables and figures for the paper. 

The code was run with RStudio (v2023.06.2) in R software (v4.3.2). 

#### Figures:

The folder "figures" contain the figures that were taken directly from R. In some of them, we made some minor changes in Powerpoint.

##### Acknowledgments:
We thank Gary Packard for triggering this interesting debate on the topic of male dimorphism and static allometry, and also for revising an early version of the manuscript. We also thank an anonymous reviewer for comments and Matthew Symonds and Emma Sherratt for considering this submission for Evolutionary Ecology.

---

## Citation

If you use anything in this repository, please cite the paper:

> Buzatto, B.A., Machado, G. & Palaoro, A.V. (2024) Sigmoid allometries generate male dimorphism in secondary sexual traits: a comment on Packard (2023). *Evolutionary Ecology* 38(4): 537–548. [https://doi.org/10.1007/s10682-024-10303-6](https://doi.org/10.1007/s10682-024-10303-6)

If you reuse the code or the archived files directly, please also cite the archive:

> Buzatto, B.A., Machado, G. & Palaoro, A.V. (2024) *Data and code for: Sigmoid allometries generate male dimorphism in secondary sexual traits — a comment on Packard (2023)* [Data set]. Zenodo. https://doi.org/ZENODO_DOI_HERE

`CITATION.cff` in this repository holds both citations in machine-readable form — GitHub's
**Cite this repository** button (top right of the repository page) will generate APA or
BibTeX from it for you.

## License

This repository is released under two licenses, because code and data are different things.

| Content | License | File |
| --- | --- | --- |
| Analysis code — everything in `sim/code/` | [MIT](https://opensource.org/licenses/MIT) | [`LICENSE`](LICENSE) |
| Data, figures and stored model output — everything in `sim/data/`, `sim/simul_data/`, and `sim/figures/` | [CC0 1.0 Universal](https://creativecommons.org/publicdomain/zero/1.0/) | [`LICENSE-DATA`](LICENSE-DATA) |

In short: do whatever you like with the data, no permission needed and no attribution
legally required; reuse the code freely as long as you keep the copyright notice. Academic
norms still apply — if the data or code are useful to you, cite the paper.

## Archiving on Zenodo

This repository is not linked to Zenodo yet. To create a permanent, citable archive:

1. Go to [zenodo.org/account/settings/github](https://zenodo.org/account/settings/github),
   log in with GitHub, and flip the switch **on** for `alexandrepalaoro/simulations-bimodality`.
   Zenodo only sees public repositories, and it only archives releases created *after*
   the switch is on.
2. Back on GitHub, go to **Releases → Draft a new release**, create a tag (`v1.0.0`),
   give it a title, and publish it. Zenodo picks the release up within a minute or two
   and mints two DOIs: one for that specific version, and one **concept DOI** that always
   points at the newest version.
3. Copy the **concept DOI**. Paste it over `ZENODO_DOI_HERE` in this README, and in
   `CITATION.cff` uncomment the `doi:` line and put it there. Delete the notice at the top
   of this file. Always use the concept DOI in papers — it never goes stale.

`.zenodo.json` in this repository tells Zenodo the title, authors, ORCIDs, keywords, and
the link to the published article, so you do not have to retype any of it in the Zenodo form.

## Reproducibility

- Everything needed to reproduce the published analyses is in this repository.
- The archived release on Zenodo is the version of record. GitHub history may move on; the
  DOI will not.
- Package versions used are listed above. If a result does not reproduce, check those first.
