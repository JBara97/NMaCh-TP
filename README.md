# NMaCh-TP
Networked Markov Chain model for Tipping Point risk assessment. Each tipping point is modelled as a two-state Markov chain in which its transition rates are affected by the state of its in-neighbours. In particular a destabilising in-neighbour (which is in the degraded state) will increase the collapse rate, while a stabilising in-neighbour (again, in its degraded state) will decrease the collapse rate.

**Methods paper**: 
Bara, J., Wunderling, N. & Barfuss, W. A risk assessment framework for interacting tipping elements. Earth Syst. Dyn. **17**, 333–352, DOI: [10.5194/esd-17-333-2026](10.5194/esd-17-333-2026) (2026)

## External data sources
Individual tipping element data (namely, critical threshold estimates)
  - Armstrong McKay et al. 2022 [1]
  - Global Tipping Points Report 2023 [2]
  - Global Tipping Points Report 2025 [3]*
Interaction data
  - Kriegler et al. 2009 [4] (probability factors)
  - Wunderling et al. 2021 [5] (maximum link strength)
  - Wunderling et al. 2024 [6]* (qualitative links)
Global warming temperature scenarios
  - Lee et al. 2021 [7] (extended SSPs)
  - Lamboll et al. 2025 [8] (PROVIDE scenarios)
*Most recent / comprehensive data used.

## References
[1] - Armstrong McKay, D. I. et al. Exceeding 1.5°C global warming could trigger multiple climate tipping points. Science 377, DOI: [10.1126/science.abn7950](10.1126/science.abn7950) (2022)
[2] - Lenton, T. M. et al. (eds.) The Global Tipping Points Report 2023 (University of Exeter, Exeter, UK, 2023).
[3] - Lenton, T. M. et al. (eds.) The Global Tipping Points Report 2025 (University of Exeter, Exeter, UK, 2025).
[4] - Kriegler, E. et al. Imprecise probability assessment of tipping points in the climate system. Proc. Natl. Acad. Sci. 106, 5041–5046, DOI: [10.1073/pnas.0809117106](10.1073/pnas.0809117106) (2009)
[5] - Wunderling, N. et al. Interacting tipping elements increase risk of climate domino effects under global warming. Earth Syst. Dyn. 12, 601–619, DOI: [10.5194/esd-12-601-2021](10.5194/esd-12-601-2021) (2021).
[6] - Wunderling, N. et al. Climate tipping point interactions and cascades: a review. Earth Syst. Dyn. 15, 41–74, DOI: [10.5194/esd-15-41-2024](10.5194/esd-15-41-2024) (2024).
[7] - Lee, J.-Y. et al. Future global climate: Scenario-based projections and near-term information. In Masson-Delmotte, V. et al. (eds.) Climate Change 2021 – The Physical Science Basis: Working Group I Contribution to the Sixth Assessment Report of the Intergovernmental Panel on Climate Change, 553–672 (Cambridge University Press, 2021)
[8] - Lamboll, R., Rogelj, J. & Schleussner, C.-F. Scenario emissions and temperature data for PROVIDE project (v1.3.0) [data set], DOI: [10.5281/zenodo.11108741](10.5281/zenodo.11108741) (2025).
