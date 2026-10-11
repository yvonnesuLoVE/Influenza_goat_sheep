# Cross-sectional Study of Highly Pathogenic Avian Influenza H5 Virus Exposure Among Goats and Sheep in Pakistan

## Overview
This repository contains the code and data used for the study:

**"Cross-sectional study of highly pathogenic avian influenza H5 virus exposure among goats and sheep in Pakistan"**

This study investigates influenza virus exposure in goats and sheep in Punjab Province, Pakistan, using a cross-sectional serological survey. Goat and sheep production is dense in regions where avian influenza viruses (AIV) are endemic, such as South Asia, which raises zoonotic concerns. While AIV circulation and persistence have been well described in domestic poultry, the exposure of small ruminants remains poorly understood.

To address this gap, we analyzed serum samples collected from goats and sheep between May and October 2023 across four districts of Punjab (Gujranwala, Kasur, Lahore and Sheikhupura). Samples were tested with a multiplex microsphere immunoassay targeting a broad panel of influenza hemagglutinin (HA) and neuraminidase (NA) antigens from avian, human and swine viruses. ELISA and hemagglutination inhibition assays were used to complement the multiplex data for selected subtypes.


## Key Findings

- A total of 781 goat and sheep serum samples (452 goats and 329 sheep) were analyzed from four districts of Punjab, Pakistan.
- High proportions of animals were MFI-positive for avian influenza HA antigens, including highly pathogenic avian influenza (HPAI) H5 clades 2.3.2.1c and 2.3.4.4b (23.9 to 34.0%), H7 (13.9 to 37.1%) and H9 (17.0 to 34.7%).
- Reactivity to human and swine H1N1/pdm09 and H3N2 HA antigens was absent or minimal.
- Multiplex immunoassay and ELISA measurements were positively correlated for H5 HA and N1 NA, indicating consistent antibody detection across assays.
- Goats and sheep are previously unrecognized hosts for viruses with zoonotic potential.

---

## Data Sources
- Goat and sheep serum samples collected from farms in Gujranwala, Kasur, Lahore and Sheikhupura districts, Punjab Province, Pakistan, between May and October 2023.
- Multiplex microsphere immunoassay (Luminex MAGPIX) targeting nine HA and three NA antigens of avian (H5, H7 and H9), human (H1N1/pdm09, seasonal H3N2 and influenza B/Victoria) and swine influenza viruses.
- ELISA and hemagglutination inhibition assays for selected subtypes.
- Analysis aimed at characterizing subtype-specific influenza exposure and cross-reactive antibody responses in small ruminants.

---

## Repository Contents
- `Goat_Sheep_AIV_serology_analysis.R`: complete analysis script. It covers data cleaning, exclusion of sera with MFI above 8,000 to more than one antigen, expectation-maximization (EM) classification of MFI values, seroprevalence by district, co-reactivity between antigens, heatmaps, statistical comparisons between goats and sheep, agreement between multiplex MFI and ELISA, and dose-response analysis.

---

## Implications
The findings identify goats and sheep as hosts with evidence of exposure to avian influenza viruses, including HPAI H5, and reveal gaps in our understanding of influenza virus transmission and adaptation. These results highlight the need to include small ruminants in influenza surveillance to enhance pandemic preparedness.

## Citation
If you use this repository in your research, please cite:

**[Your Full Citation Here]**
