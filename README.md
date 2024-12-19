# Supplementary data for the paper "Protons accumulate at the graphene-water interface"

This repository contains supplementary data supporting the findings of the paper:

"Protons accumulate at the graphene-water interface" (Xavier R. Advincula, Kara D. Fong, Angelos Michaelides, Christoph Schran)
[under revision]

## License
The content of this repository is licensed under the CC-BY-SA-4.0 license. See the file `LICENSE` for details.

## Contents
* `model`: The MACE MLP model developed and used in this work.
* `dft-inputs`: Example input files of the electronic structure settings used to produce (i) the reference DFT data, (ii) the electron density difference & Bader charges analsys.
* `aimd-inputs`: Example input files of (i) the settings of short AIMD simulations to generate training data for the development of the MACE MLP, (ii) additional AIMD simulations to generate reference data for the validation of the MLP for both bulk and nanoconfined environments. 
* `mlp-based-md`: Example input files of extensive unbiased MD simulations using the developed MLP, which form the main focus of this work.
