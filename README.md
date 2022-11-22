# flair
FLAIR (Full-Length Alternative Isoform analysis of RNA) for the correction, isoform definition, and alternative splicing analysis of noisy reads. FLAIR has primarily been used for nanopore cDNA, native RNA, and PacBio sequencing reads. 

The complete Flair manual is available via [readthedocs](https://flair.readthedocs.io/en/latest/)

## Cite FLAIR <a name="cite"></a>
If you use or discuss FLAIR, please cite the following [paper](https://www.nature.com/articles/s41467-020-15171-6):
>Tang, A.D., Soulette, C.M., van Baren, M.J. et al. Full-length transcript characterization of SF3B1 mutation in chronic lymphocytic leukemia reveals downregulation of retained introns. Nat Commun 11, 1438 (2020). https://doi.org/10.1038/s41467-020-15171-6

## Installation and usage of this fork
The script below creates an environment to execute this fork. 
```
# install flair
git clone https://github.com/ChangLabSNU/flair
cd flair
conda env create -f environment.yaml
conda activate flair_3p
pip install .

# check if flair is running properly
flair
```
In this fork, an argument `--wiggleWindow_3p`(`-3`) is added to define transcription end site correction window flank size.


