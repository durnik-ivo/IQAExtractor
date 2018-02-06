# IQAExtractor

AWK scripts for extracting binging energies and delocalization indices from IQA analysis by program AIMALL. 

### Package consist of two scripts:
- [IQAEne.awk](./IQAEne.awk) - bindin energies
- [IQADI.awk](./IQDI.awk) - delocalization indices

### Usage:
- each script comes with **USER EDIT SECTION**, which needs to be set
- in terminal, scripts can be run as:
```
awk -f IQEne.awk
```
- to save text output into a file:
```
awk -f IQEne.awk > file.txt
```
- scripts are pre-set for examplary \*.sum files included in the repository
- switching control option on (`control = 1`) produces \*.xyz files to check user fragment definition 
(Use any molecular editor, e.g. [Avogadro](https://avogadro.cc/), [VMD](http://www.ks.uiuc.edu/Research/vmd/))

### Nomenclature:
Fragment(Geometry,Vicinity)
- Fragment:

   A or B

- Geometry:

   Complex or Opt

- Vicinity:

   A or B or Free

Intra = intraatomic contributions (within atom)  
Inter = interatomic contributions (between atoms)  
Total = Intra + Inter  
