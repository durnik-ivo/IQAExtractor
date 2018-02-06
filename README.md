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

### Nomenclature:
Fragment(Geometry,Vicinity)
- Fragment
    - A
    - B
- Geometry
    - Complex
    - Opt
- Vicinity
    - A
    - B
    - Free
