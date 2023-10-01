# PlanetFormation SmallStars

Code and files published in Stefansson et al. 2023. Code developed in Miguel et al. 2020. 

# Installation and Running
To compile the code you can use any standad Fortran compiler. You can compile and run the code with the following commands:
```
gfortran -o formation Formation_Small_Stars-2020.f90
./formation
```

To obtain the results published in Stefansson et al. 2023, change the mass of the star and the disk mass according to the conditions described in the paper for each of the cases explored. 
The parameters for the distribution of mass of the disk can also change. 

All description can be found in detail on the Supplementary material of Stefansson et al. 2023.

# Notes

The `notebooks/` directory has an ipython notebook (`Plotting Figure 3 in Stefansson et al. 2023.ipynb`) to regenerate Figure 3 from Stefansson et al. 2023.
The corresponding datafiles are in the `data/` directory.

# Citation
If you use this code, please cite the following two papers:

- Miguel et al. 2020, MNRAS 491, 1998–2009, https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.1998M/abstract
- Stefansson et al. 2023, Science, https://ui.adsabs.harvard.edu/abs/2023arXiv230313321S/abstract
