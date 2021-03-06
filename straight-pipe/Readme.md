# **[Initial configuration](https://github.com/GabrielGLK/Taylor-bubble/blob/main/straight-pipe/figure/straight.pdf)**
## **Sensitivity analysis**
- [General case](https://github.com/GabrielGLK/Taylor-bubble/blob/main/straight-pipe/sensitive-analysis/taylor_2D.c): you can change the dimensionless numbers to realize the sensitivity analysis.
- [Output_taylor](https://github.com/GabrielGLK/Taylor-bubble/blob/main/straight-pipe/sensitive-analysis/output_taylor.h): some output fields.
- [Makefile](https://github.com/GabrielGLK/Taylor-bubble/blob/main/straight-pipe/sensitive-analysis/Makefile): openmpi method, you can also change the numbers of processors in OMP_NUM_THREADS=10. For MPI, the detailed information can be found in Basilisk (add the link here, can't open the website now).
- [run.sh](https://github.com/GabrielGLK/Taylor-bubble/blob/main/straight-pipe/sensitive-analysis/run.sh): bash run.sh -> run the case
- [output-data](https://github.com/GabrielGLK/Taylor-bubble/tree/main/straight-pipe/sensitive-analysis/output_data): I give some data in my thesis, but you can calculate that by yourself. It is noted in each folder the 'plot' is used for gnuplot, just copy that in gnuplot terminal.
### some other cases
Of course, the sensitivity analysis does not only contains contents above, but also timestep analysis, film thickness, bubble length... It is easy to change [General case](https://github.com/GabrielGLK/Taylor-bubble/blob/main/straight-pipe/sensitive-analysis/taylor_2D.c). Here this is not the main part of this section.

## **Experiments && Simulations**
Repeated instructions I will not write here.
- [experiment.h](https://github.com/GabrielGLK/Taylor-bubble/blob/main/straight-pipe/experiments-simulations/experiment.h): this is the new experimental data we have measured. It includes all the cases we want to simulate to compare with Chengsi's experiments data.
- [reference.h](https://github.com/GabrielGLK/Taylor-bubble/blob/main/straight-pipe/experiments-simulations/reference.h): reference data from article of Araujo (2012).
- [particle.h](https://github.com/GabrielGLK/Taylor-bubble/blob/main/straight-pipe/experiments-simulations/particle.h): liquid particle method.
- [scatter.h](https://github.com/GabrielGLK/Taylor-bubble/blob/main/straight-pipe/experiments-simulations/scatter.h): used in view.h to show liquid trajectory.
- [reference_compare.py](https://github.com/GabrielGLK/Taylor-bubble/blob/main/straight-pipe/experiments-simulations/python/reference_compare.py): post-processing.