# **Initial configuration**
The initial geometry of Taylor bubble 
## **Sensitivity_analysis**
- [General case](https://github.com/GabrielGLK/Taylor-bubble/blob/main/straight-pipe/sensitive-analysis/taylor_2D.c): you can change the dimensionless numbers to realize the sensitivity analysis.
- [Output_taylor](https://github.com/GabrielGLK/Taylor-bubble/blob/main/straight-pipe/sensitive-analysis/output_taylor.h): some output fields.
- [Makefile](https://github.com/GabrielGLK/Taylor-bubble/blob/main/straight-pipe/sensitive-analysis/Makefile): openmpi method, you can also change the numbers of processors in OMP_NUM_THREADS=10. For MPI, the detailed information can be found in Basilisk (add the link here, can't open the website now).
- [run.sh](https://github.com/GabrielGLK/Taylor-bubble/blob/main/straight-pipe/sensitive-analysis/run.sh): bash run.sh -> run the case
