######################################################################
# Running OSSAN-Euler2D for a shock diffraction problem.
#
# Follow the steps to run the shock diffraction test case.
# (Or just run this script by typing "source readme.txt".)
#
# NOTE: The following is for gfortran. If you use other compiler,
#       replace "gfortran" by your compiler (ifort, g95, etc).
#
#
# Katate Masatsuka, January 2012 (http://www.cfdbooks.com)
#####################################################################

#####################################################################
# 1. Compile the grid generation code.
#####################################################################

gfortran -O2 -o twod_rectangular_grid twod_rectangular_grid_v0.f90

#####################################################################
# 2. Run and generate grids (quads and triangles).
#    Dimensions are defined in the program: 401x401 grid.
#    It will generate the following files:
#     quad_grid_tecplot.dat - for viewing the quad grid
#     tria_grid_tecplot.dat - for viewing the tria grid
#     quad.grid             - for running OSSAN Euler code
#     tria.grid             - for running OSSAN Euler code
#     project.bcmap         - for running OSSAN Euler code
#####################################################################

./twod_rectangular_grid

#####################################################################
# 3. Compile OSSAN-Euler2D
#####################################################################

gfortran -O2 -c data_package_v0.f90
gfortran -O2 -c euler_solver_v0.f90
gfortran -O2 -c main_v0.f90
gfortran -O2 -o ossan_euler2d data_package_v0.o euler_solver_v0.o main_v0.o

#####################################################################
# 4. Rename the quad grid file
#    (NOTE: OSSAN-Euler2D reads project.grid and project.bcmap)
#####################################################################

cp quad.grid   project.grid

#####################################################################
# 5. Run OSSAN-Euler2D to compute the solution on the quad grid 
#    (NOTE: All input parametetrs are specified inside the program.)
#    It will stop at t = 0.18 (~ 1650 steps).
#####################################################################

./ossan_euler2d

#####################################################################
# 6. Save the output file
#    (NOTE: The file 'project_tecplot.dat' will be overwritten by a
#           subsequent run.)
#####################################################################

cp project_tecplot.dat project_quad_result_tecplot.dat

#####################################################################
# Move on to Step 7 if you want to run it also on the triangular grid.
#
# 7. Rename the triangular grid file
#    (NOTE: bcmap file is the same for quad and triangular grids)
#####################################################################

cp tria.grid project.grid

#####################################################################
# 8. Run OSSAN-Euler2D to compute the solution on the triangular grid
#    It will stop at t = 0.18 (>1000 steps).
#####################################################################

./ossan_euler2d

#####################################################################
# 9. Save the output file if you wish.
#####################################################################

cp project_tecplot.dat project_tria_result_tecplot.dat

