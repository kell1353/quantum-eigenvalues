!-----------------------------------------------------------------------------
!! Program: schrodinger_solution
!! Austin Keller
!!
!! This program provides numerical and analytical solutions for the energies
!! at a few different potentials. The infinite well, harmonic 
!! oscillator and the Woods-Saxon potential. It writes these solutions into
!! files for each potential for a few different states.
!! 
!! This program calculates numerically eigenenergies and eigenfunctions 
!! for the Woods Saxon potential and provides output for analysis.
!!
!-----------------------------------------------------------------------------
program schrodinger_solution 

use types
use read_write, only : read_input, write_probability_density, print_energies, write_woods_saxon_energies
use qm_solver, only: sample_box, solve_infinite_well, analytic_infinite_well, solve_harmonic_oscillator,&
    analytic_harmonic_oscillator, solve_woods_saxon
implicit none

integer :: n_points
real(dp) :: length, radius

integer, parameter :: n_energies = 3
real(dp) :: analytic_energies(1:n_energies)
real(dp), allocatable :: all_energies(:), wave_functions(:,:)
real(dp), allocatable :: sample_points(:), x_vector(:)
real(dp), parameter :: r_min = 2._dp, r_max = 10._dp

call read_input(length, n_points, radius)

!Allocate arrays
if(allocated(x_vector)) deallocate(x_vector)
    allocate(x_vector(1:n_points))
if(allocated(all_energies)) deallocate(all_energies)
    allocate(all_energies(1:n_points))
if(allocated(wave_functions)) deallocate(wave_functions)
    allocate(wave_functions(1:n_points, 1:n_points))

call sample_box(length, n_points, x_vector)
! Solving particle in a box
call solve_infinite_well(length, all_energies, wave_functions)
call analytic_infinite_well(length, n_energies, analytic_energies)
call print_energies('Infinite Well', all_energies(1:n_energies), analytic_energies)
call write_probability_density('infinite_well_wf.dat', x_vector, wave_functions)

! Solving harmonic oscillator
call solve_harmonic_oscillator(length, x_vector, all_energies, wave_functions)
call analytic_harmonic_oscillator(n_energies, analytic_energies)
call print_energies('Harmonic oscillator', all_energies(1:n_energies), analytic_energies)
call write_probability_density('harmonic_oscillator_wf.dat', x_vector, wave_functions)

! Solving Woods Saxon
call solve_woods_saxon(length, radius, x_vector, all_energies, wave_functions)
call write_probability_density('woods_saxon_wf.dat', x_vector, wave_functions)

! Woods Saxon Energies as a function of radius
call write_woods_saxon_energies('woods_saxon_energy.dat', length, x_vector, r_min, r_max)

end program schrodinger_solution