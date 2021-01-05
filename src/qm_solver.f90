!-----------------------------------------------------------------------
!Module: qm_solver
!-----------------------------------------------------------------------
!! Austin Keller
!!
!! The purpose of this module is to solve all of our analytical and numerical
!! Schrodinger equations in one place. 
!! 
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! sample_box
!! solve_infinite_well
!! analytic_infinite_well
!! solve_harmonic_oscillator
!! analytic_harmonic_oscillator
!! solve_woods_saxon
!!----------------------------------------------------------------------

module qm_solver
use types
use hamiltonian, only : construct_hamiltonian, harmonic_potential_energy, WS_potential_energy, h_bar, mass
use eigen_solver, only : solve_eigenproblem
implicit none


private
public :: solve_infinite_well, sample_box, analytic_infinite_well, solve_harmonic_oscillator, &
    analytic_harmonic_oscillator, solve_woods_saxon

contains

!-----------------------------------------------------------------------
!! Subroutine: sample_box
!-----------------------------------------------------------------------
!! Austin Keller
!!
!! This subroutine specifies sample points x within our sample box of 
!! length L. 
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! L 				real		Length of the sample box
!! N 				integer 	Number of sample points
!-----------------------------------------------------------------------
!! Output:
!!
!! x 				array		Array containing the sample points in the box
!-----------------------------------------------------------------------
subroutine sample_box(L, N, x)
    implicit none
    real(dp), intent(in) :: L
    integer, intent(in) ::  N
    real(dp), intent(out) :: x(:)

    real(dp) :: dx
    integer :: i

    dx = (2*L)/(N-1)

    x = 0._dp
    do i=1,N
        x(i) = -L + dx*(i-1)
	enddo

end subroutine sample_box

!-----------------------------------------------------------------------
!! Subroutine: solve_infinite_well
!-----------------------------------------------------------------------
!! Austin Keller
!!
!! This subroutine takes in the length of the box. It constructs the 
!! hamiltonian matrix and then numerically solves the eigen 
!! problem for the energies and wave functions for the particle in a box 
!! scenario.
!!
!! Then normalizes the wave function results.
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! length 				real		Length of the sample box
!-----------------------------------------------------------------------
!! Output:
!!
!! energies 			real		Array containing the numerical result for the 
!! 									energies of the harmonic oscillator
!! wave_functions 		real		Array containing the numerical result for the 
!! 									wave_functions of the harmonic oscillator
!-----------------------------------------------------------------------
subroutine solve_infinite_well(length, energies, wave_functions)
    implicit none
    real(dp), intent(in) :: length
    real(dp), intent(out) :: energies(:), wave_functions(:,:)

    real(dp) :: delta
    integer :: n
    real(dp), allocatable :: potential_diagonal(:), hamiltonian_diagonal(:), hamiltonian_off_diagonal(:)

    n = size(energies)

    ! Allocate Arrays
    if(allocated(potential_diagonal)) deallocate(potential_diagonal)
    allocate(potential_diagonal(1:n))

    if(allocated(hamiltonian_diagonal)) deallocate(hamiltonian_diagonal)
    allocate(hamiltonian_diagonal(1:n))

    if(allocated(hamiltonian_off_diagonal)) deallocate(hamiltonian_off_diagonal)
    allocate(hamiltonian_off_diagonal(1:n-1))

    delta = (2*length)/(n-1)
    potential_diagonal = 0._dp

    call construct_hamiltonian(delta, potential_diagonal, hamiltonian_diagonal, hamiltonian_off_diagonal)
    call solve_eigenproblem(hamiltonian_diagonal, hamiltonian_off_diagonal, energies, wave_functions)
    

    call normalize(delta, wave_functions)

end subroutine solve_infinite_well

!-----------------------------------------------------------------------
!! Subroutine: analytic_infinite_well
!-----------------------------------------------------------------------
!! Austin Keller
!!
!! This subroutine takes in the length of the box and the amount of energies.
!! It calculates analytically the n lowest energies for the particle in 
!! a box scenario.
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! L 						real	Length of the sample box
!! num_energies 			real	Parameter for the number of energies to calculate
!-----------------------------------------------------------------------
!! Output:
!!
!! analytic_energies		array	Array containing the analytic result of the n 
!! 									lowest energies for the infinite well
!-----------------------------------------------------------------------
subroutine analytic_infinite_well(L, num_energies, analytic_energies)
    implicit none
    real(dp), intent(in) :: L
    integer, intent(in) :: num_energies
    real(dp), intent(out) :: analytic_energies(:)

    integer :: i, n

    n = num_energies

    analytic_energies = 0._dp
    do i=1,n
        analytic_energies(i) = (i**2)*(((h_bar**2)*(pi**2))/(8*mass*L**2))
	enddo

    !print *, analytic_energies

end subroutine analytic_infinite_well

!-----------------------------------------------------------------------
!! Subroutine: solve_harmonic_oscillator
!-----------------------------------------------------------------------
!! Austin Keller
!!
!! This subroutine takes in the length of the box. It constructs the 
!! hamiltonian matrix and then numerically solves the eigen 
!! problem for the energies and wave functions for the harmonic oscillating
!! potential.
!!
!! Then normalizes the wave function results.
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! length 				real		Length of the sample box
!! x_vector 			real		Array of samples points within the box
!-----------------------------------------------------------------------
!! Output:
!!
!! energies 			real		Array containing the numerical result for the 
!! 									energies of the harmonic oscillator
!! wave_functions 		real		Array containing the numerical result for the 
!! 									wave_functions of the harmonic oscillator
!-----------------------------------------------------------------------
subroutine solve_harmonic_oscillator(length, x_vector, energies, wave_functions)
    implicit none
    real(dp), intent(in) :: length, x_vector(:)
    real(dp), intent(out) :: energies(:), wave_functions(:,:)

    real(dp) :: delta
    integer :: n
    real(dp), allocatable :: harmonic_potential_diagonal(:), hamiltonian_diagonal(:), hamiltonian_off_diagonal(:)

    n = size(energies)

    ! Allocate Arrays
    if(allocated(harmonic_potential_diagonal)) deallocate(harmonic_potential_diagonal)
    allocate(harmonic_potential_diagonal(1:n))

    if(allocated(hamiltonian_diagonal)) deallocate(hamiltonian_diagonal)
    allocate(hamiltonian_diagonal(1:n))

    if(allocated(hamiltonian_off_diagonal)) deallocate(hamiltonian_off_diagonal)
    allocate(hamiltonian_off_diagonal(1:n-1))

    delta = (2*length)/(n-1)

    call harmonic_potential_energy(x_vector, harmonic_potential_diagonal)
    call construct_hamiltonian(delta, harmonic_potential_diagonal, hamiltonian_diagonal, hamiltonian_off_diagonal)
    call solve_eigenproblem(hamiltonian_diagonal, hamiltonian_off_diagonal, energies, wave_functions)

    call normalize(delta, wave_functions)

end subroutine solve_harmonic_oscillator

!-----------------------------------------------------------------------
!! Subroutine: analytic_harmonic_oscillator
!-----------------------------------------------------------------------
!! Austin Keller
!!
!! This subroutine takes in the length of the box and the amount of energies.
!! It calculates analytically the n lowest energies for the harmonic 
!! oscillator potential.
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! num_energies 			array	Parameter for the number of energies to calculate
!-----------------------------------------------------------------------
!! Output:
!!
!! analytic_energies		array	Array containing the analytic result of the n 
!! 									lowest energies for the harmonic oscillator
!-----------------------------------------------------------------------
subroutine analytic_harmonic_oscillator(num_energies, analytic_energies)
    implicit none
    integer, intent(in) :: num_energies
    real(dp), intent(out) :: analytic_energies(:)

    integer :: i, n

    n = num_energies

    analytic_energies = 0._dp
    do i=1,n
        analytic_energies(i) = (i-0.5)*((h_bar**2)/mass)
	enddo

end subroutine analytic_harmonic_oscillator

!-----------------------------------------------------------------------
!! Subroutine: solve_woods_saxon
!-----------------------------------------------------------------------
!! Austin Keller
!!
!! This subroutine takes in the length of the box. It constructs the 
!! hamiltonian matrix and then numerically solves the eigen 
!! problem for the energies and wave functions for the Woods-Saxon 
!! potential.
!!
!! Then normalizes the wave function results.
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! length 				real		Length of the sample box
!! x_vector 			real		Array of samples points within the box
!-----------------------------------------------------------------------
!! Output:
!!
!! energies 			real		Array containing the numerical result for the 
!! 									energies of the Woods-Saxon
!! wave_functions 		real		Array containing the numerical result for the 
!! 									wave_functions of the Woods-Saxon
!-----------------------------------------------------------------------
subroutine solve_woods_saxon(length, radius, x_vector, energies, wave_functions)
    implicit none
    real(dp), intent(in) :: length, radius, x_vector(:)
    real(dp), intent(out) :: energies(:), wave_functions(:,:)

    real(dp) :: delta
    integer :: n
    real(dp), allocatable :: WS_potential_diagonal(:), hamiltonian_diagonal(:), hamiltonian_off_diagonal(:) 

    n = size(energies)

    ! Allocate Arrays
    if(allocated(WS_potential_diagonal)) deallocate(WS_potential_diagonal)
    allocate(WS_potential_diagonal(1:n))

    if(allocated(hamiltonian_diagonal)) deallocate(hamiltonian_diagonal)
    allocate(hamiltonian_diagonal(1:n))

    if(allocated(hamiltonian_off_diagonal)) deallocate(hamiltonian_off_diagonal)
    allocate(hamiltonian_off_diagonal(1:n-1))

    delta = (2*length)/(n-1)

    call WS_potential_energy(radius, x_vector, WS_potential_diagonal)
    call construct_hamiltonian(delta, WS_potential_diagonal, hamiltonian_diagonal, hamiltonian_off_diagonal)
    call solve_eigenproblem(hamiltonian_diagonal, hamiltonian_off_diagonal, energies, wave_functions)

    call normalize(delta, wave_functions)

end subroutine solve_woods_saxon



!-----------------------------------------------------------------------
!! Subroutine: normalize
!-----------------------------------------------------------------------
!! Austin Keller
!!
!! This subroutine takes in the delta and wave functions and calculates
!! the normalized constant for the wave functions. Then applies the constant
!! to the values of the wave function to normalize the results.
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! dx 					real 	the delta value for the numerical integration
!! wave_functions 		real 	array containing the numerical result for the 
!! 								wave_functions 
!-----------------------------------------------------------------------
!! Output:
!!
!! wave_functions 		real	array containing the normalization of the 
!! 								wave functions
!-----------------------------------------------------------------------
subroutine normalize(dx, wave_functions)
    implicit none
    real(dp), intent(inout) :: wave_functions(:,:)
    real(dp), intent(in) :: dx
    integer :: i, n(1:2) 
    real(dp), allocatable :: A, sum_psi(:), vectors(:), normalized_constants(:)


    n = shape(wave_functions)
    allocate(normalized_constants(1:n(1)))
    allocate(sum_psi(1:n(1)))

    !print *, sum(dx*(wave_functions(:,1))**2)
    do i=1,n(1)
        sum_psi(i) = (sum(dx*(wave_functions(:,i))**2))

        normalized_constants(i) = sqrt(1/(sum_psi(i)))
        wave_functions(:,i) = wave_functions(:,i)*(normalized_constants(i))
    enddo

end subroutine normalize


end module qm_solver