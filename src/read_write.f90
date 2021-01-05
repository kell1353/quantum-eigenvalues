!-----------------------------------------------------------------------
!Module: read_write
!-----------------------------------------------------------------------
!! Rodrigo Navarro Perez
!!
!! Contains subroutines and functions related to reading input from the
!! user and  writing output into a text file
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! read_input
!! read_advanced_input
!! print_energies
!! write_probability_density
!! write_woods_saxon_energies
!!----------------------------------------------------------------------
!! Included functions:
!!
!! read_real
!! read_integer
!-----------------------------------------------------------------------
module read_write
use types
use qm_solver, only: solve_woods_saxon
implicit none

private
public :: read_input, write_probability_density, print_energies, write_woods_saxon_energies

contains

!-----------------------------------------------------------------------
!! Subroutine: read_input
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! Displays a message describing what the program does and the expected
!! input. After that it uses the `read_real` and `read_integer`
!! functions to assign values to the different parameters.
!!----------------------------------------------------------------------
!! Output:
!!
!! n_points     integer     number of grid points the discretized wave function
!! length       real        length of the box
!! radius       real        radius of the Woods-Saxon potential
!-----------------------------------------------------------------------
subroutine read_input(size_box, num_lattice_points, WS_radius)
    implicit none
    real(dp), intent(out) :: size_box, WS_radius
    integer, intent(out) ::  num_lattice_points


    print *, 'This program takes an individual real and integer inputs and calculates the energies and'
    print *, 'wave functions using numerical sampling techniques to solve the Schrodinger Equation'
    print *, 'at different potentials.'

    size_box = read_real('Size of the Box, L')
    num_lattice_points = read_integer('number of points on your lattice, N')
    WS_radius = read_real('radius of Woods-Saxon potential, R')

end subroutine read_input

!-----------------------------------------------------------------------
!! Function: read_real
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! When reading real input from a user, checks have to be made to make 
!! sure that the user provided the correct type of input. 
!! 
!! We enclose the input reading inside an infinite loop that can only
!! be exited when a correct input is given.
!! 
!! The first check is to make sure that the user input is positive and non-zero.
!!
!! The second check is to make sure that the string is not empty 
!! (i.e. the user simply pressed the enter key)
!! 
!! The third check is made by using the 'read' statement to convert
!! the string into a number, if that is not possible iostat gives an
!! error code different from zero.
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! name     character   A string with a brief description of the value being asked for
!!----------------------------------------------------------------------
!! Variables:
!!
!! string   character   The intial user input converted to a string
!! ierror   integer     The integer value that represents if there is an error
!!----------------------------------------------------------------------
!! Output:
!!
!! x        real        A positive non negative number given by the user
!! string   character   The intial user input converted to a string
!! ierror   integer     The integer value that represents if there is an error
!-----------------------------------------------------------------------
real(dp) function read_real(name) result(x)
    implicit none
    character(len=*), intent(in) :: name
    character(len=120) :: string
    integer :: ierror

    print *, 'Provide a nonzero positive value for the '//trim(name)//':'

    do
        read(*,'(a)',iostat=ierror) string
        if(string /= '') then
            read(string, *, iostat=ierror) x
            if (ierror == 0 .and. x > 0 ) exit
            print *, "'"//trim(string)//"'"//' is not a valid number, please provide a positive real number'
        else
            print *, 'that was an empty input, please provide a positive real number'
        endif
    enddo
end function read_real

!-----------------------------------------------------------------------
!! Function: read_integer
!-----------------------------------------------------------------------
!! By: Austin Keller
!!
!! When reading integer input from a user, checks have to be made to make 
!! sure that the user provided the correct type of input. 
!! 
!! We enclose the input reading inside an infinite loop that can only
!! be exited when a correct input is given.
!! 
!! The first check is to make sure that the user input is positive and non-zero.
!!
!! The second check is to make sure that the string is not empty 
!! (i.e. the user simply pressed the enter key)
!! 
!! The third check is made by using the 'read' statement to convert
!! the string into a number, if that is not possible iostat gives an
!! error code different from zero.
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! name     character   A string with a brief description of the value being asked for
!!----------------------------------------------------------------------
!! Variables:
!!
!! string   character   The intial user input converted to a string
!! ierror   integer     The integer value that represents if there is an error
!!----------------------------------------------------------------------
!! Output:
!!
!! x        integer     A positive non negative number given by the user
!-----------------------------------------------------------------------
integer function read_integer(name) result(x)
    implicit none
    character(len=*), intent(in) :: name
    character(len=120) :: string
    integer :: ierror

    print *, 'Provide a nonzero positive integer for the '//trim(name)//':'

    do
        read(*,'(a)',iostat=ierror) string
        if(string /= '') then
            read(string, *, iostat=ierror) x
            if (ierror == 0 .and. x > 0 ) exit
            print *, "'"//trim(string)//"'"//' is not a valid number, please provide a positive real number'
        else
            print *, 'that was an empty input, please provide a positive real number'
        endif
    enddo
end function read_integer


!-----------------------------------------------------------------------
!! Subroutine: write_probability_density
!-----------------------------------------------------------------------
!! Austin Keller
!!
!! This function has inputs of the file name, x sample points, and wave
!! functions. Writes the sample points, grounded state, 1st excited state
!! and 2nd excited state probability densities to the specified file name.
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! file_name        character   file for the values to be printed to 
!! x_vector         real        array of samples points within the box
!! wave_functions   real        array containing the numerical results for the 
!!                              wave functions
!!----------------------------------------------------------------------
subroutine write_probability_density(file_name, x_vector, wave_functions)
    implicit none
    character(len=*), intent(in) :: file_name
    real(dp), intent(in) :: x_vector(:), wave_functions(:,:)
    integer :: unit, i, n

    n = size(x_vector)

    open(newunit=unit,file=trim(file_name))
    write(unit,'(4a28)') 'x', 'ground state', '1st excited', '2nd excited'
    
    do i=1,n
        write(unit, *) x_vector(i), wave_functions(i,1)**2, wave_functions(i,2)**2, wave_functions(i,3)**2
    enddo

    close(unit)

    print *, 'the states of the wave function were written in ', file_name
    
end subroutine write_probability_density


!-----------------------------------------------------------------------
!! Subroutine: print_energies
!-----------------------------------------------------------------------
!! Austin Keller
!!
!! This subroutine prints the numerical and anlytic energy 
!! solutions. 
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! name             character   Title of the energies being printed
!-----------------------------------------------------------------------
!! Output:
!!
!! numerical        array       Array cointing the numerical energies
!! analytic         array       Array containing the analytical energies
!!
!-----------------------------------------------------------------------
subroutine print_energies(name, numerical, analytic)
    implicit none
    character(len=*), intent(in) :: name
    real(dp), intent(out) :: numerical(:), analytic(:)
    integer :: i, n

    if (size(numerical) /= size(analytic)) then
        print*, "arrays size don't match in print_energies"
        stop
    endif

    n = size(numerical)

    print *, '' !Just to space out prints in the terminal
    print*, 'Comparing numerical and anlytic solutions in'
    print*, trim(name)
    
    print'(a9,2a15)', 'number', 'numerical', 'analytic'

    do i=1,n
        print*, i, numerical(i), analytic(i)
    enddo
    print *, ''
end subroutine print_energies


!-----------------------------------------------------------------------
!! Subroutine: write_woods_saxon_energies
!-----------------------------------------------------------------------
!! Austin Keller
!!
!! This subroutine takes in the inputs specified below, allocates the arrays
!! for the eigenvalues and eigenvectors. 
!!
!! Initiates a loop to solve for the Woods-Saxon wave functions for radii 
!! between the minimum and maximum inputs provided and writes the results
!! into a file.
!!
!!----------------------------------------------------------------------
!! Input:
!!
!! file_name             character   file for the values to be printed to
!! length                real        length of our sample box
!! x_vector              real        array of samples points within the box
!! r_min                 real        the minimum value for the radius of the WS potential
!! r_max                 real        the maximum value for the radius of the WS potential
!!----------------------------------------------------------------------
subroutine write_woods_saxon_energies(file_name, length, x_vector, r_min, r_max)
    implicit none
    character(len=*), intent(in) :: file_name
    real(dp), intent(in) :: length, x_vector(:), r_min, r_max

    real(dp) :: r, r_step
    integer :: unit, n
    real(dp), allocatable :: energies(:), wave_functions(:,:)

    n = size(x_vector)

    ! Allocate arrays
    if(allocated(energies)) deallocate(energies)
    allocate(energies(1:n))

    if(allocated(wave_functions)) deallocate(wave_functions)
    allocate(wave_functions(1:n, 1:n))

    ! Write to file
    open(newunit=unit,file=trim(file_name))
    write(unit,'(4a28)') 'radius', 'ground state', '1st excited', '2nd excited'
    
    r = r_min
    r_step = .5_dp
    do 
        if(r > r_max) exit
        call solve_woods_saxon(length, r, x_vector, energies, wave_functions)
        write(unit, *) r, energies(1), energies(2), energies(3)
        r = r + r_step
    enddo

    close(unit)

    print *, 'the three lowest energies per given radius of the Woods-Saxon potential written here ', file_name
end subroutine write_woods_saxon_energies


end module read_write
