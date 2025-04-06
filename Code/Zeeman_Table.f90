program magnetic_interactions
    implicit none
  
    ! Constants
    real(16), parameter :: FINE_STRUCTURE_CONST = 1.0q0 / 137.0q0
    real(16), parameter :: H_GROUND_STATE_ENERGY = -0.5q0
    real(16), parameter :: PROTON_MASS_AU = 1836.15267q0
    real(16), parameter :: NUCLEAR_MAGNETON = 1.0q0 / (2.0q0 * PROTON_MASS_AU)
    real(16), parameter :: NUCLEAR_G_FACTOR = -4.25509q0
    real(16), parameter :: BOHR_MAGNETON = 1.0q0 / 2.0q0
    real(16), parameter :: ELECTRONIC_G_FACTOR = 2.0q0
    real(16), parameter :: NUCLEAR_CHARGE = 2.0q0
  
    ! Local variables
    real(16) :: B_scaled
    real(16) :: ratio_val, ordinary_val, new_val
    integer :: i
    real(16), dimension(4) :: fields
    character(len=50), dimension(4) :: Labels
  
    ! Initialize Labels and fields
    Labels = [ '$m_I=-\\frac{1}{2}, m_s=-\\frac{1}{2}$', &
               '$m_I=-\\frac{1}{2}, m_s=+\\frac{1}{2}$', &
               '$m_I=+\\frac{1}{2}, m_s=-\\frac{1}{2}$', &
               '$m_I=+\\frac{1}{2}, m_s=+\\frac{1}{2}$' ]
  
    fields = [1.45q0, 5.7q0, 11.7q0, 45.5q0]
  
    ! Loop over each field value
    do i = 1, 4
       B_scaled = fields(i) / 2.3505175q5
       print *, 'Field (T): ', fields(i)
  
       ! Print Ratios
       print *, 'RATIOS'
       print *, ratio(B_scaled, 0.5q0, 0.5q0)
       print *, ratio(B_scaled, 0.5q0, -0.5q0)
       print *, ratio(B_scaled, -0.5q0, 0.5q0)
       print *, ratio(B_scaled, -0.5q0, -0.5q0)
  
       ! Print Ordinary
       print *, 'ORDINARY'
       print *, ordinary(B_scaled, 0.5q0, 0.5q0)
       print *, ordinary(B_scaled, 0.5q0, -0.5q0)
       print *, ordinary(B_scaled, -0.5q0, 0.5q0)
       print *, ordinary(B_scaled, -0.5q0, -0.5q0)
  
       ! Print New
       print *, 'NEW'
       print *, new(B_scaled, 0.5q0, 0.5q0)
       print *, new(B_scaled, 0.5q0, -0.5q0)
       print *, new(B_scaled, -0.5q0, 0.5q0)
       print *, new(B_scaled, -0.5q0, -0.5q0)
  
       ! Print Difference
       print *, 'DIFFERENCE'
       print *, new(B_scaled, 0.5q0, 0.5q0) - ordinary(B_scaled, 0.5q0, 0.5q0)
  
       ! Contribution of New Part
       print *, 'CONTRIBUTION OF NEW PART'
       print *, 0.5q0 * FINE_STRUCTURE_CONST**2 * B_scaled**3 * NUCLEAR_CHARGE**(-7.0q0/2.0q0) * 0.5q0
       print *, 0.5q0 * FINE_STRUCTURE_CONST**2 * B_scaled**3 * NUCLEAR_CHARGE**(-7.0q0/2.0q0) * -0.5q0
    end do
  
  contains
  
    ! Function to calculate the ratio
    real(16) function ratio(Magnetic_Field, electron_spin, nuclear_spin)
      real(16), intent(in) :: Magnetic_Field, electron_spin, nuclear_spin
      ratio = (FINE_STRUCTURE_CONST**2 * Magnetic_Field**2 * PROTON_MASS_AU * electron_spin) / &
              (NUCLEAR_G_FACTOR * nuclear_spin * NUCLEAR_CHARGE**(-7.0q0/2.0q0))
    end function ratio
  
    ! Function to calculate the ordinary value
    real(16) function ordinary(Magnetic_Field, electron_spin, nuclear_spin)
      real(16), intent(in) :: Magnetic_Field, electron_spin, nuclear_spin
      ordinary = H_GROUND_STATE_ENERGY + (1.0q0 / (2.0q0 * PROTON_MASS_AU)) * NUCLEAR_G_FACTOR * nuclear_spin * Magnetic_Field + &
                 0.5q0 * ELECTRONIC_G_FACTOR * electron_spin * Magnetic_Field
    end function ordinary
  
    ! Function to calculate the new value with relativistic contribution
    real(16) function new(Magnetic_Field, electron_spin, nuclear_spin)
      real(16), intent(in) :: Magnetic_Field, electron_spin, nuclear_spin
      real(16) :: relativistic_contribution
      relativistic_contribution = 0.5q0 * FINE_STRUCTURE_CONST**2 * Magnetic_Field**3 * NUCLEAR_CHARGE**(-7.0q0/2.0q0) * electron_spin
      new = ordinary(Magnetic_Field, electron_spin, nuclear_spin) + relativistic_contribution
    end function new
  
  end program magnetic_interactions
  