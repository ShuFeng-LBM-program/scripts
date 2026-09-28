import math

class UnitConverter:
    def __init__(self, dx: float, rho_ref_phy: float, gamma: float):
        """
        Initializes conversion factors between Physical and Lattice units.
        
        :param dx: Physical cell size (lattice spacing) [m]
        :param rho_ref_phy: Reference physical density [kg/m^3]
        :param gamma: Heat capacity ratio
        """
        # Lattice Boltzmann Constants
        self.gamma = gamma
        self.cs_lb = math.sqrt(self.gamma / 3.0)
        self.cs2_lb = self.cs_lb**2
        self.cs_phy = 340.0  # Speed of sound in air [m/s]

        # Basic Factors
        self.factor_length = dx                          # C_l [m]
        self.factor_rho = rho_ref_phy                    # C_rho [kg/m^3]
        self.factor_vel = self.cs_phy / self.cs_lb       # C_u [m/s]
        
        # Derived Factors
        self.factor_time = self.factor_length / self.factor_vel  # C_t [s]
        self.factor_mass = self.factor_rho * (self.factor_length**3) # C_m [kg]
        
        # Physics Factors
        self.factor_kinematic_vis = (self.factor_length**2) / self.factor_time # C_nu [m^2/s]
        self.factor_pressure = self.factor_mass / (self.factor_length * (self.factor_time**2)) # C_p [Pa]
        self.factor_omega = self.factor_vel / self.factor_length

        print(f"UnitConverter initialized:")
        print(f"  C_l   = {self.factor_length:.6e}")
        print(f"  C_rho = {self.factor_rho:.6e}")
        print(f"  C_vel = {self.factor_vel:.6e}")
        print(f"  C_t   = {self.factor_time:.6e}")

    def phys_to_lb_velocity(self, u_phys):
        return u_phys / self.factor_vel

    def lb_to_phys_velocity(self, u_lb):
        return u_lb * self.factor_vel
    
    def step_to_time(self, step):
        return step * self.factor_time
    
    def time_to_step(self, time):
        return time / self.factor_time
    
    def densityLB_to_pressureLB(self, density_lb):
        return self.cs2_lb * density_lb
    
    def lb_to_phys_pressure(self, pressure_lb):
        return pressure_lb * self.factor_pressure
    
    def phys_to_lb_pressure(self, pressure_phys):
        return pressure_phys / self.factor_pressure
    
    def lb_to_phys_length(self, length_lb):
        return length_lb * self.factor_length
    
    def phys_to_lb_length(self, length_phys):
        return length_phys / self.factor_length

    def lb_to_phys_omega(self, omega_lb):
        return omega_lb * self.factor_omega

    def phys_to_lb_omega(self, omega_phys):
        return omega_phys / self.factor_omega

    # You can add similar helpers for pressure, viscosity, etc.
