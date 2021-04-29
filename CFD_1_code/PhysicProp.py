"""
Class PhysicProp.
- In charge of: storing the physical properties of the working fluid, typically dry air.
Physical units in International System.
"""


class PhysicProp:
    def __init__(self, a_fluid: str):
        # Physical properties of isentropic flows.
        if a_fluid == 'air':
            self.gamma = 1.4  # Gamma = cp/cv.
            self.R = 287.05  # Specific gas constant, in J/kg/K.
            self.cp = self.gamma * self.R * (self.gamma - 1)  # Specific heat coefficient, in J/kg/K.
            self.P0 = 101325  # Std. atmospheric pressure, in Pascal.
            self.T0 = 288.15  # Std. atmospheric temperature, in Kelvin.
            self.rho0 = self.P0 / self.R / self.T0  # Std. atmospheric density, in kg/m3.
            self.u0 = 10  # Flow x component velocity, in m/s.

        # Physical properties of non compressible flows.
        if a_fluid == 'water':
            pass
