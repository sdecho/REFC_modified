import numpy as np


class REFC:
    def __init__(self, dM_re, dM_e, dM_x, dM_cc, N):
        self.dM_re = dM_re  # mass of recharge magma at each step, >0
        self.dM_e = dM_e  # mass of erupted magma at each step, <0
        self.dM_x = dM_x  # mass of crystallization at each step, <0
        self.dM_cc = dM_cc  # mass of crustal simulation at each step, >0
        self.N = N  # number of steps of simulation
        self.M_ch, self.M_x, self.M_cc, self.M_re, self.M_e = self.total_mass()

    def total_mass(self):
        M_ch = np.ones((self.N, 1))
        M_x = np.zeros((self.N, 1))
        M_e = np.zeros((self.N, 1))
        M_re = np.zeros((self.N, 1))
        M_cc = np.zeros((self.N, 1))
        for i in range(1, self.N):
            M_ch[i] = M_ch[i - 1] + self.dM_re + self.dM_x + self.dM_cc + self.dM_e
            M_x[i] = M_x[i-1] + self.dM_x
            M_cc[i] = M_cc[i-1] +self.dM_cc
            M_re[i] = M_re[i - 1] + self.dM_re
            M_e[i] = M_e[i - 1] + self.dM_e
        return M_ch, M_x, M_cc, M_re, M_e

    def incompatible(self, c_re, c_cc, D, c0):
        c_ch = np.ones((self.N, 1))
        c_ch[0] = c0
        for i in range(1, self.N):
            c_x = c_ch[i-1] * D
            c_ch[i] = (c_re * self.dM_re + c_x * self.dM_x + c_cc * self.dM_cc + c_ch[i-1] * self.dM_e +
                    c_ch[i-1] * self.M_ch[i-1]) / self.M_ch[i]

        return c_ch

    def water_solubility(self, solubility, c_re, c_cc, c0):
        c_melt = np.ones((self.N, 1))  #water concentration in the melt
        c_fluid = np.zeros((self.N, 1)) #water mass fraction in the instant cumulate+vapor
        m_fluid = np.zeros((self.N, 1)) #accumulated water mass

        if c0 >= solubility:
            c_melt[0] = solubility
            c_fluid[0] = (c0 * (self.dM_e + self.dM_x + self.dM_cc + self.dM_re) - \
                         self.dM_re * c_re - self.dM_e * c_melt[0] - self.dM_cc * c_cc)/self.dM_x
        else:
            c_melt[0] = c0
            c_fluid[0] = 0
        m_fluid[0] = c_fluid[0] * self.dM_x/100

        for i in range(1, self.N):
            c_ch_try = (c_re * self.dM_re + c_fluid[i - 1] * self.dM_x + c_cc * self.dM_cc
                        + c_melt[i - 1] * self.dM_e +
                        c_melt[i - 1] * self.M_ch[i - 1]) / self.M_ch[i]
            if c_ch_try >= solubility:
                c_melt[i] = solubility
                c_fluid[i] = (c_melt[i] * (self.dM_e + self.dM_x + self.dM_cc + self.dM_re) - self.dM_re * c_re -
                          self.dM_e * c_melt[i] - self.dM_cc * c_cc)/self.dM_x
            else:
                c_melt[i] = c_ch_try
                c_fluid[i] = 0

            m_fluid[i] = m_fluid[i - 1] + c_fluid[i] * self.dM_x/100

        return c_melt, c_fluid, m_fluid

    def volatile_partition(self, kd, c_re, c_cc, c0, cfluid_x, m_fluid):
        """
        kd: volatile partition coefficient
        cfluid_x: mass fraction of fluid/(cumulate + fluid)
        m_fluid: accumulate mass of the fluid
        """
        c_ch = np.ones((self.N, 1)) #S concentration in the melt
        c_x = np.ones((self.N, 1)) #S concentration in the instant vapor+cumulate
        cS_fluid = np.ones((self.N, 1)) #S concentration in the accumulate fluid
        D_bulk = cfluid_x * kd/100  # cfluid_x fluid concentration in the cumulate+fluid
        c_ch[0] = c0
        c_x[0] = D_bulk[0] * c_ch[0]  # sulfur concentration in instant the cumulate + fluid
        if cfluid_x[0] == 0:
            cS_fluid[0] = 0  # sulfur concentration in the fluid
        else:
            cS_fluid[0] = kd * c_ch[0]/10000
        for i in range(1, self.N):
            c_ch[i] = (c_re * self.dM_re + c_x[i-1] * self.dM_x + c_cc * self.dM_cc + c_ch[i - 1] * self.dM_e +
                       c_ch[i - 1] * self.M_ch[i - 1]) / self.M_ch[i]
            c_x[i] = D_bulk[i] * c_ch[i]
            if m_fluid[i] ==0:
                cS_fluid[i] = 0
            else:
                cS_fluid[i] = (cS_fluid[i - 1] * m_fluid[i - 1] + c_x[i] * self.dM_x/10000) / m_fluid[i]
        return c_ch, c_x, cS_fluid

    def solubility_control(self, solubility, c_re, c_cc, c0):
        c_ch = np.ones((self.N, 1))
        c_x = np.ones((self.N, 1))
        c_ch[0] = c0
        if c0 == solubility:
            c_x[0] = (c0*(self.dM_e+self.dM_x+self.dM_cc+self.dM_re) - self.dM_re * c_re - self.dM_e * c0 -self.dM_cc * c_cc)/self.dM_x
        else:
            c_x[0] = 0
        for i in range(1, self.N):
            c_ch_try = (c_re * self.dM_re + c_x[i-1] * self.dM_x + c_cc * self.dM_cc + c_ch[i-1] * self.dM_e +
                        c_ch[i-1] * self.M_ch[i-1]) / self.M_ch[i]
            if c_ch_try >= solubility:
                c_ch[i] = solubility
                c_x[i] = (c_ch[i]*(self.dM_e+self.dM_x+self.dM_cc+self.dM_re) - self.dM_re * c_re -
                          self.dM_e * c_ch[i] - self.dM_cc * c_cc)/self.dM_x
            else:
                c_ch[i] = c_ch_try
                c_x[i] = 0
        return c_ch, c_x

    def chalcophile(self, Dsf, c_re, c_cc, c0, cs_x):
        c_ch = np.ones((self.N, 1))
        c_x = np.ones((self.N, 1))
        D_bulk = cs_x * Dsf/363636
        c_ch[0] = c0
        c_x[0] = D_bulk[0] * c_ch[0]
        for i in range(1, self.N):
            c_ch[i] = (c_re * self.dM_re + c_x[i-1] * self.dM_x + c_cc * self.dM_cc + c_ch[i - 1] * self.dM_e +
                       c_ch[i - 1] * self.M_ch[i - 1]) / self.M_ch[i]
            c_x[i] = D_bulk[i] * c_ch[i]
        return c_ch, c_x

    def ferric_iron(self, DFe3, c_re, c_cc, c0, cs_mt):
        c_ch = np.ones((self.N, 1))
        c_x = np.ones((self.N, 1))
        D_bulk = np.ones((self.N, 1))

        c_ch[0] = c0
        D_bulk[0] = cs_mt[0] * DFe3[0]
        c_x[0] = D_bulk[0] * c_ch[0]

        for i in range(1, self.N):
            D_bulk[i] = cs_mt[i] * DFe3[i]
            c_x[i] = D_bulk[i] * c_ch[i-1]
            c_ch[i] = (c_re * self.dM_re + c_x[i-1] * self.dM_x + c_cc * self.dM_cc + c_ch[i - 1] * self.dM_e +
                       c_ch[i - 1] * self.M_ch[i - 1]) / self.M_ch[i]
            c_x[i] = D_bulk[i] * c_ch[i]
        return c_ch, c_x









