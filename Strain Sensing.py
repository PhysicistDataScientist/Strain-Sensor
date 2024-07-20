# Comments:
'''
None.
'''

# Libraries:
from os import listdir
from os.path import join, exists
from chardet import detect
from numpy import ones, array, linspace, exp, log10, tan, pi
from scipy.signal import find_peaks
from pandas import read_csv, read_excel, DataFrame, ExcelWriter
from matplotlib.pyplot import subplots, show

# Classes:
class Capillary:
    # Initialize the object:
    def __init__(self, file, material_core, material_clad):
        self.file = file
        self.band = file.split('-')[0].replace('Band ','').rstrip()
        self.material_core, self.material_clad = material_core, material_clad
        # Select the refractive index of the core:
        if self.material_core == 'air':
            self.n_core = 1.000
        else:
            self.n_core = 1.333
        # Select the refractive index dispersion relation:
        if self.material_clad == 'PMMA':
            self.n = lambda wv: (1 + 1.1819 * wv**2 / (wv**2 - 0.011313))**0.5
        else:
            self.n = lambda wv: (1 + (0.6961663 * wv**2 / (wv**2 - 0.0684043**2)) + (0.4079426 * wv**2 / (wv**2 - 0.1162414**2)) + (0.8974794 * wv**2 / (wv**2 - 9.896161**2)))**0.5
        # Function of theoretical minimum wavelengths:
        self.wv_min = lambda wv, m, t: 2 * t / m * (self.n(wv)**2 - self.n_core**2)**0.5 # Thickness and wavelength shall be in um.
        # Initializations:
        self.df_t, self.df_spec = DataFrame(), DataFrame()
        self.D_ext, self.t, self.R_core, self.L, self.wv_cal, self.m_cal = 0, 0, 0, 0, 0, 0
    # See the thickness variation:
    def thickness_variation(self):
        self.df_t = read_excel(join(sample_folder, f'Band {self.band} - Table - Transverse Section.xlsx'))
        self.df_t['t_mean [um]'] = round(self.df_t.loc[:, 't1 [um]' : 't4 [um]'].sum(axis=1) / 4, 2)
        self.df_t['ut [um]'] = self.df_t.loc[:, 't1 [um]':'t4 [um]'].std(ddof=1, axis=1)
        t_mean_mean = self.df_t['t_mean [um]'].mean()
        s_mean_mean = self.df_t.loc[:, 't_mean [um]'].std(ddof=1, axis=0)
        fig, ax = subplots(nrows=1, ncols=1, layout='constrained', figsize=(12,4))
        ax.set_xlabel('Position [cm]', loc='center', fontsize=14)
        ax.set_ylabel('Thickness [um]', loc='center', fontsize=14)
        ax.errorbar(self.df_t['z [cm]'], self.df_t['t_mean [um]'], yerr=self.df_t['ut [um]'], marker='.', ms=20, mfc='gray', mec='black', linestyle='none', lw=2, ecolor='black', capsize=5, label='Data')
        ax.plot(array([min(self.df_t['z [cm]']), max(self.df_t['z [cm]'])]), t_mean_mean * ones(2), c='red', linestyle='dashed', lw=2, label=r'$\overline{t}$ = ' + f'({round(t_mean_mean, 2)}' + r' $\pm$ '+ f'{round(s_mean_mean, 2)})' + r'$\mu$m')
        ax.legend(loc='best', fontsize=12)
        show()
        if save:
            fig.savefig(join(sample_folder, f'Band {self.band} - Graphic - Thickness Variation.png'))
    # Ask the fiber parameters:
    def get_parameters(self):
        # Ask the fiber length:
        while True:
            question = input('Fiber length [cm]: [8.5] ')
            try:
                if not question.replace('.', '').isdigit():
                    raise ValueError('Not a valid number!')
                L = float(question)
                self.L = L * 1e4 # Convert cm -> um.
            except ValueError as ve:
                print(ve)
            else:
                break
        # Ask the external diameter:
        while True:
            question = input('External diameter [um]: [200] ')
            try:
                if not question.replace('.', '').isdigit():
                    raise ValueError('Not a valid number!')
                self.D_ext = float(question)
            except ValueError as ve:
                print(ve)
            else:
                break
        # Ask the cladding thickness:
        while True:
            question = input('Cladding thickness [um]: [24.74] ')
            try:
                if not question.replace('.', '').isdigit():
                    raise ValueError('Not a valid number!')
                self.t = float(question)
            except ValueError as ve:
                print(ve)
            else:
                break
    # Find a valley for calibration:
    def calibrate_valley(self):
        self.df_spec = read_csv(join(sample_folder, self.file))
        self.df_spec = self.df_spec[self.df_spec['Wavelength [nm]'] >= 450] 
        self.df_spec = self.df_spec[self.df_spec['Transmission [dB]'] >= -100]
        # See the complete spectrum:
        print(f'Experimental spectrum of band {self.band}:')
        fig, ax = subplots(nrows=1, ncols=1, layout='constrained', figsize=(12,4))
        ax.set_xlabel('Wavelength [nm]', loc='center', fontsize=14)
        ax.set_ylabel('Transmission [dB]', loc='center', fontsize=14)
        ax.scatter(self.df_spec['Wavelength [nm]'], self.df_spec['Transmission [dB]'], s=20, c='black')
        show()
        if save:
            fig.savefig(join(sample_folder, f'Band {self.band} - Graphic - Experimental Spectrum.png'))
        # Find a valley:
        satisfied = False
        while not satisfied:
            # Ask the width of the valleys and the endpoints:
            while True:
                question1 = input('Choose the width of the valleys [nm]: [8] ')
                question2 = input('Pass the wavelength endpoints [nm]: [1400,1450] ')
                try:
                    # Width:
                    if not question1.replace('.','').isdigit():
                        raise ValueError('Not a number!')
                    width = float(question1)
                    # Endpoints:
                    if question2.count(',') != 1:
                        raise ValueError('Invalid number of parameters!')
                    aux_list = question2.split(',')
                    endpoints = list() # nm.
                    for aux in aux_list:
                        if not aux.replace('.','').isdigit():
                            raise ValueError('Not a valid number!')
                        endpoints.append(float(aux))
                    if endpoints[0] < 350:
                        raise ValueError('It must be greater than 350nm!')
                    if endpoints[1] > 1750:
                        raise ValueError('It must be smaller than 1750nm!')
                    # Select the data from the data frame:
                    df_aux = self.df_spec[self.df_spec['Wavelength [nm]'] >= endpoints[0]]
                    df_aux = df_aux[df_aux['Wavelength [nm]'] <= endpoints[1]]
                    # Find the valleys:
                    valleys, _ = find_peaks(-df_aux['Transmission [dB]'], width=width)
                    if len(valleys) == 0:
                        raise ValueError('No valley found!')
                except ValueError as ve:
                    print(ve)
                else:
                    break
            # Get the valleys:
            wavelengths = df_aux.loc[:, 'Wavelength [nm]'].values # nm.
            minimums = wavelengths[valleys] # nm.
            # Visualize the selected region:
            fig, ax = subplots(nrows=1, ncols=1, layout='constrained', figsize=(12,4))
            ax.set_xlabel('Wavelength [nm]', loc='center', fontsize=14)
            ax.set_ylabel('Transmission [dB]', loc='center', fontsize=14)
            ax.scatter(df_aux['Wavelength [nm]'], df_aux['Transmission [dB]'], s=20, c='black')
            for minimum in minimums:
                ax.plot(minimum * ones(2), array([min(df_aux['Transmission [dB]']), max(df_aux['Transmission [dB]'])]), lw=2, c='red', linestyle='dashed')
            show()
            print('You can choose on from the following valleys for calibration:')
            [print(f'{minimum}nm') for minimum in minimums]
            print()
            # Ask the valley to be chosen:
            while True:
                question = input('Pass the valley for calibration [nm]: [1419.6] ')
                try:
                    if not question.replace('.', '').isdigit():
                        raise ValueError('Not a valid number!')
                    minimum = float(question)
                    if minimum not in minimums:
                        raise ValueError('Not a valley!')
                    # Visualize this choice:
                    fig, ax = subplots(nrows=1, ncols=1, layout='constrained', figsize=(12,4))
                    ax.set_xlabel('Wavelength [nm]', loc='center', fontsize=14)
                    ax.set_ylabel('Transmission [dB]', loc='center', fontsize=14)
                    ax.scatter(df_aux['Wavelength [nm]'], df_aux['Transmission [dB]'], s=20, c='black')
                    ax.plot(minimum * ones(2), array([min(df_aux['Transmission [dB]']), max(df_aux['Transmission [dB]'])]), lw=2, c='red', linestyle='dashed')
                    show()
                except ValueError as ve:
                    print(ve)
                else:
                    break   
            # Ask if the user is satisfied with the choice for calibration:
            while True:
                question = input(f'Are you satisfied with the choice of {minimum}nm? [y/n]')
                try:
                    if question not in ['y', 'n']:
                        raise ValueError('Invalid answer!')
                    if question == 'y':
                        satisfied = True
                except ValueError as ve:
                    print(ve)
                else:
                    break
        minimum *= 1e-3 # Convert nm -> um.
        self.wv_cal = minimum
        # Save the graphic:
        if save:
            fig.savefig(join(sample_folder, f'Band {self.band} - Graphic - Calibration, wv = {round(self.wv_cal * 1e3, 1)}nm.png'))
    # Analyse the data:
    def analyse_data(self):
        # Calculate the m-th order of the valley given the thickness:
        self.m_cal = round(2 * self.t / self.wv_cal * (self.n(self.wv_cal)**2 - self.n_core**2)**0.5) # m for calibration.
        satisfied = False
        while not satisfied:
            ## Ask the nearby m orders to calculate the valleys:
            while True:
                question = input('Pass the number of nearby m orders (m_left,m_right): [5,5] ')
                try:
                    if question.count(',') != 1:
                        raise ValueError('Invalid number of parameters!')
                    aux_list = question.split(',')
                    m_limits = list()
                    for aux in aux_list:
                        if not aux.isdigit():
                            raise ValueError('Not an integer!')
                        m_limits.append(int(aux))
                except ValueError as ve:
                    print(ve)
                else:
                    break
            m_left, m_right = m_limits
            m_list = [m for m in range(self.m_cal - m_left, self.m_cal + m_right + 1)]
            # Find the theoretical peaks:
            minimums_the = list()
            N = 20 # Number of iterations of the numerical equation solver.
            for m in m_list:
                # Solve numerically:
                minimum = self.wv_cal # Initial guess [um].
                for i in range(N):
                    minimum = self.wv_min(minimum, m, self.t)
                minimums_the.append(minimum)
            minimums_the = array(minimums_the, dtype=float)
            df_min = DataFrame(data={'Valley [nm]': minimums_the * 1e3}) # Convert um -> nm.
            # The entire spectrum in wavelength:
            fig, ax = subplots(nrows=1, ncols=1, layout='constrained', figsize=(12,4))
            ax.set_xlabel('Wavelength [nm]', loc='center', fontsize=12)
            ax.set_ylabel('Transmission [dB]', loc='center', fontsize=12)
            ax.scatter(self.df_spec['Wavelength [nm]'], self.df_spec['Transmission [dB]'], c='black', s=30)
            ind_min = 0
            for i, (m, minimum) in enumerate(zip(m_list, minimums_the)):
                if m > self.m_cal:
                    color = 'blue'
                elif m < self.m_cal:
                    color = 'red'
                else:
                    color = 'green'
                    ind_min = i
                ax.plot(minimum * 1e3 * ones(2), array([min(self.df_spec['Transmission [dB]']), max(self.df_spec['Transmission [dB]'])]), c=color, linewidth=2, linestyle='dashed')
            show()
            # Select a region of the spectrum:
            ## Ask the endpoints of the spectrum in wavelength:
            while True:
                question = input('Pass the wavelength endpoints [nm]: [1300,1700] ')
                try:
                    if question.count(',') != 1:
                        raise ValueError('Invalid number of parameters!')
                    aux_list = question.split(',')
                    x_limits = list() # nm.
                    for aux in aux_list:
                        if not aux.replace('.','').isdigit():
                            raise ValueError('Not a valid number!')
                        x_limits.append(float(aux))
                    if x_limits[0] < 350:
                        raise ValueError('It must be greater than 350nm!')
                    if x_limits[1] > 1750:
                        raise ValueError('It must be smaller than 1750nm!')
                except ValueError as ve:
                    print(ve)
                else:
                    break
            # Select the data from the data frame:
            df_aux = self.df_spec[self.df_spec['Wavelength [nm]'] >= x_limits[0]]
            df_aux = df_aux[df_aux['Wavelength [nm]'] <= x_limits[1]]
            ## Ask the graphic y limits:
            while True:
                question = input('Pass the y limits (bottom, top) [nm]: [-80,-55] ')
                try:
                    if question.count(',') != 1:
                        raise ValueError('Invalid number of parameters!')
                    aux_list = question.split(',')
                    y_limits = list() # nm.
                    for aux in aux_list:
                        if not aux.replace('-','').replace('.','').isdigit():
                            raise ValueError('Not a valid number!')
                        y_limits.append(float(aux))
                    if y_limits[0] > y_limits[1]:
                        raise ValueError('The bottom must be smaller than the top!')
                except ValueError as ve:
                    print(ve)
                else:
                    break
            ## Plot the selected region:
            fig1, ax1 = subplots(nrows=1, ncols=1, layout='constrained', figsize=(12,4))
            ax1.set_xlabel('Wavelength [nm]', loc='center', fontsize=12)
            ax1.set_ylabel('Transmission [dB]', loc='center', fontsize=12)
            ax1.set_xlim(left=x_limits[0], right=x_limits[1])
            ax1.set_ylim(bottom=y_limits[0], top=y_limits[1])
            ax1.scatter(df_aux['Wavelength [nm]'], df_aux['Transmission [dB]'], c='black', s=30)
            for i, minimum in enumerate(minimums_the):
                if i > ind_min:
                    color = 'blue' 
                elif i < ind_min:
                    color = 'red'
                else:
                    color = 'green'
                ax1.plot(minimum * 1e3 * ones(2), array([min(df_aux['Transmission [dB]']), max(df_aux['Transmission [dB]'])]), c=color, linewidth=2, linestyle='dashed')
            show()
            ## Change wavelength to frequency:
            c = 3e14 # Speed of the light in vaccumm [um/s]. 
            df_aux['Frequency [Hz]'] = c / (df_aux['Wavelength [nm]'] * 1e-3) # The wavelength must be in um.
            minimum_freq_list = [c / minimum for minimum in minimums_the]
            f_beg, f_end = c / (x_limits[1] * 1e-3), c / (x_limits[0] * 1e-3) # Endpoints of frequency [Hz].
            ## Plot the select region in frequency:
            fig2, ax2 = subplots(nrows=1, ncols=1, layout='constrained', figsize=(12,4))
            ax2.set_xlabel('Frequency [THz]', loc='center', fontsize=12)
            ax2.set_ylabel('Transmission [dB]', loc='center', fontsize=12)
            ax2.set_xlim(left=f_beg * 1e-12, right=f_end * 1e-12)
            ax2.set_ylim(bottom=y_limits[0], top=y_limits[1])
            ax2.scatter(df_aux['Frequency [Hz]'] * 1e-12, df_aux['Transmission [dB]'], c='black', s=30)
            for i, minimum_freq in enumerate(minimum_freq_list):
                if i > ind_min:
                    color = 'blue' 
                elif i < ind_min:
                    color = 'red'
                else:
                    color = 'green'
                ax2.plot(minimum_freq * ones(2) * 1e-12, array([min(df_aux['Transmission [dB]']), max(df_aux['Transmission [dB]'])]), c=color, linewidth=2, linestyle='dashed')
            show() 
            # Ask if the user is satisfied with the plots:
            while True:
                question = input(f'Are you satisfied with the plots? [y/n]')
                try:
                    if question not in ['y', 'n']:
                        raise ValueError('Invalid answer!')
                    if question == 'y':
                        satisfied = True
                except ValueError as ve:
                    print(ve)
                else:
                    break
        # Save the graphics:
        if save:
            with ExcelWriter(join(sample_folder, f'Band {self.band} - Table - Valleys Formula, wv_cal = {round(self.wv_cal * 1e3, 1)}nm, m = {self.m_cal}, t = {round(self.t, 2)}um.xlsx'), engine='openpyxl') as writer:
                df_min.to_excel(writer, index=False)
            fig1.savefig(join(sample_folder, f'Band {self.band} - Graphic - Spectrum Wavelength, wv = {round(self.wv_cal * 1e3, 1)}nm, m = {self.m_cal}, t = {round(self.t, 2)}um.png'))
            fig2.savefig(join(sample_folder, f'Band {self.band} - Graphic - Spectrum Frequency, wv = {round(self.wv_cal * 1e3, 1)}nm, m = {self.m_cal}, t = {round(self.t, 2)}um.png'))
    # Given a set of m-th for the calibration valley, calculate the thickness the capillary should have for each order:
    def thickness_from_order(self):
        # Ask the m-th orders to be analysed:
        print(f'Current m = {self.m_cal}')
        while True:
            question = input(f'm-th orders for the calibration valley of {round(self.wv_cal * 1e3, 1)}nm: [20,38,40] ')
            try:
                m_list = list()
                if ',' not in question:
                    if len(question) != 0:
                        if not question.isdigit():
                            raise ValueError('Not a non negative integer!')
                        m = int(question)
                        if m == 0:
                            raise ValueError('Null order!')
                        m_list.append(m)
                else:
                    aux_list = question.split(',')
                    for aux in aux_list:
                        if not aux.isdigit():
                            raise ValueError('Not a non negative integer!')
                        m = int(aux)
                        if m == 0:
                            raise ValueError('Null order!')
                        m_list.append(m)
            except ValueError as ve:
                print(ve)
            else:
                break
        # Calculte the corresponding thickness for the calibration wavelength:
        for m in m_list:
            t = m * self.wv_cal / (2 * (self.n(self.wv_cal)**2 - self.n_core**2)**0.5)
            print(f'm = {m}: t = {round(t, 2)}um')
        print()
    # Simulate the transmittance spectrum:
    def simulate(self):
        print('Let\'s get the spectrum by simulation:')
        # Ask for the wavelength endpoints:
        while True:
            question = input('Wavelength endpoints [nm] and number of points: [1300,1700,1000] ')
            try:
                if question.count(',') != 2:
                    raise ValueError('Incompatible number of parameters!')
                list_aux = question.split(',')
                parameters = list()
                for aux in list_aux:
                    if not aux.replace('.','').isdigit():
                        raise ValueError('Invalid parameter!')
                    parameters.append(float(aux))
                if parameters[2] == 0:
                    raise ValueError('Null number of points!')
                parameters[2] = int(parameters[2])
            except ValueError as ve:
                print(ve)
            else:
                break
        wavelength = linspace(parameters[0], parameters[1], parameters[2]) * 1e-3 # Convert nm -> um.
        # Get the valleys by simulating the transmittance:
        ## Parameters:
        j01 = 2.405 # First bessel function zero.
        self.R_core = self.D_ext / 2 - self.t # Core radius [um].
        ## Pre calculations:
        k0 = 2 * pi / wavelength # Wavenumber [um^-1].
        n_clad = self.n(wavelength) # Cladding refractive index.
        phi = k0 * self.t * (n_clad**2 - self.n_core**2)**0.5
        ## Eliminate multiples of 2 * pi:
        mult = phi // (2 * pi) # Multiple of 2 * pi.
        phi = phi - mult * (2 * pi) # Angle in [0, 2 * pi].
        eps = (n_clad / self.n_core)**2
        ## Calculations:
        alpha = ((1 + (1 / tan(phi))**2) / (eps - 1)) * (j01**3 / (k0**3 * self.n_core**3 * self.R_core**4)) * ((eps**2 + 1) / 2) # [um^-1].
        P_dB = 10 * log10(exp(-alpha * self.L)) # Power [dB].
        ## Valleys:
        valleys, _ = find_peaks(-P_dB)
        minimums = wavelength[valleys] # Minimum transmittance wavelength [um].
        df = DataFrame(data = {'Valley [nm]': minimums * 1e3}) # Convert um -> nm for the valleys.
        ## Plot:
        y_bottom, y_top = -0.1, 0
        fig, ax = subplots(nrows=1, ncols=1, layout='constrained', figsize=(12,4))
        ax.set_xlabel('Wavelength [nm]', loc='center', fontsize=14)
        ax.set_ylabel('Power [dB]', loc='center', fontsize=14)
        ax.set_ylim(y_bottom, y_top)
        ax.plot(wavelength * 1e3, P_dB, c='black', lw=2)
        for minimum in minimums:
            ax.plot(minimum * 1e3 * ones(2), array([y_bottom, y_top]), c='red', linestyle='dashed')
        show()
        # Save the graphic and the data frame:
        if save:
            fig.savefig(join(sample_folder, f'Band {self.band} - Graphic - Simulation Transmittance Valleys.png'))
            with ExcelWriter(join(sample_folder, f'Band {self.band} - Table - Valleys Simulation, wv_cal = {round(self.wv_cal * 1e3, 1)}nm, m = {self.m_cal}, t = {round(self.t, 2)}um.xlsx'), engine='openpyxl') as writer:
                df.to_excel(writer, index=False)
        # Apply strain to change the transmittance:
        ## Ask the axial strains to analyse:
        while True:
            question = input('Axial strains to be analysed [%]: [1,2,3] ')
            try:
                strain_axial_list = list()
                if ',' not in question:
                    if not question.replace('.','').isdigit():
                        raise ValueError('Not a valid number!')
                    strain_axial = float(question) / 100 # Convert % -> decimal.
                    strain_axial_list.append(strain_axial)
                else:
                    aux_list = question.split(',')
                    for aux in aux_list:
                        if not aux.replace('.','').isdigit():
                            raise ValueError('Not a valid number!')
                        strain_axial = float(aux) / 100 # Convert % -> decimal.
                        strain_axial_list.append(strain_axial)
            except ValueError as ve:
                print(ve)
            else:
                break
        # Choose the formula to calculate the variations due strain:
        print('You can choose one of the following formulas to calculare variation due strain:')
        print('- Complete: dD = -D * [1 - (1 + s_a)^-poisson]')
        print('- Approximate: dD = - D * poisson * s_a')
        while True:
            question = input('Choose the formula to calculate variations: [complete/approximate] ')
            try:
                if question not in ['complete', 'approximate']:
                    raise ValueError('Invalid answer!')
                formula = question
            except ValueError as ve:
                print(ve)
            else:
                break
        ## Transmittance:
        fig1, ax1 = subplots(nrows=1, ncols=1, layout='constrained', figsize=(12,4))
        ax1.set_xlabel('Wavelength [nm]', loc='center', fontsize=14)
        ax1.set_ylabel('Power [dB]', loc='center', fontsize=14)
        ax1.set_ylim(y_bottom, y_top)
        ax1.plot(wavelength * 1e3, P_dB, lw=2, label='0.0')
        ## Refractive index:
        fig2, ax2 = subplots(nrows=1, ncols=1, layout='constrained', figsize=(12,4))
        ax2.set_xlabel('Wavelength [nm]', loc='center', fontsize=14)
        ax2.set_ylabel('Cladding refractive index', loc='center', fontsize=14)
        ax2.plot(wavelength * 1e3, n_clad, lw=2, label='0.0')
        ## For each strain:
        for s_a in strain_axial_list:
            n_clad = self.n(wavelength)
            # Get the material strain parameters:
            if self.material_clad == 'PMMA':
                poisson = 0.34
                p11, p12 = 0.298, 0.294
            else:
                poisson = 0.15
                p11, p12 = 0.121, 0.270
            # Calculate the variations due to the axial strain:
            if formula == 'complete':
                dt = -self.t * (1 - (1 + s_a)**(-poisson))
                dR_core = -self.R_core * (1 - (1 + s_a)**(-poisson))
            else:
                dt = -self.t * poisson * s_a
                dR_core = -self.R_core * poisson * s_a
            dn = - n_clad**3 * s_a / 2 * ((1 - poisson) * p12 - poisson * p11)
            # Update the variables:
            t = self.t + dt
            R_core = self.R_core + dR_core
            n_clad += dn
            phi = k0 * t * (n_clad**2 - self.n_core**2)**0.5
            mult = phi // (2 * pi)
            phi = phi - mult * (2 * pi)
            alpha = ((1 + (1 / tan(phi))**2) / (eps - 1)) * (j01**3 / (k0**3 * self.n_core**3 * R_core**4)) * ((eps**2 + 1) / 2)
            P_dB = 10 * log10(exp(-alpha * self.L))
            # Plot:
            ax1.plot(wavelength * 1e3, P_dB, lw=2, label=f'{s_a * 1e2}')
            ax2.plot(wavelength * 1e3, n_clad, lw=2, label=f'{s_a * 1e2}')
        ax1.legend(title=r'$\epsilon_{axial}$' + ' [%]', loc='best', fontsize=12, title_fontsize=12)
        ax2.legend(title=r'$\epsilon_{axial}$' + ' [%]', loc='best', fontsize=12, title_fontsize=12)
        show()
        if save:
            fig1.savefig(join(sample_folder, f'Band {self.band} - Graphic - Simulation Transmittance Strain, D_ext = {round(self.D_ext, 2)}um, t = {round(self.t, 2)}um, {formula.capitalize()}.png'))
            fig2.savefig(join(sample_folder, f'Band {self.band} - Graphic - Simulation Refractive Index Strain, D_ext = {round(self.D_ext, 2)}um, t = {round(self.t, 2)}um, {formula.capitalize()}.png'))
        # Analyse the variations due the different formulas:
        ## Ask the strains:
        while True:
            question = input('Final strain [%] and the number of points: [10,1000] ')
            try:
                if question.count(',') != 1:
                    raise ValueError('Incorrect number of parameters!')
                aux_list = question.split(',')
                if not aux_list[0].replace('.','').isdigit():
                    raise ValueError('Invalid limit!')
                s_a_end = float(aux_list[0])
                if not aux_list[1].isdigit():
                    raise ValueError('Invalid number of points!')
                N = int(aux_list[1])
            except ValueError as ve:
                print(ve)
            else:
                break
        # Create the strain points:
        strain_points = linspace(0, s_a_end, N)
        strain_points *= 1e-2 # Convert % to decimal.
        # Define the formulas to be analysed:
        f_complete = lambda D, s_a: -D * (1 - (1 + s_a)**-poisson)
        f_approx = lambda D, s_a: -D * poisson * s_a
        # Create the the variations:
        dD_ext_complete, dD_ext_approx = f_complete(self.D_ext, strain_points), f_approx(self.D_ext, strain_points) 
        dt_complete, dt_approx = f_complete(self.t, strain_points), f_approx(self.t, strain_points) 
        # Plot:
        ## External diameter:
        fig1, ax1 = subplots(nrows=1, ncols=1, layout='constrained', figsize=(12,4))
        ax1.set_xlabel('Strain [%]', loc='center', fontsize=14)
        ax1.set_ylabel(r'$\Delta D_{ext}$ $[\mu m]$', loc='center', fontsize=14)
        ax1.plot(strain_points * 1e2, dD_ext_complete, lw=2, label='Complete')
        ax1.plot(strain_points * 1e2, dD_ext_approx, lw=2, label='Approximate')
        ax1.legend(title=r'$D_{ext}$ = ' + str(round(self.D_ext, 2)) + r'$\mu m$', loc='best', fontsize=12, title_fontsize=12)
        show()
        ## Thickness:
        fig2, ax2 = subplots(nrows=1, ncols=1, layout='constrained', figsize=(12,4))
        ax2.set_xlabel('Strain [%]', loc='center', fontsize=14)
        ax2.set_ylabel(r'$\Delta t$ $[\mu m]$', loc='center', fontsize=14)
        ax2.plot(strain_points * 1e2, dt_complete, lw=2, label='Complete')
        ax2.plot(strain_points * 1e2, dt_approx, lw=2, label='Approximate')
        ax2.legend(title=f't = {round(self.t, 2)}' + r'$\mu m$', loc='best', fontsize=12, title_fontsize=12)
        show()
        ## Save the graphics:
        if save:
            fig1.savefig(join(sample_folder, f'Band {self.band} - Graphic - Simulation Strain Diameter, D_ext = {round(self.D_ext, 2)}um, t = {round(self.t, 2)}um.png'))
            fig2.savefig(join(sample_folder, f'Band {self.band} - Graphic - Simulation Strain Thickness, D_ext = {round(self.D_ext, 2)}um, t = {round(self.t, 2)}um.png'))
        # Calibration valley displacement:
        if formula == 'complete':
            dt = f_complete(self.t, strain_points)
        else:
            dt = f_approx(self.t, strain_points)
        t = self.t + dt
        n_clad0 = self.n(self.wv_cal)
        dn = - n_clad0**3 * strain_points / 2 * ((1 - poisson) * p12 - poisson * p11)
        n_clad = n_clad0 + dn
        ## Valley displacements:
        wv_m_t = 2 * t / self.m_cal * (n_clad0**2 - self.n_core**2)**0.5
        wv_m_n = 2 * self.t / self.m_cal * (n_clad**2 - self.n_core**2)**0.5
        wv_m_tn = 2 * t / self.m_cal * (n_clad**2 - self.n_core**2)**0.5
        # Plot:
        legend_title = r'$\lambda_{%d}$' % self.m_cal
        legend_title += f' = {round(self.wv_cal * 1e3, 2)}nm'
        fig, ax = subplots(nrows=1, ncols=1, layout='constrained', figsize=(12,4))
        ax.set_xlabel('Strain [%]', loc='center', fontsize=14)
        ax.set_ylabel(r'$\lambda_{%d}$ [nm]' % self.m_cal, loc='center', fontsize=14)
        ax.plot(strain_points * 1e2, wv_m_t * 1e3, lw=2, label=r'$\Delta t$ effect')
        ax.plot(strain_points * 1e2, wv_m_n * 1e3, lw=2, label=r'$\Delta n$ effect')
        ax.plot(strain_points * 1e2, wv_m_tn * 1e3, lw=2, label='Both effects')
        ax.legend(title=legend_title, loc='best', fontsize=12, title_fontsize=12)
        show()
        if save:
            fig.savefig(join(sample_folder, f'Band {self.band} - Graphic - Simulation Strain Valley, D_ext = {round(self.D_ext, 2)}um, t = {round(self.t, 2)}um.png'))
    # Show the capillary info:
    def show_info(self):
        print('---------------------------------')
        print(f'Band: {self.band}')
        print(f'Core: {self.material_core}')
        print(f'n_core = {self.n_core}')
        print(f'R_core = {round(self.R_core, 2)}um')
        print(f'Cladding: {self.material_clad}')
        print(f'L = {round(self.L * 1e-4, 1)}cm')
        print(f'D_ext = {round(self.D_ext, 2)}um')
        print(f't = {round(self.t, 2)}um')
        print(f'Calibration: lambda_{self.m_cal} = {round(self.wv_cal * 1e3, 1)}nm')
        print('---------------------------------\n')

# Function:
## Introduction to the program:
def intro():
    print('--------------------------------')
    print('Welcome to "Capillary Analyser"!')
    print('--------------------------------')
    print()
    print('Here you will be able to anaylse the spectrum of capillaries.')
    print()
## Find the path to the sample's files: 
def get_sample_folder():
    global folder
    while True:
        with open("Path Manager.txt", mode='rt') as f:
            folder = f.read()
            f.close()
        try: 
            if not exists(folder):
                raise FileExistsError('There is no folder with such path!')
            else:
                break
        except FileExistsError as fee:
            print(fee)
            print('Go to the "Path Manager.txt" and change its content!')
            input('Type anything when you are done: ')
    # Get the sample's folder path:
    global n_sample; global sample_folder
    while True:
        # Get the sample's number:
        while True:
            n_sample = input('Sample\'s number: ')
            try:
                if not n_sample.lstrip('-').isdigit():
                    raise ValueError('It is not a number!')
                elif '-' in n_sample:
                    raise ValueError('Negative integers are not allowed!')
                elif ',' in n_sample:
                    raise ValueError('Not integer numbers are not allowed!')
                n_sample = int(n_sample)
                if n_sample == 0:
                    raise ValueError('There is no fiber label 0!')
                else:
                    break
            except ValueError as ve:
                print(ve)
        # Join the folder with the sample identification suffix:
        sample_folder = join(folder, f'E - {n_sample}')
        try:
            if not exists(sample_folder):
                raise FileExistsError("Invalid fiber type or sample's number!")
            else:
                break
        except FileExistsError as fee:
            print(fee)
            print('Try again.')
# Ask if the user wants to save the graphics and data frames:
def ask_save():
    global save
    while True:
        question = input('Do you want to save graphics and data frames?: [y/n] ')
        try:
            if question not in ['y', 'n']:
                raise ValueError('Invalid answer!')
            elif question == 'y': 
                save = True
            else: 
                save = False 
        except ValueError as ve: 
            print(ve)
        else: 
            break
## Data file modification:
def modify_files():
    # Access files withing sample folder:
    files_list = [file for file in listdir(sample_folder) if 'Data - Spectrum' in file]
    possible_bands = [file.split('-')[0].replace('Band ','').rstrip() for file in files_list]
    print('You can analyse the following files:')
    [print(band) for band in possible_bands]
    print()
    # Ask the bands to be analysed:
    chosen_bands = list()
    while True:
        question = input('List the bands to be analysed: [A,B,C] ')
        try:
            if ',' not in question:
                if question not in possible_bands:
                    raise ValueError('Invalid band!')
                chosen_bands.append(question)
            else:
                aux_list = question.split(',')
                for aux in aux_list:
                    if aux not in possible_bands:
                        raise ValueError('Invalid band!') 
                    chosen_bands.append(aux)
            if len(chosen_bands) == 0:
                raise ValueError('No band was selected!')
        except ValueError as ve:
            print(ve)
        else:
            break
    print('You chose to analyse the following bands:')
    [print(band) for band in chosen_bands]
    print()
    # Select the paths to analyse:
    new_files_list = list()
    for file in files_list:
        band = file.split('-')[0].replace('Band ','').rstrip()
        if band in chosen_bands:
            new_files_list.append(file)
    files_list = new_files_list
    # Modify the spectrums:
    header = 'Wavelength [nm],Transmission [dB]\n'
    for file in files_list:
        data_path = join(sample_folder, file)
        file_name = data_path.split('E - ')[1].split('\\')[1]
        # Get the file enconding:
        with open(data_path, mode='rb') as f:
            detection = detect(f.read())
            f.close()
        # Read the data:
        with open(data_path, mode='rt', encoding=detection['encoding']) as f:
            rows = f.readlines()
            f.close()
        ## The file need to be modified:
        if '00\t\n' in rows:
            # Find the header position:
            ind = [i for i, row in enumerate(rows) if row == '00\t\n']
            ind = ind[0]
            # Cut off the data above the header and substitute by the new header:
            data = str()
            for i, row in enumerate(rows):
                if 'CTRWL' in row:
                    break
                elif i > ind:
                    data += row.replace(' ', '')
            data = header + data
            # Modify the old file:
            with open(data_path, mode='wt') as g:
                g.write(data)
                g.close()
            print(f'Modification on "{file_name}" is done.')
        ## The file was already modified:
        else:
            print(f'"{file_name}" is already modified.')
        print()
    return files_list
## Analyse the refractive index options:
def ask_material():
    # Refractive index dispersions:
    n_PMMA = lambda wv: (1 + 1.1819 * wv**2 / (wv**2 - 0.011313))**0.5 
    n_silica = lambda wv: (1 + (0.6961663 * wv**2 / (wv**2 - 0.0684043**2)) + (0.4079426 * wv**2 / (wv**2 - 0.1162414**2)) + (0.8974794 * wv**2 / (wv**2 - 9.896161**2)))**0.5
    # Ask the core material:
    print('Available core materials:')
    print(f'n_air = 1.000')
    print(f'n_water = 1.333')
    print()
    while True:
        question = input('Choose the material for the capillary core: [air/water] ')
        try:
            if question not in ['air', 'water']:
                raise ValueError('Unavailable material!')
            else:
                material_core = question
        except ValueError as ve:
            print(ve)
        else:
            break
    # Compare the dispersions for the cladding materials:
    print('You can choose one of the following dispersion relations:')
    wavelength = linspace(0.450, 1.750, 1000) # Supercontinuoum laser spectrum range.
    fig, ax = subplots(nrows=1, ncols=1, layout='constrained', figsize=(12,4))
    ax.set_xlabel('Wavelength [nm]', loc='center', fontsize=14)
    ax.set_ylabel('Refractive index', loc='center', fontsize=14)
    ax.plot(wavelength * 1e3, n_PMMA(wavelength), c='blue', lw=2, linestyle='solid', label='PMMA')
    ax.plot(wavelength * 1e3, n_silica(wavelength), c='green', lw=2, linestyle='solid', label='Silica')
    ax.legend(loc='best', fontsize=12)
    show()
    if save:
        fig.savefig(join(sample_folder, f'Graphic - Comparison Refractive Index.png'))
    print()
    # Ask the cladding material:
    while True:
        question = input('Choose the material for the capillary cladding: [PMMA/silica] ')
        try:
            if question not in ['PMMA', 'silica']:
                raise ValueError('Unavailable material!')
            else:
                material_clad = question
        except ValueError as ve:
            print(ve)
        else:
            break
    # Build the refractive index values:
    if material_clad == 'PMMA':
        color = 'blue'
        label = material_clad
        n_mat = n_PMMA(wavelength)
    else:
        color = 'green'
        label = material_clad.capitalize()
        n_mat = n_silica(wavelength)
    print(f'You chose {label}.')
    # Plot of the chosen material:
    fig, ax = subplots(nrows=1, ncols=1, layout='constrained', figsize=(12,4))
    ax.set_xlabel('Wavelength [nm]', loc='center', fontsize=14)
    ax.set_ylabel('Refractive index', loc='center', fontsize=14)
    ax.plot(wavelength * 1e3, n_mat, c=color, lw=2, linestyle='solid')
    show()
    if save:
        fig.savefig(join(sample_folder, f'Graphic - {label} Refractive Index.png'))
    return material_core, material_clad
## Ask if the user wants to stop the program: 
def ask_stop():
    while True:
        # Aks about the stop:
        question = input('Do you want to stop the program? [y/n]')
        try:
            if question not in ['y', 'n']: # Invalid answer.
                raise ValueError('Invalid answer!')
            elif question == 'y':  # Will stop.
                stop = True
            else: # Won't stop.
                stop = False
        except ValueError as ve: # It occured an error.
            print(ve)
        else: # It ran smoothly.
            break
    return stop
## Main control:
def main():
    # Welcome the user to the program:
    intro()
    # Main loop:
    stop = False
    while not stop:
        # Get the folder containg the samples, the sample to be analysed and its folder:
        get_sample_folder()
        # Ask if the user wants to save the graphics and data frames:
        ask_save()
        # Modify the spectrum data and get their path location:
        file_list = modify_files()
        # Ask the material of the claddings:
        material_core, material_clad = ask_material()
        # Create the capillaries:
        for file in file_list:
            capillary = Capillary(file=file, material_core=material_core, material_clad=material_clad)
            capillary.thickness_variation()
            capillary.get_parameters()
            capillary.calibrate_valley()
            capillary.analyse_data()
            capillary.thickness_from_order()
            capillary.simulate()
            capillary.show_info()
        # Ask if the user wants to stop the program:
        stop = ask_stop()
# Run the program:
main()