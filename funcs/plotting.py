import plotly.graph_objects as go
import numpy as np
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy.integrate import quad

# def nasa9(input_array):
#     data = input_array.astype(float)
    
#     # Generate temp range
#     x_data = np.arange(100, 1501, 100).astype(float)
#     # Slice the last 15 rows for fitting
#     y_data = data[-15:].astype(float)   # Assuming data is a 17x1 array

#     # Perform the curve fitting
#     popt, pcov = curve_fit(cp, x_data, y_data)
    
#     a0, a1, a2, a3, a4, a5, a6 = popt
    
#     a7 = h0(data[0], a0, a1, a2, a3, a4, a5, a6)
    
#     a8 = s0(data[1], a0, a1, a2, a3, a4, a5, a6)
    
#     return(popt, a7, a8)
    

class nasa9fit:
    def __init__(self, input_array):
        self.input_array = input_array
        self.href = input_array[0]
        self.sref = input_array[1]
        self.popt = self.fit_nasa9()

    def fit_nasa9(self):
        ### Define the polynomial function for cp, unit cal/K/mol.
        ### Cp/R
        def cp(x, a0, a1, a2, a3, a4, a5, a6):
            return (a0 * x**(-2) + a1 * x**(-1) + a2 + a3 * x + a4 * x**2 + a5 * x**3 + a6 * x**4) * 1.987204 
    
        data = self.input_array.astype(float)
        x_data = np.arange(100, 1501, 100).astype(float)
        y_data = data[-15:].astype(float)
        popt, pcov = curve_fit(cp, x_data, y_data)
        a0, a1, a2, a3, a4, a5, a6 = popt
               
        return popt
 
    ### Define the polynomial function for h0, unit kcal/mol
    ### eval temp = 298.15
    def h0(x, href, a0, a1, a2, a3, a4, a5, a6):
        a7 = []     
        a7 = (href/(x*1.987204258*1000) - (-a0 * x**(-2) + a1 * np.log(x)/x + a2 + a3/2 * x + a4/3 * x**2 + a5/4 * x**3 + a6/5 * x**4) ) *x
        return a7

    ### Cp evaluator at diff temperatures        
    def C(self, x):    ### Sanity check for H, passed
        a0, a1, a2, a3, a4, a5, a6 = self.popt
        return (a0 * x**(-2) + a1 * x**(-1) + a2 + a3 * x + a4 * x**2 + a5 * x**3 + a6 * x**4) * 1.987204 


    # ### Define the polynomial function for s0, unit kcal/mol
    # ### eval temp = 298.15
    # def s0(x, sref, a0, a1, a2, a3, a4, a5, a6):
    #     a8=[]
    #     a8 = sref/1.987204 - (-a0/2 * x**(-2) - a1 * x**(-1) + a2*np.log(x) + a3 * x + a4/2 * x**2 + a5/3 * x**3 + a6/4 * x**4)
    #     return a8
       
    ### H evaluator at diff temperatures
    def H(self, x):
        T_min = 298.15
        T_max = x
        # Perform numerical integration to evaluate H(T)
        H_T, _ = quad(self.C, T_min, T_max)
        # Add the reference enthalpy
        H = self.href + H_T/1000.
        return H
    
    ### S evaluator at diff temperatures
    def S(self, x):
        T_min = 298.15
        T_max = x
        # Define the integrand function for entropy
        def integrand(T):
            return ( self.C(T) / T  )
        
        # Perform numerical integration to evaluate ΔS
        delta_S, _ = quad(integrand, T_min, T_max)

        # Calculate entropy at temperature T
        S_T = self.sref + delta_S
        return S_T  
    
    # def H(self, x):    ### Sanity check for H, passed
    #     a0, a1, a2, a3, a4, a5, a6 = self.popt
    #     a7 = self.h0(x, self.href, a0, a1, a2, a3, a4, a5, a6)
    #     return( ((-a0 * x**(-2) + a1 * np.log(x)/x + a2 + a3/2 * x + a4/3 * x**2 + a5/4 * x**3 + a6/5 * x**4) + a7/x) * (x*1.987204258*1000) )

    # ### H evaluator at diff temperatures        
    # def S(self, x):    ### Sanity check for H, passed
    #     a0, a1, a2, a3, a4, a5, a6 = self.popt
    #     a8 = self.s0(x, self.sref, a0, a1, a2, a3, a4, a5, a6)
    #     return(   (-a0/2 * x**(-2) - a1 * x**(-1) + a2*np.log(x) + a3 * x + a4/2 * x**2 + a5/3 * x**3 + a6/4 * x**4 + a8) * 1.987204  )




# # Create a scatter plot for the discrete data points
    ### Cp plotter across T[100,1500] K
    def plotC(self):
        # Generate temperature points
        temperature_range = np.linspace(100, 1500, 100)
        # Calculate heat capacity at each temperature point
        capacity_values = []
        for T in temperature_range:
            capacity_values.append( self.C(T) )
        # Plot the results using Plotly
        fig = go.Figure()
        
        # Add enthalpy plot
        fig.add_trace(go.Scatter(x=temperature_range, y=capacity_values, mode='lines', name='Heat Capacity Cp'))
        
        
        x_values2 = list(range(100, 1501, 100))
        y_values2 = self.input_array.astype(float)[-15:]
        
        # Add second series
        fig.add_scatter(x=x_values2, y=y_values2, mode='markers', name='METE with 1st contributions and 2nd corrections')

        # Update layout
        fig.update_layout(title='Heat Capacity Cp vs Temperature',
                          title_x=0.5,
                          xaxis_title='Temperature (K)',
                          yaxis_title='Value',
                          height=800,  # Set canvas height in pixels
                          width=1000,    # Set canvas width in pixels
                          # Update the layout to move the legend inside the canvas to the bottom right
                          legend=dict(x=1, y=0.08, xanchor='right', yanchor='bottom'))

        # Show the plot
        fig.show()

    ### H plotter across T[298.15,1500] K
    def plotH(self):
        # Generate temperature points
        temperature_range = np.linspace(298.15, 1500, 100)
        # Calculate enthalpy at each temperature point
        enthalpy_values = []
        for T in temperature_range:
            enthalpy_values.append( self.H(T) )
        # Plot the results using Plotly
        fig = go.Figure()
        
        # Add enthalpy plot
        fig.add_trace(go.Scatter(x=temperature_range, y=enthalpy_values, mode='lines', name='Enthalpy'))
        # Update layout
        fig.update_layout(title='Enthalpy vs Temperature',
                          title_x=0.5,
                          xaxis_title='Temperature (K)',
                          yaxis_title='Value',
                          height=800,  # Set canvas height in pixels
                          width=1000    # Set canvas width in pixels
                          )

        # Show the plot
        fig.show()
        
    ### S plotter across T[298.15,1500] K    
    def plotS(self):
        
        # Create a new figure
        fig = go.Figure()
        # Generate temperature points
        temperature_range = np.linspace(298.15, 1500, 100)
        # Calculate entropy at each temperature point
        entropy_values = []
        for T in temperature_range:
            entropy_values.append( self.S(T) )
        # Plot the results using Plotly
        fig = go.Figure()
        
        # Add enthalpy plot
        fig.add_trace(go.Scatter(x=temperature_range, y=entropy_values, mode='lines', name='Entropy'))
        # Update layout
        fig.update_layout(title='Entropy vs Temperature',
                          title_x=0.5,
                          xaxis_title='Temperature (K)',
                          yaxis_title='Value',
                          height=800,  # Set canvas height in pixels
                          width=1000    # Set canvas width in pixels
                          )

        # Show the plot
        fig.show()






