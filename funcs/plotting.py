import plotly.graph_objects as go
import numpy as np

# Define the function- NASA 9
def custom_function(x, arr):
    a0 = arr[0]
    a1 = arr[1]
    a2 = arr[2]
    a3 = arr[3]
    a4 = arr[4]
    a5 = arr[5]
    a6 = arr[6]
    return a0 * x**(-2) + a1 * x**(-1) + a2 + a3 * x + a4 * x**2 + a5 * x**3 + a6 * x**4

# Generate a range of x values
x_values = np.linspace(100, 1500, 100)
# Calculate the corresponding y values using the function

arr1 = [(108.27573020999999, 25.850120253, -173325.23306, 3861.8333912999997, -25.272258364, 0.16347287100000002, -0.0001603, 8.07e-08, -1.628e-11)]
arr1 = arr1[0][2:]

y_values = custom_function(x_values, arr1)

# Create a Plotly figure
fig = go.Figure()

###GA1#####
# Add the function plot
fig.add_trace(go.Scatter(x=x_values, y=y_values, mode='lines', name='NASA9 fitted'))

# Define scatter points
x_scatter = np.linspace(100, 1500, 15)
arr2 = [([-72.338    ,  28.1785   ,  11.4558   ,  17.438    ,  22.9934   ,  28.8075   ,  34.0245   ,  38.4188   ,  42.1035   ,  45.2248   ,  47.8965   ,
        50.1967   ,  52.1827   ,  53.9024   ,  55.3935   ,  56.6892   ,  57.8183224])]
y_scatter = arr2[0][2:]
y_scatter = custom_function(x_scatter, arr1)
# Add scatter points
fig.add_trace(go.Scatter(x=x_scatter, y=y_scatter, mode='markers', name='DFT Data'))

###GA2#####
#################################################################################################
x_scatter = np.linspace(100, 1500, 15)
arr2 = ([108.25, 25.98 ,10.89 ,16.7 ,22.36 ,27.91, 32.65 ,36.49, 39.62, 42.2, 44.36, 46.19, 47.73 ,49.05 ,50.18, 51.15, 51.99])
y_scatter = arr2[2:]
# Add scatter points
fig.add_trace(go.Scatter(x=x_scatter, y=y_scatter, mode='markers', name='DFT Data w/ GA2 correction'))


# Update layout with custom canvas size
fig.update_layout(
    title="Plot of Cps between [100,1500] K",
    xaxis_title="x",
    yaxis_title="y",
    template="plotly_white",
    height=1000,  # Set canvas height in pixels
    width=800    # Set canvas width in pixels
)

# Show the plot
fig.show()