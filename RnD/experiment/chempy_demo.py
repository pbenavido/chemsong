import numpy as np
from chempy import Reaction, Substance
from chempy.kinetics.ode import get_odesys
from chempy.kinetics.rates import MassAction

# Define substances
benzene = Substance.from_formula("C6H6")
ethene = Substance.from_formula("C2H4")
ethylbenzene = Substance.from_formula("C8H10")

# Define the reaction (simplified as a one-step reaction for demonstration)
reaction = Reaction(
    {"C6H6": 1, "C2H4": 1},  # Reactants
    {"C8H10": 1},  # Products
    MassAction(4.0e10),  # Rate constant (hypothetical)
)

# Define the initial conditions
initial_conditions = {"C6H6": 1.0, "C2H4": 2.0, "C8H10": 0}  # Molar concentrations
time = np.linspace(0, 10, 100)  # Time from 0 to 10 seconds

# Get the ordinary differential equation system
odesys, extra = get_odesys(reaction)

# Integrate the ODE system
result = odesys.integrate(time, initial_conditions)

# Create a result summary
summary = {
    "Time (s)": time,
    "Concentration of Benzene (mol/L)": result.yout[:, extra["species"].index("C6H6")],
    "Concentration of Ethene (mol/L)": result.yout[:, extra["species"].index("C2H4")],
    "Concentration of Ethylbenzene (mol/L)": result.yout[
        :, extra["species"].index("C8H10")
    ],
}

# Print the summary for specific time points (e.g., every 10 seconds)
for i in range(0, len(time), 10):
    print({key: value[i] for key, value in summary.items()})
