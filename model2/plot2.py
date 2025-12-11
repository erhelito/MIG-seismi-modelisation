from result1D import Result1D
from oneD_CFAULT_RNS import NdParamMec, ParamComp
import matplotlib.pyplot as plt

path = 'model2/Results1D/test/asp_test_simulation_2.pkl'

data = Result1D.load_results(path)
data.theta_evolution() # Plot state variable evolution
data.slip_rate_evolution() # Plot slip rate evolution
data.slip_rate_evolution_by_iterations()
data.slip_rate_evolution_3D()
data.slip_rate_evolution_map()

plt.show()