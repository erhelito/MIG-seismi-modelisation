import os
import pickle
import matplotlib.pyplot as plt
import numpy as np

class Result1D:
    "Results storage"

    def __init__(self, x, T, V, Nu, Phi, pnd, pc, filename = ''):
        self.x=x
        self.T=T
        self.V=V
        self.Nu=Nu
        self.Phi=Phi
        self.pnd=pnd
        self.pc=pc

        if filename=='' or filename is None:
            # Build a robust default filename using available scalar attributes.
            try:
                # pnd.a may be an array: take its mean as representative
                a_val = float(np.mean(pnd.a)) if hasattr(pnd, 'a') else 0.0
            except Exception:
                a_val = 0.0

            eta_val = float(getattr(pnd, 'eta', 0.0))
            h_val = float(getattr(pnd, 'h', 0.0))

            self.filename = f"a{a_val:.2e}_eta{eta_val:.2e}_h{h_val:.2e}.pkl"
        else:
            self.filename = filename

    def save_results(self, folder_name=''):
        # Save results under the package-relative `model2/Results1D/<folder_name>/`
        base_dir = os.path.dirname(__file__)
        results_dir = os.path.join(base_dir, 'Results1D')

        if folder_name:
            results_dir = os.path.join(results_dir, folder_name)

        os.makedirs(results_dir, exist_ok=True)

        outpath = os.path.join(results_dir, self.filename)
        
        with open(outpath, 'wb') as f:
            pickle.dump(self, f)

    @staticmethod
    def load_results(path=''):
        """
        Load previously saved results from a pickle file.

        Parameters
        ----------
        filename : str
            The name of the pickle file to load.

        Returns
        -------
        data : Result
            The loaded results.
        """
        
        import oneD_CFAULT_RNS as prmod

        with open(path, 'rb') as f:
            data = pickle.load(f)

        return data

    def theta_evolution(self, save = False, path = '', name=''):
        """
        Plot the theta evolution.

        Parameters
        ----------
        save : bool
            If True, save the figure in the given path.
        path : str
            The path where the figure will be saved.
        name : str
            The name of the figure to save.
        """
        plt.figure('Theta evolution')

        plt.plot(self.T,self.Nu,'-k')
        plt.xlabel('Time (ND)')
        plt.ylabel('Log of state variable (ND)')
        plt.title(r'$\theta$ evolution')
        plt.grid()

        if save:
            if name != '':
                plt.savefig(f"{path}/theta_evolution_{name}.pdf")
            else:
                plt.savefig(f"{path}/theta_evolution_{self.filename}.pdf")

    def slip_rate_evolution(self, save = False, path = '', name=''):
        """
        Plot the slip rate evolution.

        Parameters
        ----------
        save : bool
            If True, save the figure in the given path.
        path : str
            The path where the figure will be saved.
        name : str
            The name of the figure to save.
        """
        plt.figure('Slip rate evolution')

        plt.contour(self.x,self.T,np.log10(self.V),30)

        plt.xlabel('Position along the fault (ND)')
        plt.ylabel('Time (ND)')
        plt.grid()
        plt.colorbar(label='Log10 Slip rate (ND)')

        if save:
            if name != '':
                plt.savefig(f"{path}/slip_rate_evolution_{name}.pdf")
            else:
                plt.savefig(f"{path}/slip_rate_evolution_{self.filename}.pdf")

    def slip_rate_evolution_by_iterations(self, save = False, path = '', name=''):
        """
        Plot the slip rate evolution, by iterations.

        Parameters
        ----------
        save : bool
            If True, save the figure in the given path.
        path : str
            The path where the figure will be saved.
        name : str
            The name of the figure to save.
        """
        plt.figure('Slip rate evolution by iterations')

        plt.contour(self.x,np.arange(len(self.T)),np.log10(self.V),30)

        plt.xlabel('Position along the fault (ND)')
        plt.ylabel('Iteration (ND)')
        plt.grid()
        plt.colorbar(label='Log10 Slip rate (ND)')

        if save:
            if name != '':
                plt.savefig(f"{path}/slip_rate_evolution_{name}.pdf")
            else:
                plt.savefig(f"{path}/slip_rate_evolution_{self.filename}.pdf")

    def slip_rate_evolution_3D(self, save=False, path='', name=''):
        """
        Plot the slip rate evolution in 3D.

        Parameters
        ----------
        save : bool
            If True, save the figure in the given path.
        path : str
            The path where the figure will be saved.
        name : str
            The name of the figure to save.
        """
        
        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

        fig = plt.figure('Slip rate 3D evolution')
        ax = fig.add_subplot(projection='3d')

        # Create meshgrid matching the contour call: X (positions) repeated along rows,
        # Y (times) repeated along columns so shapes match V (T x x).
        X, Y = np.meshgrid(self.x, self.T)
        Z = np.log10(self.V)

        surf = ax.plot_surface(X, Y, Z, cmap='viridis')

        ax.set_title('Slip rate evolution 3D')

        ax.set_xlabel('Position along the fault (ND)')
        ax.set_ylabel('Time (ND)')

        fig.colorbar(surf, label='Log10 Slip rate (ND)')

        if save:
            if name != '':
                plt.savefig(f"{path}/slip_rate_evolution_3D_{name}.pdf")
            else:
                plt.savefig(f"{path}/slip_rate_evolution_3D_{self.filename}.pdf")

    def slip_rate_evolution_map(self, save=False, path='', name=''):
        """
        Plot the slip rate evolution in a colormap.

        Parameters
        ----------
        save : bool
            If True, save the figure in the given path.
        path : str
            The path where the figure will be saved.
        name : str
            The name of the figure to save.
        """
        plt.figure("Slip rate colormap evolution")

        plt.pcolormesh(self.x, self.T, np.log10(self.V), cmap='viridis', shading='auto')
        plt.xlabel('Position along the fault (ND)')
        plt.ylabel('Time (ND)')
        plt.grid()
        plt.colorbar(label='Log10 Slip rate (ND)')

        if save:
            if name != '':
                plt.savefig(f"{path}/slip_rate_evolution_map_{name}.pdf")
            else:
                plt.savefig(f"{path}/slip_rate_evolution_map_{self.filename}.pdf")