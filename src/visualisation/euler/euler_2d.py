import itertools
import json
import os

import click
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DEFAULTS = {
    'results_dir': '.',  # Current working directory
    'output_dir': '.',  # Current working directory
    'flow_var': 'rho',  # Density
    'colour_map': 'viridis',  # Matplotlib colour map
    'figure_size': '8,8',  # Figure size in inches, assumes a 100dpi
    'limits': None,  # Range of the flow variable to assign to colour map extent
    'cells': None,  # Total number of cells in x and y direction for the domain extent
    'contours': False,  # Whether to utilise a filled contour plot instead of direct flow values
    'num_contours': 30  # How many contour levels to utilise if contours are selected
}

RESULTS_FILE_HEADER_ROWS = 8  # Number of rows of meta data in results files


def extract_figsize(figsize):
    """
    Extracts the figure size in inches (at 100 dpi) from the
    comma-separated string of floats.

    Parameters
    ----------
    figsize : `str`
        The comma-separated string of floats, e.g. '8,8'

    Returns
    -------
    `tuple(float)`
        The extracted two-tuple of floats, e.g. (8.0, 8.0)
    """
    try:
        figsize_tup = tuple(map(float, figsize.split(',')))
    except:
        raise ValueError(
            "Unable to extract appropriate figure size via '%s'." % figsize
        )
    else:
        return figsize_tup

def extract_limits(limits):
    """
    Extracts the range of the flow variable to clip the plotting
    to from the comma-separated string of floats.

    Parameters
    ----------
    figsize : `str`
        The comma-separated string of floats, e.g. '0.0,5.0'

    Returns
    -------
    `tuple(float)`
        The extracted two-tuple of floats, e.g. (0.0, 5.0)
    """    
    if limits is None:
        return None
    else:
        try:
            limits_tup = tuple(map(float, limits.split(',')))
        except:
            raise ValueError(
                "Unable to extract limits from limit string: '%s'" % limits
            )
        else:
            return limits_tup

def extract_cells(cells):
    """
    Extracts the number of cells in each of the x and y directions
    for the total domain extent from the comma-separated
    string of ints.

    This is used to plot the entire encompassing area of the domain,
    even if it is comprised of multiple individual blocks.

    Parameters
    ----------
    figsize : `str`
        The comma-separated string of ints, e.g. '100,100'

    Returns
    -------
    `tuple(int)`
        The extracted two-tuple of ints, e.g. (100, 100)
    """    
    try:
        cells_tup = tuple(map(int, cells.split(',')))
    except:
        raise ValueError(
            "Unable to extract cells from cell string: '%s'" % cells
        )
    else:
        return cells_tup

def load_sim_json_dict(sim_json):
    """
    Loads the simulation JSON file and parses it.

    Parameters
    ----------
    sim_json : `str`
        The filepath to the simulation JSON file.

    Returns
    -------
    `dict`
        The dictionary representing the simulation specification.
    """
    if sim_json is None:
        raise ValueError(
            'Simulation JSON file has not been specified.'
        )
    return json.load(open(sim_json, 'r'))

def load_results_names(results_dir):
    """
    Loads all of the results files in a particular directory
    and sorts them.

    Parameters
    ----------
    results_dir : `str`
        The results directory string.

    Returns
    -------
    `list[str]`
        The list of all the results files in the results directory.
    """
    return sorted(
        [
            res for res in os.listdir(results_dir)
            if res.startswith('results_')
        ]
    )


class EulerBlockStructuredVisualisation2D(object):
    """
    Provides two-dimensional visualisation plots of compressible
    Euler flow variables or their 'schlieren'  gradients derived
    from HyperFlow's native results output files.

    This can be used to generate all frames of an animation for an
    unsteady simulation (or an animation that plots convergence
    to steady state).

    Parameters
    ----------
    sim_json : `dict`
        Parsed JSON simulation parameters dictionary
    results_dir : `str`
        Directory to find the HyperFlow output results files
    output_dir : `str`
        Directory to store the animation frames of the visualisation
    flow_var : `str`
        The particular flow variable to utilise for plotting
    colour_map : `str`
        The Matplotlib colour map to use for plotting
    figsize : `tuple(float)`
        A tuple of the figure size in inches at 100dpi
    limits : `tuple(float)`
        A tuple (e.g. (0.0, 5.0)) of the range limits of the flow variable
        used to scale the Matplotlib colour map
    cells : `tuple(int)`
        The number of cells for the domain in the x and y direction
    contours : `boolean`
        Whether to display a contour plot
    num_contours : `int`
        Number of contours to use if a contour plot is chosen
    """

    # Parameters used in the 'schlieren-like' normalised gradient
    SCHLIEREN_BETA = 1.0
    SCHLIEREN_LAMBDA = 15.0    

    def __init__(
        self,
        sim_json,
        results_dir=DEFAULTS['results_dir'],
        output_dir=DEFAULTS['output_dir'],
        flow_var=DEFAULTS['flow_var'],
        colour_map=DEFAULTS['colour_map'],
        figsize=DEFAULTS['figure_size'],
        limits=DEFAULTS['limits'],
        cells=DEFAULTS['cells'],
        contours=DEFAULTS['contours'],
        num_contours=DEFAULTS['num_contours']
    ):
        self.sim_json = sim_json
        self.results_dir = results_dir
        self.output_dir = output_dir
        self.flow_var = flow_var
        self.colour_map = colour_map
        self.figsize = figsize
        self.limits = limits
        self.cells = cells
        self.contours = contours
        self.num_contours = num_contours

        self.extents = self._load_extents_from_json()

    def _normalised_gradient_exponential_schlieren(
        self,
        flow_data
    ):
        """
        Calculate a 'schlieren-like' normalised exponentiated gradient
        for a particular flow variable.

        This is particularly useful for highlighting shocks/discontinuities
        in the flow.

        Parameters
        ----------
        flow_data : `numpy.array`
            Multidimensional array of flow variable data by x, y coordinates
        
        Returns
        -------
        `numpy.array`
            Multidimensional array of the 'schlieren-like' normalised
            exponentiated gradient.
        """

        # Calculate total domain width, height
        # NOTE: Assumes homogeneous dx, dy
        width = self.extents[1] - self.extents[0]
        height = self.extents[3] - self.extents[2]
        
        # Calculate domain step size in each direction
        dx = width / self.cells[0] 
        dy = height / self.cells[1]

        # Pre-allocate normalised gradient schlieren array
        norm_grad = np.zeros(flow_data.shape)
        max_norm_grad = -1e6
        min_norm_grad = 1e6

        # Iterate over all cells 
        for i in range(1, flow_data.shape[0] - 1):
            for j in range(1, flow_data.shape[1]- 1):
                divx = 0.0
                divy = 0.0

                # Depending upon whether at the boundary or an internal
                # cell utilise a one-sided or central difference respectively
                # in the x-direction
                if (i == 0):
                    divx = (flow_data[i+1][j] - flow_data[i][j]) / dx
                elif (i == flow_data.shape[0] - 1):
                    divx = (flow_data[i][j] - flow_data[i-1][j]) / dx
                else:
                    divx = (flow_data[i+1][j] - flow_data[i-1][j]) / (2.0 * dx)

                # Depending upon whether at the boundary or an internal
                # cell utilise a one-sided or central difference respectively
                # in the y-direction
                if (j == 0):
                    divy = (flow_data[i][j+1] - flow_data[i][j]) / dy
                elif (j == flow_data.shape[1] - 1):
                    divy = (flow_data[i][j] - flow_data[i][j-1]) / dy
                else:
                    divy = (flow_data[i][j+1] - flow_data[i][j-1]) / (2.0 * dy)

                # Calculate the normalised gradient and store the
                # current maximum normalised gradient
                try:
                    norm_grad[i][j] = np.sqrt(divx**2 + divy**2)
                    if (norm_grad[i][j] > max_norm_grad):
                        max_norm_grad = norm_grad[i][j]
                    if (norm_grad[i][j] < min_norm_grad):
                        min_norm_grad = norm_grad[i][j]
                except:
                    pass

        # Calculate the 'schlieren-like' exponentiated normalised gradient for
        for i in range(0, flow_data.shape[0]):
            for j in range(0, flow_data.shape[1]):
                norm_grad[i][j] = self.SCHLIEREN_BETA * np.exp(
                    -self.SCHLIEREN_LAMBDA * norm_grad[i][j] / max_norm_grad
                )

        return norm_grad

    def _load_extents_from_json(
        self
    ):
        """
        Parses the JSON simulation file for the domain extents
        and adds them to a list of floats for convenience.

        Parameters
        ----------
        None

        Returns
        -------
        `list[float]`
            The list of coordinate extents of the enclosing domain.

        """
        domain = self.sim_json['domain']
        return [
            domain['x_left'],
            domain['x_right'],
            domain['y_bottom'],
            domain['y_top']
        ]

    def _set_axis_no_gaps(
        self,
        axis
    ):
        """
        Turns off the plot axis so that the imshow call covers the
        entire extent of the figure area.

        Parameters
        ----------
        axis : `matplotlib.pyplot.axis`
            The axis entity used to control the plot.

        Returns
        -------
        None
        """
        axis.xaxis.set_major_locator(plt.NullLocator())
        axis.yaxis.set_major_locator(plt.NullLocator())
        axis.set_axis_off()
        axis.margins(0,0)

    def _plot_frame(
        self,
        frame_filepath,
        flow_data,
        x_coords,
        y_coords
    ):
        """
        Creates a figure and axis, then plots the data via a contour plot
        or imshow depending upon the chosen plot type.

        Parameters
        ----------
        frame_filepath : `str`
            The full path to the output frame
        flow_data : `numpy.array`
            Multidimensional array of flow data
        x_coords : `numpy.array`
            Single-dimensional array of all cell-centre x-coordinate locations
        y_coords : `numpy.array`
            Single-dimensional array of all cell-centre y-coordinate locations

        Returns
        -------
        None
        """
        fig, axis = plt.subplots(1, 1, figsize=(self.figsize[0], self.figsize[1]))
        
        fig.subplots_adjust(left=0, bottom=0, right=1, top=1, hspace=0, wspace=0)
        self._set_axis_no_gaps(axis)

        if self.contours:
            contour_levels = np.linspace(np.nanmin(flow_data), np.nanmax(flow_data), self.num_contours)
            if all(x == contour_levels[0] for x in contour_levels):
                contour_levels = np.array([contour_levels[0], contour_levels[0] + 1e-12])

            cs = axis.contourf(
                x_coords, y_coords, 
                np.swapaxes(flow_data, 0, 1),
                contour_levels, cmap=self.colour_map, vmin=self.limits[0], vmax=self.limits[1]
            )
            axis.contour(cs, colors='black', linewidths=0.5)
        else:
            flow_data = np.swapaxes(flow_data, 0, 1)
            axis.imshow(
                flow_data,
                origin='lower',
                cmap=self.colour_map,
                vmin=self.limits[0],
                vmax=self.limits[1],
                extent=self.extents
            )

        plt.savefig(frame_filepath, dpi=100)
        plt.close()

    def plot_figures(
        self,
        results
    ):
        """
        Iterates through all of the result file names and generates
        a single frame for each file.

        The flow variable can itself be plotted or its schlieren/gradient
        with either a full 'imshow' or as a filled contour plot.

        Parameters
        ----------
        """
        for frame_idx, result_file in enumerate(results):
            print("Plotting Frame %06d..." % frame_idx)
            
            # Obtain the full file paths to the frame output and results data
            frame_filepath = os.path.join(self.output_dir, 'frame_%06d.png' % frame_idx)
            sim_filepath = os.path.join(self.results_dir, result_file)

            # Skip all of the meta data when loading the results file
            sim_data = pd.read_csv(sim_filepath, skiprows=RESULTS_FILE_HEADER_ROWS).dropna()
            
            x = sim_data['x'].to_numpy()
            y = sim_data['y'].to_numpy()
            sim_data = sim_data.set_index(['x', 'y'])      

            # Obtain arrays of the unique cell-centre x and y coordinates
            x_coords = np.array(sim_data.index.levels[0].unique())
            y_coords = np.array(sim_data.index.levels[1].unique())
            sim_data = sim_data.reindex(itertools.product(x_coords, y_coords))

            # Reshape the data frame into the appropriate format for plotting
            flow_data = sim_data[
                self.flow_var.replace('-gradient', '')
            ].values.reshape(self.cells[0], self.cells[1])

            # Plot the schlieren-like gradient plot
            if '-gradient' in self.flow_var:
                flow_data = self._normalised_gradient_exponential_schlieren(
                    flow_data, self.cells[0], self.cells[1], self.extents
                )
            
            self._plot_frame(
                frame_filepath,
                flow_data,
                x_coords,
                y_coords
            )


@click.command()
@click.option('--sim-json', 'sim_json', default=None, help='Simulation JSON file full location')
@click.option('--results-dir', 'results_dir', default=DEFAULTS['results_dir'], help='System location from which to read results files')
@click.option('--output-dir', 'output_dir', default=DEFAULTS['output_dir'], help='System location to store frame output plots')
@click.option('--flow-var', 'flow_var', default=DEFAULTS['flow_var'], help='Flow variables to plot from "rho", "u", "v", "p" or "rho-gradient", "p-gradient" etc.')
@click.option('--colour-map', 'colour_map', default=DEFAULTS['colour_map'], help='Matplotlib colour map for plot')
@click.option('--figsize', 'figsize', default=DEFAULTS['figure_size'], help='Plot figure size in inches')
@click.option('--limits', 'limits', default=DEFAULTS['limits'], help='Flow variable plot limits, e.g. 0.0,10.0')
@click.option('--cells', 'cells', default=DEFAULTS['cells'], help='Total cells in x and y direction of domain, e.g. 100,100')
@click.option('--contours', 'contours', is_flag=True, default=DEFAULTS['contours'], help='Whether to plot contours instead of imshow')
@click.option('--num-contours', 'num_contours', default=DEFAULTS['num_contours'], help='Whether to plot contours instead of imshow')
def cli(sim_json, results_dir, output_dir, flow_var, colour_map, figsize, limits, cells, contours, num_contours):
    sim_json = load_sim_json_dict(sim_json)
    figsize = extract_figsize(figsize)
    limits = extract_limits(limits)
    cells = extract_cells(cells)
    results = load_results_names(results_dir)

    euler_viz_2d = EulerBlockStructuredVisualisation2D(
        sim_json,
        results_dir=results_dir,
        output_dir=output_dir,
        flow_var=flow_var,
        colour_map=colour_map,
        figsize=figsize,
        limits=limits,
        cells=cells,
        contours=contours
    )
    euler_viz_2d.plot_figures(results)


if __name__ == "__main__":
    cli()