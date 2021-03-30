from enum import Enum
import json

import click
import numpy as np


class BoundaryCondition(Enum):
    """
    Enumerates all of the possible boundary condition
    indices utilised by the mesh importer.
    """
    INTERNAL = 0
    SOLID = -1
    TRANSMISSIVE = -2
    INLET = -3


class BlockStructuredMeshGenerator:
    """
    Generates the cell coordinates for a block structured
    mesh. Also generates the list of integers corresponding
    to each edge that specify the particular cell's boundary
    condition.

    The mesh importer requires a JSON specification file that
    provides a list of all box bounding coordinates, cells
    in the x,y directions and boundary conditions for each
    of the four sides of the quadrilateral block.

    All block coordinates are specified in millimetres and
    converted internally by this class to metres for output
    to HyperFlow.

    Parameters
    ----------
    block_json_file : `string`
        Full path to the JSON file specifying the blocks.
    """

    def __init__(self, block_json_file):
        self.block_json_file = block_json_file

    def _parse_block_json_file(self):
        """
        Parse the provided JSON file into a dictionary of block records.

        Returns
        -------
        `dict`
            The individual block specification records.
        """
        return json.load(open(self.block_json_file, 'r'))

    def _xy_scale(self, block, eps, nu):
        """
        Map the cell in the unit square (itself with-1/+1 corners)
        to the cell in quadrilateral physical (x,y) space.

        Parameters
        ----------
        block : `dict`
            Dictionary record for an individual block
        eps : `float`
            Horizontal coordinate in grid-aligned coordinate system
        nu : `float`
            Vertical coordinate in grid-aligned coordinate system

        Returns
        -------
        `tuple(float)`
            Individual x, y coordinates in physical space for a particular
            block coordinate.
        """
        N = np.array([
            0.25 * (1.0 - eps) * (1.0 - nu),
            0.25 * (1.0 + eps) * (1.0 - nu),
            0.25 * (1.0 + eps) * (1.0 + nu),
            0.25 * (1.0 - eps) * (1.0 + nu)
        ])

        # Obtain separate arrays of all x and y physical coordinates
        x = np.array([coord_pair[0] for coord_pair in block['coordinates']])
        y = np.array([coord_pair[1] for coord_pair in block['coordinates']])

        # Scale the unit square (with corners -1/+1) to the quadrilateral
        return N.dot(x), N.dot(y)

    def _calculate_cell_vertices(self, block, eps, nu, d_eps, d_nu):
        """
        Creates the bottom left, bottom right, top right and top left
        coordinate values in the physical (x, y) coordinate system
        from the grid-aligned coordinate system.

        Parameters
        ----------
        block : `dict`
            Dictionary record for an individual block
        eps : `float`
            Horizontal coordinate in grid-aligned coordinate system
        nu : `float`
            Vertical coordinate in grid-aligned coordinate system
        d_eps : `float`
            Increment in horizontal coordinate in grid-aligned coordinate system
        d_nu : `float`
            Increment in vertical coordinate in grid-aligned coordinate system

        Returns
        -------
        `list(tuple(floats))`
            A list comprising the four corner coordinates in physical space
            of each of the blocks. [Bottom left, bottom right, top right, top left]
        """
        bl = self._xy_scale(block, eps, nu)
        br = self._xy_scale(block, eps + d_eps, nu)
        tr = self._xy_scale(block, eps + d_eps, nu + d_nu)
        tl = self._xy_scale(block, eps, nu + d_nu)
        return [bl, br, tr, tl]

    def _calculate_edge_boundary_indices(
        self, boundary_conditions, eps_cells, nu_cells, cell_eps_idx, cell_nu_idx
    ):
        """
        Calculate the appropriate boundary conditions for a particular
        cell within the block.

        Parameters
        ----------
        boundary_conditions: `list[integer]`
            The list of south, east, north and west boundary condition integers
        eps_cells : `int`
            Total number of cells in the epsilon direction
        nu_cells : `int`
            Total number of cells in the nu direction
        cell_eps_idx : `int`
            Index of the current epsilon cell
        cell_nu_idx : `int`
            Index of the current nu cell

        Returns
        -------
        `list[integer]`
            The appropriate boundary conditions for a particular cell.
        """
        edge_boundary_indices = [0] * 4 

        # South
        if (cell_nu_idx == 0):
            edge_boundary_indices[0] = boundary_conditions[0]

        # East
        if (cell_eps_idx == (eps_cells - 1)):
            edge_boundary_indices[1] = boundary_conditions[1]

        # North
        if (cell_nu_idx == (nu_cells - 1)):
            edge_boundary_indices[2] = boundary_conditions[2]

        # West
        if (cell_eps_idx == 0):
            edge_boundary_indices[3] = boundary_conditions[3]

        return edge_boundary_indices

    def _discretise_domain_into_cells(self, blocks):
        """
        Convert all blocks into individual cells with specified
        physical (x, y) coordinates and integers for boundary
        conditions.

        Parameters
        ----------
        blocks : `dict`
            The JSON-parsed block specifications.

        Returns
        -------
        `tuple(list, list)`
            The cell vertices and corresponding edge
            boundary conditions integers list.
        """
        cell_vertices_list = []
        edge_boundary_list = []

        for block in blocks:
            eps_cells = block['cells'][0]
            nu_cells = block['cells'][1]

            # Calculate spatial steps in block space
            # for a unit square with corners -1/+1 in (eps, nu)
            d_eps = 2.0 / eps_cells
            d_nu = 2.0 / nu_cells

            for cell_eps_idx in range(eps_cells):
                for cell_nu_idx in range(nu_cells):
                    # Calculate (eps, nu) coordinates in block-space for
                    # a particular cell
                    eps = cell_eps_idx * d_eps - 1.0
                    nu = cell_nu_idx * d_nu - 1.0
                    
                    # Append the cell (x, y) coordinates in physical space
                    cell_vertices_list.append(
                        self._calculate_cell_vertices(
                            block, eps, nu, d_eps, d_nu
                        )
                    )
                    
                    # Append the cell boundary condition integers
                    edge_boundary_list.append(
                        self._calculate_edge_boundary_indices(
                            block['boundary_conditions'], eps_cells, nu_cells, cell_eps_idx, cell_nu_idx
                        )
                    )

        return cell_vertices_list, edge_boundary_list

    def _output_mesh(self, output_file, cell_vertices_list, edge_boundary_list):
        """
        Output the cell coordinates and edge boundary integers to disk.

        Parameters
        ----------
        output_file : `string`
            The name of the output mesh file.
        cell_vertices_list : `list`
            The list of cell coordinates in physical (x, y) space.
        edge_boundary_list : `list`
            The list of integers representing the cell boundary condition types.
        """
        with open(output_file, 'w') as outfile:
            # Output the cell vertices for each cell
            outfile.write('CELLS\n')
            outfile.write('%d\n' % len(cell_vertices_list))
            for cell_vertices in cell_vertices_list:
                for i, cell_vertex in enumerate(cell_vertices):
                    outfile.write("%s %s" % cell_vertex)
                    if i < len(cell_vertices) - 1:
                        outfile.write(" ")
                outfile.write("\n")

            # Output the edge integers for each cell
            outfile.write('EDGES\n')
            outfile.write('%d\n' % (len(edge_boundary_list) * 4))
            for edge_boundary_indices in edge_boundary_list:
                for edge_index in edge_boundary_indices:
                    outfile.write("%d\n" % int(edge_index))

    def generate_mesh(self, output_file):
        """
        Generate and output the mesh file from the JSON
        block specification.

        Parameters
        ----------
        output_file : `string`
            The name of the output mesh file.
        """
        blocks = self._parse_block_json_file()
        cell_vertices_list, edge_boundary_list = self._discretise_domain_into_cells(blocks)
        self._output_mesh(output_file, cell_vertices_list, edge_boundary_list)


@click.command()
@click.option('--block-json-file', 'block_json_file', help='Full path to the block JSON file specifying the mesh blocks.')
@click.option('--output-file', 'output_file', help='Full path to the mesh output file.')
def cli(block_json_file, output_file):
    meshgen = BlockStructuredMeshGenerator(block_json_file)
    meshgen.generate_mesh(output_file) 	
 	

if __name__ == "__main__":
    cli() 

