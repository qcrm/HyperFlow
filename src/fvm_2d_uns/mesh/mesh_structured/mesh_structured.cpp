#ifndef __MESH_STRUCTURED_CPP
#define __MESH_STRUCTURED_CPP

#include "mesh_structured.h"

namespace HyperFlow {

/* Constructor */
LoadStructuredMesh::LoadStructuredMesh()
{} 

/* Constructor with filename */
LoadStructuredMesh::LoadStructuredMesh(const std::string& _filename)
: 
    filename(_filename)
{} 

/* Destructor */
LoadStructuredMesh::~LoadStructuredMesh()
{}

/* Load the HyperFlor native format structured mesh from disk */
std::shared_ptr<Mesh> LoadStructuredMesh::load_structured_mesh()
{
    std::ifstream hmsh;
    
    hmsh.open(filename);

    if (!hmsh) {
        std::cout << "Unable to open HyperFlow mesh file. Exiting." << std::endl;
        exit(1); // terminate with error
    }

    /* Skip header info */
    std::string str;
    std::getline(hmsh, str);
    std::getline(hmsh, str);
    unsigned int cell_count = std::stoul(str, nullptr, 0);
    std::cout << "CELL COUNT: " << cell_count << std::endl;

    /* Add the cell information */
    CellVec1D cells;

    while (std::getline(hmsh, str)) {
        if (str.compare("EDGES") == 0) {
            break;
        }

        std::istringstream iss(str);

        std::vector<std::string> string_nums;
        std::copy(std::istream_iterator<std::string>(iss),
                  std::istream_iterator<std::string>(),
                  std::back_inserter(string_nums));

        double v1x = std::stod(string_nums[0]);
        double v1y = std::stod(string_nums[1]);
        Vec1D vertex1 = {v1x, v1y};

        double v2x = std::stod(string_nums[2]);
        double v2y = std::stod(string_nums[3]);
        Vec1D vertex2 = {v2x, v2y};

        double v3x = std::stod(string_nums[4]);
        double v3y = std::stod(string_nums[5]);
        Vec1D vertex3 = {v3x, v3y};

        double v4x = std::stod(string_nums[6]);
        double v4y = std::stod(string_nums[7]);
        Vec1D vertex4 = {v4x, v4y};

        Vec2D cell_vertices {vertex1, vertex2, vertex3, vertex4};
        Cell cell(cell_vertices);
        cells.push_back(cell);
    }

    /* Load in edge boundary types */
    std::getline(hmsh, str);
    unsigned int edge_count = std::stoul(str, nullptr, 0);

    std::cout << "EDGE COUNT: " << edge_count << std::endl;
    std::vector<std::vector<int> > cell_edge_indices;
    EdgeVec1D edges;
    std::vector<int> edge_types;

    while (std::getline(hmsh, str)) {
        int edge_type = std::stoi(str, nullptr, 0);
        edge_types.push_back(edge_type);
    }

    /* Generate edges */
    int outer_edge_idx = 0;
    for (int cell_idx=0; cell_idx<cells.size(); cell_idx++) {
        std::cout << "CELL INDEX " << cell_idx + 1 << " OF " << cells.size() << std::endl;
        unsigned int edge_size = cells[cell_idx].get_edges().size();

        for (unsigned int edge_idx=0; edge_idx<edge_size; edge_idx++) {
            Edge edge = cells[cell_idx].get_edges()[edge_idx];
            std::vector<int> cell_edge_idx = {cell_idx, edge_types[outer_edge_idx], outer_edge_idx};

            for (unsigned int nbrc_idx=0; nbrc_idx<cells.size(); nbrc_idx++) {
                unsigned int nbrc_edge_size = cells[nbrc_idx].get_edges().size();    

                for (unsigned int nbrc_edge_idx=0; nbrc_edge_idx<nbrc_edge_size; nbrc_edge_idx++) {
                    if (edge.opposite_edges_equal(cells[nbrc_idx].get_edges()[nbrc_edge_idx])) {
                        cell_edge_idx[1] = nbrc_idx;
                    }
                }
            }

            cell_edge_indices.push_back(cell_edge_idx);
            edges.push_back(edge);
            outer_edge_idx++;
        }
    } 

    hmsh.close();

    std::shared_ptr<Mesh> mesh = std::make_shared<Mesh>(edges, cells, cell_edge_indices);
    std::cout << "COMPLETED MESH" << std::endl;
    return mesh;
}

}

#endif
