#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <mpi.h>
#include <numeric>

class DSU {
public:
    DSU(int num_elements) {
        parent.resize(num_elements);
        std::iota(parent.begin(), parent.end(), 0);
        rank.resize(num_elements, 0);
    }

    int find_set(int vertex) const {
        if(parent[vertex] == vertex)
            return vertex;
        
        // Path compression
        parent[vertex] = find_set(parent[vertex]);
        return parent[vertex];
    }

    void union_set(int node1, int node2) {
        int parent1 = find_set(node1);
        int parent2 = find_set(node2);

        if(parent1 == parent2)
            return;

        // Union by rank
        if(rank[parent1] > rank[parent2])
            std::swap(parent1, parent2);

        parent[parent1] = parent2;
        rank[parent2]++;
    }

private:
    mutable std::vector<int> parent;
    std::vector<int> rank;
};

class WeightedGraph {
public:
    struct Edge {
        int from;
        int to;
        int weight;
    };

    WeightedGraph(int num_vertices, int num_edges) 
        : vertices(num_vertices), edges(num_edges) {
        edgeList.resize(num_edges);
    }

    static WeightedGraph readFromFile(const std::string& filename) {
        std::ifstream file(filename);
        if(!file)
            throw std::runtime_error("Couldn't open input file");

        int num_vertices, num_edges;
        file >> num_vertices >> num_edges;
        WeightedGraph graph(num_vertices, num_edges);

        for(int i = 0; i < num_edges; i++) {
            file >> graph.edgeList[i].from >> graph.edgeList[i].to >> graph.edgeList[i].weight;
        }

        return graph;
    }

    void prettyPrint() const {
        const int width = 15; 

        std::cout << std::string(50, '-') << "\n";
        std::cout << std::left
                << std::setw(width) << "From"
                << std::setw(width) << "To"
                << std::setw(width) << "Weight" << "\n";
        std::cout << std::string(50, '-') << "\n";

        for(const auto& edge : edgeList) {
            std::cout << std::left
                    << std::setw(width) << edge.from
                    << std::setw(width) << edge.to
                    << std::setw(width) << edge.weight << "\n";
        }
        std::cout << std::string(50, '-') << "\n";
    }


    int getVertices() const { return vertices; }
    int getEdges() const { return edges; }
    const std::vector<Edge>& getEdgeList() const { return edgeList; }
    std::vector<Edge>& getEdgeList() { return edgeList; }

private:
    int vertices;
    int edges;
    std::vector<Edge> edgeList;
};

class BoruvkaMST {
public:
    static WeightedGraph findMST(const WeightedGraph& graph) {
        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        bool parallel = size != 1;

        // Broadcast graph information
        int vertices = graph.getVertices();
        int edges = graph.getEdges();
        MPI_Bcast(&edges, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&vertices, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // Scatter edges for parallel processing
        int edgesPart = (edges + size - 1) / size;
        std::vector<WeightedGraph::Edge> edgeListPart(edgesPart);
        if (parallel) {
            scatterEdgeList(graph.getEdgeList(), edgeListPart, edges, edgesPart, size, rank);
        } else {
            edgeListPart = graph.getEdgeList();
        }

        // Initialize MST
        WeightedGraph mst(vertices, vertices - 1);
        DSU fragment(vertices);
        int edgesMST = 0;

        // Data structures for finding closest edges
        std::vector<WeightedGraph::Edge> closestEdge(vertices);
        std::vector<WeightedGraph::Edge> closestEdgeReceived;
        if(parallel) {
            closestEdgeReceived.resize(vertices);
        }

        // Main Boruvka's algorithm loop
        for(int i = 1; i < vertices && edgesMST < vertices - 1; i *= 2) {
            // Reset closest edges
            for(auto& edge : closestEdge) {
                edge.weight = std::numeric_limits<int>::max();
            }

            // Find closest edges
            for(int j = 0; j < edgesPart; j++) {
                const auto& currentEdge = edgeListPart[j];
                int set1 = fragment.find_set(currentEdge.from);
                int set2 = fragment.find_set(currentEdge.to);

                if(set1 != set2) {
                    for (int k = 0; k < 2; k++) {
                        int canonicalElement = k == 0 ? set1 : set2;
                        if (closestEdge[canonicalElement].weight > currentEdge.weight) {
                            closestEdge[canonicalElement] = currentEdge;
                        }
                    }
                }
            }

            if(parallel) {
                gatherClosestEdges(closestEdge, closestEdgeReceived, vertices, size, rank);
            }

            // Add new edges to MST
            for(int j = 0; j < vertices; j++) {
                if(closestEdge[j].weight != std::numeric_limits<int>::max()) {
                    int from = closestEdge[j].from;
                    int to = closestEdge[j].to;

                    if(fragment.find_set(from) != fragment.find_set(to)) {
                        if(rank == 0) {
                            mst.getEdgeList()[edgesMST] = closestEdge[j];
                        }
                        edgesMST++;
                        fragment.union_set(from, to);
                    }
                }
            }
        }

        return mst;
    }

private:
    static void scatterEdgeList(const std::vector<WeightedGraph::Edge>& edgeList,
                               std::vector<WeightedGraph::Edge>& edgeListPart,
                               int elements, int& elementsPart,
                               int size, int rank) {
        MPI_Scatter(edgeList.data(), elementsPart * sizeof(WeightedGraph::Edge),
                   MPI_BYTE, edgeListPart.data(), elementsPart * sizeof(WeightedGraph::Edge),
                   MPI_BYTE, 0, MPI_COMM_WORLD);

        if(rank == size - 1 && elements % elementsPart != 0) {
            elementsPart = elements % elementsPart;
        }

        if(elements / 2 + 1 < size && elements != size) {
            if(rank == 0) {
                std::cerr << "Unsupported size/process combination, exiting!\n";
            }
            MPI_Finalize();
            std::exit(EXIT_FAILURE);
        }
    }

    static void gatherClosestEdges(std::vector<WeightedGraph::Edge>& closestEdge,
                                  std::vector<WeightedGraph::Edge>& closestEdgeReceived,
                                  int vertices, int size, int rank) {
        MPI_Status status;
        for(int step = 1; step < size; step *= 2) {
            if(rank % (2 * step) == 0) {
                int from = rank + step;
                if(from < size) {
                    MPI_Recv(closestEdgeReceived.data(), vertices * sizeof(WeightedGraph::Edge),
                            MPI_BYTE, from, 0, MPI_COMM_WORLD, &status);

                    for(int i = 0; i < vertices; i++) {
                        if(closestEdgeReceived[i].weight < closestEdge[i].weight) {
                            closestEdge[i] = closestEdgeReceived[i];
                        }
                    }
                }
            } else if(rank % step == 0) {
                int to = rank - step;
                MPI_Send(closestEdge.data(), vertices * sizeof(WeightedGraph::Edge),
                        MPI_BYTE, to, 0, MPI_COMM_WORLD);
            }
        }
        MPI_Bcast(closestEdge.data(), vertices * sizeof(WeightedGraph::Edge),
                  MPI_BYTE, 0, MPI_COMM_WORLD);
    }
};

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    try {
        WeightedGraph graph(0, 0);
        if(rank == 0) {
            if(argc != 2) {
                throw std::runtime_error("Usage: " + std::string(argv[0]) + " <input_file>");
            }
            graph = WeightedGraph::readFromFile(argv[1]);
            std::cout << "Original Graph:\n";
            graph.prettyPrint();
        }

        double startTime = MPI_Wtime();
        WeightedGraph mst = BoruvkaMST::findMST(graph);
        
        if(rank == 0) {
            std::cout << "Minimum Spanning Tree (Boruvka):\n";
            mst.prettyPrint();

            long weightMST = 0;
            for(const auto& edge : mst.getEdgeList()) {
                weightMST += edge.weight;
            }

            std::cout << "MST weight: " << weightMST << "\n";
            std::cout << "Time elapsed: " << MPI_Wtime() - startTime << " s\n";
        }
    }
    catch(const std::exception& e) {
        if(rank == 0) {
            std::cerr << "Error: " << e.what() << "\n";
        }
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}
