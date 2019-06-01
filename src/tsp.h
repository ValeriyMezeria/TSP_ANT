#include <iostream>
#include <vector>

const double INF = 999999;
const size_t ITERATIONS = 100;
const double Q = 2.5; // level of "greedy"
const double P = 0.9; // level of "herdness"
const double PHEROMONE_AMOUNT = 100.0; // pheromone amount coefficient secreted by each ant
const double PHEROMONE_VAPORIZATION = 0.15; // speed of pheromone evaporation
const double PHEROMONE_INITIAL = 50.0; // initial amount of pheromone on each transition

enum class TYPE { NONE = -1, TSP, ATSP };
enum class EDGE_WEIGHT_TYPE { NONE = -1, EXPLICIT, EUC_2D, ATT };
enum class EDGE_WEIGHT_FORMAT { NONE = -1, FULL_MATRIX };
enum class SECTION { NONE = -1, EDGE_WEIGHT_SECTION, NODE_COORD_SECTION };

class TSP
{
public:
    std::string getName() const { return m_name; }
    std::string getDescription() const { return m_comment; }
    int getSize() const { return m_size; }
    void showMatrix() { for (auto i : m_matrix) { std::cout << std::endl; for (auto j : i) std::cout << j << " "; } std::cout << std::endl; }
    bool readFromFile(const std::string& filename);
    void solve(bool verbose = true);

private:
    std::string m_name;
    std::string m_comment;

    TYPE m_type;
    EDGE_WEIGHT_TYPE m_edgeWeightType;
    EDGE_WEIGHT_FORMAT m_edgeWeightFormat;

    size_t m_size;

    double m_record;

    std::vector<std::vector<double>> m_matrix;
    std::vector<std::vector<double>> m_pheromone;
    std::vector<double> m_path;

    size_t m_iterations;

    //parsing
    bool isSection(const std::string& str);
    bool isEOF(const std::string& str);
    bool parseSection(const std::string& str);
    bool parseParam(const std::string& str, const std::string& name, std::string& result);
    bool readMatrix(std::ifstream& ifs);
    bool readNodeCoord(std::ifstream& ifs);
    double dist(const std::pair<double, double>& a, const std::pair<double, double>& b);

    //solving
    double getLenght(const std::vector<double>& path);
    void updatePheromone(const std::vector<double>& path);
    double transitionProbability(const std::vector<bool>& ant_memory,size_t from, size_t to);
    void ANT();
};
