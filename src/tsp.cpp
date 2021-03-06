#include "tsp.h"

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <chrono>
#include <random>

namespace {

std::string ltrim(const std::string& str)
{
  const std::string pattern = " \f\n\r\t\v";
  return str.substr(str.find_first_not_of(pattern));
}

std::string rtrim(const std::string& str)
{
  const std::string pattern = " \f\n\r\t\v";
  return str.substr(0,str.find_last_not_of(pattern) + 1);
}

std::string trim(const std::string& str)
{
  return ltrim(rtrim(str));
}

int rand(int a, int b)
{
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<int> dist(a, b);

    return dist(mt);
}

TYPE str2type(const std::string& str)
{
    if (str == "TSP")
        return TYPE::TSP;
    else if (str == "ATSP")
        return TYPE::ATSP;
    else
        return TYPE::NONE;
}

EDGE_WEIGHT_TYPE str2edgeWeightType(const std::string& str)
{
    if (str == "EXPLICIT")
        return EDGE_WEIGHT_TYPE::EXPLICIT;
    else if (str == "EUC_2D")
        return EDGE_WEIGHT_TYPE::EUC_2D;
    else if (str == "ATT")
        return EDGE_WEIGHT_TYPE::ATT;
    else
        return EDGE_WEIGHT_TYPE::NONE;
}

EDGE_WEIGHT_FORMAT str2edgeWeightFormat(const std::string& str)
{
    if (str == "FULL_MATRIX")
        return EDGE_WEIGHT_FORMAT::FULL_MATRIX;
    else
        return EDGE_WEIGHT_FORMAT::NONE;
}

SECTION str2section(const std::string& str)
{
    if (str == "EDGE_WEIGHT_SECTION")
        return SECTION::EDGE_WEIGHT_SECTION;
    else if (str == "NODE_COORD_SECTION")
        return SECTION::NODE_COORD_SECTION;
    else
        return SECTION::NONE;
}

}

bool TSP::readFromFile(const std::string& filename)
{
    std::ifstream ifs(filename.c_str());

    if (!ifs.is_open())
        return false;

    std::string tmp;

    while (getline(ifs, tmp))
    {
        std::string val;

        if (parseParam(tmp, "NAME", val))
            m_name = val;
        else if (parseParam(tmp, "COMMENT", val))
            m_comment = val;
        else if (parseParam(tmp, "DIMENSION", val))
            m_size = std::stoi(val);
        else if (parseParam(tmp, "EDGE_WEIGHT_TYPE", val))
            m_edgeWeightType = str2edgeWeightType(val);
        else if (parseParam(tmp, "TYPE", val))
            m_type = str2type(val);
        else if (parseParam(tmp, "EDGE_WEIGHT_FORMAT", val))
            m_edgeWeightFormat = str2edgeWeightFormat(val);


        if (isEOF(tmp))
            break;

        if (isSection(tmp))
        {
            if (str2section(tmp) == SECTION::EDGE_WEIGHT_SECTION)
                readMatrix(ifs);
            else if (str2section(tmp) == SECTION::NODE_COORD_SECTION)
                readNodeCoord(ifs);
        }
    }

    return true;
}

void TSP::solve(bool verbose)
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    int time = 0;

    for (size_t i = 0; i < m_size; i++)
        m_pheromone.push_back(std::vector<double>(m_size, PHEROMONE_INITIAL));

    start = std::chrono::system_clock::now();

    ANT();

    end = std::chrono::system_clock::now();
    time = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();

    if (verbose)
    {
        std::cout << "Elapsed time: " << time << " ms" << std::endl;
        std::cout << "Record lenght: " << m_record << std::endl;
        std::cout << "Path: ";
        for (auto i : m_path) std::cout << i << " ";
        std::cout << std::endl << std::endl;
    }
}

bool TSP::isSection(const std::string& str)
{
    return str2section(str) != SECTION::NONE;
}

bool TSP::isEOF(const std::string& str)
{
    return str == "EOF";
}

bool TSP::parseSection(const std::string& str)
{
    return str2section(str) != SECTION::NONE;
}

bool TSP::parseParam(const std::string& str, const std::string& name, std::string& result)
{
    size_t pos;

    if (str.find(name) == std::string::npos)
        return false;

    if ((pos = str.find(':')) == std::string::npos)
        return false;

    std::string val(str, pos + 1, str.size());
    result = trim(val);

    return true;
}

bool TSP::readMatrix(std::ifstream& ifs)
{
    for (size_t i = 0; i < m_size; i++)
    {
        m_matrix.push_back(std::vector<double>());
        for (size_t j = 0; j < m_size; j++)
        {
            double tmp;
            ifs >> tmp;

            if (!ifs)
                return false;

            m_matrix.back().push_back(tmp);
        }
    }

    return true;
}

bool TSP::readNodeCoord(std::ifstream& ifs)
{
    std::vector<std::pair<double, double> > coord;
    int n;

    for (size_t i = 0; i < m_size; i++)
    {
        std::pair<double, double> tmp;
        ifs >> n;
        ifs >> tmp.first;
        ifs >> tmp.second;

        coord.push_back(tmp);
    }

    for (size_t i = 0; i < m_size; i++)
    {
        m_matrix.push_back(std::vector<double>());
        for (size_t j = 0; j < m_size; j++)
            if (i == j)
                m_matrix.back().push_back(INF);
            else
                m_matrix.back().push_back(dist(coord[i], coord[j]));
    }

    return true;
}

double TSP::dist(const std::pair<double, double>& a, const std::pair<double, double>& b)
{
    if (m_edgeWeightType == EDGE_WEIGHT_TYPE::ATT)
    {
        double r = sqrt(((b.first - a.first) * (b.first - a.first) + (b.second - a.second) * (b.second - a.second)) / 10.0);
        return round(r) < r ? round(r) + 1.0 : round(r);
    }

    return sqrt((b.first - a.first) * (b.first - a.first) + (b.second - a.second) * (b.second - a.second));
}

double TSP::getLenght(const std::vector<double>& path)
{
    if (path.empty())
        return INF;

    double result = 0;

    for (size_t i = 0; i < m_size - 1; i++)
        result += m_matrix[path[i]][path[i + 1]];

    result += m_matrix[path[m_size - 1]][path[0]];

    return result;
}

void TSP::updatePheromone(const std::vector<double>& path)
{
    double pheromone = PHEROMONE_AMOUNT / getLenght(path);

    for (size_t i = 0; i < m_size; i++)
        for (size_t j = 0; j < m_size; j++)
            m_pheromone[i][j] *= (1 - PHEROMONE_VAPORIZATION);

    for (size_t i = 0; i < path.size() - 1; i++)
        m_pheromone[path[i]][path[i+1]] += pheromone;

    m_pheromone[path[m_size - 1]][path[0]] += pheromone;
}

double TSP::transitionProbability(const std::vector<bool>& ant_memory,size_t from, size_t to)
{
    double sum = 0;

    for (size_t i = 0; i < m_size; i++)

        sum += (ant_memory[i] ? 0 : pow(1.0 / (m_matrix[from][i] == 0 ? 1.0 : m_matrix[from][i]), Q) * pow(m_pheromone[from][i], P));

    return pow(1.0 / (m_matrix[from][to] == 0 ? 1.0 : m_matrix[from][to]), Q) * pow(m_pheromone[from][to], P) / sum;
}

void TSP::ANT()
{
    for (size_t it = 0; it < 1 + ITERATIONS / m_size; it++)
    {
        for (size_t v = 0; v < m_size; v++)
        {
            size_t cur = v;
            std::vector<bool> antMemory(m_size, false);
            std::vector<double> path;

            antMemory[cur] = true;
            path.push_back(cur);

            for (size_t i = 0; i < m_size - 1; i++)
            {
                double die = rand(0, 99) / 100.0;
                double probability = 0;

                for (size_t j = 0; j < m_size; j++)
                {
                    if (!antMemory[j])
                    {
                        if (die >= probability && die <= probability + transitionProbability(antMemory, cur, j))
                        {
                            cur = j;
                            path.push_back(j);
                            antMemory[j] = true;
                            break;
                        }
                        probability += transitionProbability(antMemory, cur, j);
                    }
                }
            }

            updatePheromone(path);

            if (getLenght(path) < getLenght(m_path))
            {
                m_path = path;
                m_record = getLenght(path);
            }
        }
    }
}
