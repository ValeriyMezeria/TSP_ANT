#include "tsp.h"

int main(int argc, char** argv)
{
    TSP a;

    if (argc != 2)
    {
       std::cout << "Wrong arguments number:" << argc << std::endl;
       return -1;
    }

    if (!a.readFromFile(argv[1]))
    {
        std::cout << "Task read failed!" << std::endl;
        return -1;
    }

    std::cout << "Name: " << a.getName() << std::endl;
    std::cout << "Description: " << a.getDescription() << std::endl;
    std::cout << "Size: " << a.getSize() << std::endl;

    a.solve();

    return 0;
}
