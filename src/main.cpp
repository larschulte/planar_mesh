#include "MeshObject/Application.hpp"
#include "MeshObject/InteractiveViewer.hpp"

using InputPointT = VilensPointT;
int main()
{
    std::srand(30); // Fixed seed
    // test by print
    std::cout << std::rand() << std::endl;

    // application
    Application<InputPointT> app;

    // interactive viewer
    InteractiveViewer<InputPointT> iviewer(app);

   return 0;
}