#include "MeshObject/Application.hpp"
#include "MeshObject/InteractiveViewer.hpp"

#include "point_type/VilensPointT.hpp"
#include "point_type/BagPointT.hpp"
#include "MeshObject/Settings.hpp"

void headless_mode()
{
    // application
    Application<VilensPointT> app;

    // process
    app.load_point_cloud();
    app.process_the_rest();
    app.write_mesh();
}

void display_mode()
{
    // application
    Application<VilensPointT> app;

    // interactive viewer
    InteractiveViewer<VilensPointT> iviewer(app);
}


// using InputPointT = VilensPointT;
using InputPointT = BagPointT;
int main()
{
    // Fixed seed
    std::srand(30); 
    std::cout << std::rand() << std::endl;

    // mode
    Settings settings;
    if (settings.headless_mode)
    {
        headless_mode();
    }
    else
    {
        display_mode();
    }

   return 0;
}