#include "MeshObject/Application.hpp"
#include "MeshObject/InteractiveViewer.hpp"

#include "point_type/VilensPointT.hpp"
#include "point_type/BagPointT.hpp"
#include "MeshObject/Settings.hpp"

// using InputPointT = VilensPointT;
using InputPointT = pcl::PointXYZ;

void headless_mode()
{
    // application
    Application<InputPointT> app;
    Settings settings;

    // process
    app.load_point_cloud();
    if (settings.num_scans == -1)
    {
        app.process_the_rest();
    }
    else
    {
        for (int i = 0; i < settings.num_scans - 1; i++)
        {
            app.loop();
        }
    }
    app.write_mesh();
}

void display_mode()
{
    // application
    Application<InputPointT> app;

    // interactive viewer
    InteractiveViewer<InputPointT> iviewer(app);
}

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
