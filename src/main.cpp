#include "MeshObject/Application.hpp"
#include "MeshObject/InteractiveViewer.hpp"

#include "point_type/VilensPointT.hpp"
#include "point_type/BagPointT.hpp"
#include "MeshObject/Settings.hpp"

// using InputPointT = VilensPointT;
using InputPointT = BagPointT;

void headless_mode()
{
    // application
    Application<InputPointT> app;
    Settings settings;

    // process
    if (settings.num_scans == -1)
    {
        app.process_the_rest();
    }
    else
    {
        for (int i = 0; i < settings.num_scans; i++)
        {
            app.load_pointcloud_from_dataloader();
            app.process_pointcloud();
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