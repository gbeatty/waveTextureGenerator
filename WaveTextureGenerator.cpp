
#include "tessendorf.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <strstream>
#include <iostream>
#include <iomanip>
#include <ppl.h> 
#include <atomic>

int main()
{
    const unsigned size = 1024;
    std::atomic<int> count = 0;

    Concurrency::parallel_for(0, 200,
        [&size, &count](int i)
    {
        double time = i;

        std::cout << "count = " << count++ << std::endl;
        
        tessendorf tess(0.000001, 20, Cartesian3(0.0, -1.0, 0.0), .02, time / 4.0, 50.0, size, size, 4000, 4000, .1, 1.0);

        Cartesian3* vertices = tess.simulate();

        std::vector<uint8_t> data;
        for (unsigned y = 0; y < size; ++y)
        {
            for (unsigned x = 0; x < size; ++x)
            {
                size_t index = (y * size) + x;
                double height = vertices[index].z / 2.0;
                uint8_t adjustedHeight = static_cast<uint8_t>(height + 128.0);
                data.push_back(adjustedHeight);
            }
        }

        std::stringstream fileName;
        fileName << "E:/waves/WaveTextureGenerator/x64/Release/HeightField" << std::setfill('0') << std::setw(3) << time << ".png";
        std::string name = fileName.str();

        stbi_write_png(name.c_str(), size, size, 1, &data[0], 0);
        //stbi_write_bmp(name.c_str(), size, size, 1, &data[0]);
    });

    /*for(double time = 0.0; time<200; ++time)
    {
        std::cout << "time = " << time << std::endl;
        tessendorf tess(0.000001, 20, Cartesian3(0.0, -1.0, 0.0), .02, time / 4.0, 50.0, size, size, 4000, 4000, .1, 1.0);

        Cartesian3* vertices = tess.simulate();

        std::vector<uint8_t> data;
        for (unsigned y = 0; y<size; ++y)
        {
            for (unsigned x = 0; x<size; ++x)
            {
                size_t index = (y * size) + x;
                double height = vertices[index].z / 2.0;
                uint8_t adjustedHeight = static_cast<uint8_t>(height + 128.0);
                data.push_back(adjustedHeight);
            }
        }

        std::stringstream fileName;
        fileName << "E:/waves/WaveTextureGenerator/x64/Release/HeightField" << std::setfill('0') << std::setw(3) << time << ".png";
        std::string name = fileName.str();

        stbi_write_png(name.c_str(), size, size, 1, &data[0], 0);
        //stbi_write_bmp(name.c_str(), size, size, 1, &data[0]);
    }*/
    
    return 0;
}

