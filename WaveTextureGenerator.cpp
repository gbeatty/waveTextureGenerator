
#include "tessendorf.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <strstream>
#include <iostream>
#include <iomanip>
#include <ppl.h> 
#include <atomic>
#include <filesystem>

using namespace std::experimental;

int main()
{
    const auto outputDir = filesystem::current_path().append("WaveTextures");
    filesystem::remove_all(outputDir);
    filesystem::create_directory(outputDir);

    std::atomic<int> count = 0;
	std::mutex mutex;

	const unsigned size = 1024;
	const int period = 20;
	const int samplesPerSecond = 10;
	const int numberOfSamples = period * samplesPerSecond;

#ifdef _DEBUG
    for (int i = 0; i < numberOfSamples; ++i)
#else
    Concurrency::parallel_for(0, numberOfSamples, [&](int i)
#endif
    {
		const double timePerSample = 1.0 / static_cast<double>(samplesPerSecond);
        const double time = i * timePerSample;

		std::unique_lock<std::mutex> lock(mutex);
        std::cout << "count = " << count++ << std::endl;
		lock.unlock();
        
        tessendorf tess(0.000001, 20, Cartesian3(0.0, -1.0, 0.0), .02, time, period, size, size, 4000, 4000, .1, 1.0);

        std::vector<VertexData> vertices = tess.simulate();

        std::vector<uint8_t> heightFieldPixels;
        std::vector<uint8_t> normalVectorPixels;
        for (unsigned y = 0; y < size; ++y)
        {
            for (unsigned x = 0; x < size; ++x)
            {
                size_t index = (y * size) + x;
                double height = vertices[index].Vertex.z;
                uint8_t adjustedHeight = static_cast<uint8_t>(height + 128.0);
                heightFieldPixels.push_back(adjustedHeight);

                Cartesian3 scaledNormal = vertices[index].Normal * 128.0;
                Cartesian3 shiftedScaledNormal = scaledNormal + 128.0;
                normalVectorPixels.push_back(shiftedScaledNormal.x);
            }
        }

        std::stringstream heightfieldFileName;
        heightfieldFileName << "HeightField_" << std::setfill('0') << std::setw(3) << i << ".png";
        auto heightfieldPath = outputDir / heightfieldFileName.str();

        std::stringstream normalFileName;
        normalFileName << "NormalVector_" << std::setfill('0') << std::setw(3) << i << ".png";
        auto normalPath = outputDir / normalFileName.str();

        stbi_write_png(heightfieldPath.string().c_str(), size, size, 1, &heightFieldPixels[0], 0);
        stbi_write_png(normalPath.string().c_str(), size, size, 1, &normalVectorPixels[0], 0);
    }
#ifndef _DEBUG
    );
#endif
	
    auto ffmpegPath = filesystem::current_path() / "ffmpeg.exe";

    std::stringstream heightfieldCmd;
    heightfieldCmd << ffmpegPath.string() << " -y -start_number 0 -framerate " << samplesPerSecond << " -i ./WaveTextures/HeightField_%03d.png -c:v libx264 -vf \"fps = " << samplesPerSecond << ", format = yuv420p\" HeightField.mp4";
    std::cout << heightfieldCmd.str() << std::endl;
    std::system(heightfieldCmd.str().c_str());

    std::stringstream normalVectorCmd;
    normalVectorCmd << ffmpegPath.string() << " -y -start_number 0 -framerate " << samplesPerSecond << " -i ./WaveTextures/NormalVector_%03d.png -c:v libx264 -vf \"fps = " << samplesPerSecond << ", format = yuv420p\" NormalVector.mp4";
    std::system(normalVectorCmd.str().c_str());

    return 0;
}

