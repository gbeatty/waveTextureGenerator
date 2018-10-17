#pragma once

#include "Cartesian3.h"

class VertexData
{
public:

    Cartesian3 Vertex;
    Cartesian3 Normal;

    VertexData() = default;

    VertexData(Cartesian3 vertex, Cartesian3 normal) :
        Vertex(vertex),
        Normal(normal)
    {
        
    }
};