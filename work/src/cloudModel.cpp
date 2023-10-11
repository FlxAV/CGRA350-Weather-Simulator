
// std
#include <iostream>
#include <string>
#include <chrono>

// glm
#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>

// project
#include "application.hpp"
#include "cgra/cgra_geometry.hpp"
#include "cgra/cgra_gui.hpp"
#include "cgra/cgra_image.hpp"
#include "cgra/cgra_shader.hpp"
#include "cgra/cgra_wavefront.hpp"

#include "cloudModel.hpp"
#include <glm/gtc/noise.hpp>
#include <random>
#include <glm/gtc/random.hpp>

#include <vector>
#include <glm/glm.hpp>


using namespace std;
using namespace cgra;
using namespace glm;


CloudModel::CloudModel(mesh_builder meshb, int resolution, glm::vec3 localUp) {  //GLuint cloudDensityTextureID
    this->meshb = meshb;
    this->resolution = resolution;
    this->localUp = localUp;
    this->axisA = vec3(localUp.y, localUp.z, localUp.x);
    this->axisB = cross(localUp, axisA);



}

vec3 computeNormalCloud(const std::vector<float>& noiseMap, int x, int y, int resolution) {
    // Edge cases
    float hL = x > 0 ? noiseMap[(x - 1) + y * resolution] : noiseMap[x + y * resolution];
    float hR = x < resolution - 1 ? noiseMap[(x + 1) + y * resolution] : noiseMap[x + y * resolution];
    float hU = y > 0 ? noiseMap[x + (y - 1) * resolution] : noiseMap[x + y * resolution];
    float hD = y < resolution - 1 ? noiseMap[x + (y + 1) * resolution] : noiseMap[x + y * resolution];

    // Tangents
    vec3 tangentRight = vec3(1, hR - noiseMap[x + y * resolution], 0);
    vec3 tangentUp = vec3(0, hU - noiseMap[x + y * resolution], 1);

    // Cross product
    vec3 normal = cross(tangentUp, tangentRight);
    return normalize(normal);
}


void CloudModel::createMesh_v4(float threshold, float gradualFactor, float amp, float freq, float nHeight) {
    int vertexCount = resolution * resolution;
    std::vector<vec3> vertices(vertexCount);
    std::vector<vec3> normals(vertexCount); // Store normals
    std::vector<vec2> uvs(vertexCount); // Store UVs
    std::vector<unsigned int> triangles(6 * (resolution - 1) * (resolution - 1));

    int triIndex = 0;

    perlinMap_v3(amp, freq, nHeight);

    for (int y = 0; y < resolution; y++) {
        for (int x = 0; x < resolution; x++) {
            int i = x + y * resolution;
            vec2 percent = vec2(static_cast<float>(x), static_cast<float>(y));
            percent.x /= static_cast<float>(resolution - 1);
            percent.y /= static_cast<float>(resolution - 1);

            vec3 pointOnUnitCube = localUp + (percent.x - 0.5f) * 2 * axisA + (percent.y - 0.5f) * 2 * axisB;

            vertices[i] = pointOnUnitCube;
          
            // Compute the normal using the updated computeNormal function
            normals[i] = computeNormalCloud(noiseMap, x, y, resolution);

            // Compute UVs based on percent
            uvs[i] = percent;

            if (x != resolution - 1 && y != resolution - 1) {
                triangles[triIndex] = i;
                triangles[triIndex + 1] = i + resolution + 1;
                triangles[triIndex + 2] = i + resolution;

                triangles[triIndex + 3] = i;
                triangles[triIndex + 4] = i + 1;
                triangles[triIndex + 5] = i + resolution + 1;
                triIndex += 6;
            }
        }
    }
    meshb.clear();

    float maxHeight = -20;

    // Adjust the Y value based on noiseMap and threshold
    for (int i = 0; i < vertices.size(); i++) {
        float noiseValue = noiseMap[i];
        if (noiseValue >= threshold) {
            vertices[i].y = maxHeight;
        }
        else if (noiseValue >= threshold - gradualFactor) {
            // Gradually increase Y value for vertices close to the threshold
            float blendFactor = (threshold - noiseValue) / gradualFactor;
            vertices[i].y = maxHeight * (1.0 - blendFactor);
        }
    }

    // Add the new vertices, normals, and UVs to the mesh_builder
    for (int i = 0; i < vertexCount; i++) {
        meshb.push_vertex(cgra::mesh_vertex{ vertices[i], normals[i], uvs[i] });
    }

    for (int i = 0; i < triangles.size(); i += 3) {
        meshb.push_indices({ triangles[i], triangles[i + 1], triangles[i + 2] });
    }
    mesh = meshb.build();
}

void CloudModel::perlinMap_v3(float amp,float freq, float nHeight) {
    // Create a Perlin noise map
    noiseMap = std::vector<float>(resolution * resolution);

    // Seed the random number generator based on the current time
    unsigned int seed = static_cast<unsigned int>(std::time(0));
    std::srand(seed);

    for (int y = 0; y < resolution; y++) {
        for (int x = 0; x < resolution; x++) {
            int i = x + y * resolution;
            float amplitude = amp;  // Increase the initial amplitude
            float frequency = freq;    // Decrease the initial frequency
            float noiseHeight = nHeight;

            for (int o = 0; o < 6; o++) {  // Increase the number of octaves
                float sampleX = x * frequency;
                float sampleY = y * frequency;

                float perlinValue = glm::perlin(vec2(sampleX, sampleY)) * 2 - 1;
                noiseHeight -= perlinValue * amplitude;

                amplitude *= 0.5f;  // You can use a factor smaller than 1 to decrease amplitude
                frequency *= 2.0f;  // You can use a factor larger than 1 to increase frequency
            }

            noiseMap[i] = noiseHeight;
        }
    }
}




void CloudModel::draw(const glm::mat4& view, const glm::mat4& proj) {

    mat4 modelview = view * modelTransform;

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


    glUseProgram(shader); // load shader and variables
    glUniformMatrix4fv(glGetUniformLocation(shader, "uProjectionMatrix"), 1, false, value_ptr(proj));
    glUniformMatrix4fv(glGetUniformLocation(shader, "uModelViewMatrix"), 1, false, value_ptr(modelview));
    glUniform3fv(glGetUniformLocation(shader, "uColor"), 1, value_ptr(color));

    float currentTime = std::chrono::duration<float, std::chrono::seconds::period>(std::chrono::steady_clock::now().time_since_epoch()).count();
    glUniform1f(glGetUniformLocation(shader, "uTime"), currentTime);

    mesh.draw(); // draw

    glDisable(GL_BLEND);
	
}


void CloudModel::createMesh_v2(float threshold, float gradualFactor) {
    int vertexCount = resolution * resolution;
    std::vector<vec3> vertices(vertexCount);
    std::vector<vec3> normals(vertexCount); // Store normals
    std::vector<vec2> uvs(vertexCount); // Store UVs
    std::vector<unsigned int> triangles(6 * (resolution - 1) * (resolution - 1));

    int triIndex = 0;

    for (int y = 0; y < resolution; y++) {
        for (int x = 0; x < resolution; x++) {
            int i = x + y * resolution;
            vec2 percent = vec2(static_cast<float>(x), static_cast<float>(y));
            percent.x /= static_cast<float>(resolution - 1);
            percent.y /= static_cast<float>(resolution - 1);

            vec3 pointOnUnitCube = localUp + (percent.x - 0.5f) * 2 * axisA + (percent.y - 0.5f) * 2 * axisB;

            vertices[i] = pointOnUnitCube;

            // Compute the normal by taking the cross product of the axisA and axisB vectors
            // and use it as the normal for each vertex
            normals[i] = normalize(cross(axisA, axisB));

            // Compute UVs based on percent
            uvs[i] = percent;

            if (x != resolution - 1 && y != resolution - 1) {
                triangles[triIndex] = i;
                triangles[triIndex + 1] = i + resolution + 1;
                triangles[triIndex + 2] = i + resolution;

                triangles[triIndex + 3] = i;
                triangles[triIndex + 4] = i + 1;
                triangles[triIndex + 5] = i + resolution + 1;
                triIndex += 6;
            }
        }
    }
    meshb.clear();

    perlinMap();

    //float threshold = -10.0;

    float maxHeight = 20;
    //float gradualFactor = 10;

    // Adjust the Y value based on noiseMap and threshold
    for (int i = 0; i < vertices.size(); i++) {
        float noiseValue = noiseMap[i];
        if (noiseValue >= threshold) {
            vertices[i].y = maxHeight;
        }
        else if (noiseValue >= threshold - gradualFactor) {
            // Gradually increase Y value for vertices close to the threshold
            float blendFactor = (threshold - noiseValue) / gradualFactor;
            vertices[i].y = maxHeight * (1.0 - blendFactor);
        }
    }

    // Add the new vertices, normals, and UVs to the mesh_builder
    for (int i = 0; i < vertexCount; i++) {
        meshb.push_vertex(cgra::mesh_vertex{ vertices[i], normals[i], uvs[i] });
    }

    for (int i = 0; i < triangles.size(); i += 3) {
        meshb.push_indices({ triangles[i], triangles[i + 1], triangles[i + 2] });
    }
    mesh = meshb.build();
}

void CloudModel::perlinMap() {
    // Create a Perlin noise map
    noiseMap = std::vector<float>(resolution * resolution);

    // Seed the random number generator based on the current time
    unsigned int seed = static_cast<unsigned int>(std::time(0));
    std::srand(seed);

    for (int y = 0; y < resolution; y++) {
        for (int x = 0; x < resolution; x++) {
            int i = x + y * resolution;
            float amplitude = 10;  // Increase the initial amplitude
            float frequency = 0.01f;    // Decrease the initial frequency
            float noiseHeight = 0;

            for (int o = 0; o < 6; o++) {  // Increase the number of octaves
                float sampleX = x * frequency;
                float sampleY = y * frequency;

                float perlinValue = glm::perlin(vec2(sampleX, sampleY)) * 2 - 1;
                noiseHeight += perlinValue * amplitude;

                amplitude *= 0.5f;  // You can use a factor smaller than 1 to decrease amplitude
                frequency *= 2.0f;  // You can use a factor larger than 1 to increase frequency
            }

            noiseMap[i] = noiseHeight;
        }
    }
}

void CloudModel::printNoiseMap() {
    perlinMap();

    for (int i = 0; i < resolution; i++) {
        cout << noiseMap[i] << "  |  ";
    }

}

void CloudModel::createMesh_v1(int numPoints, float size, int seed, float verticalSpread) {
    // Create a vector to hold the vertices
    std::vector<glm::vec3> vertices;
    std::vector<unsigned int> indices;

    // Use a random number generator with the given seed
    std::mt19937 mt(seed);
    std::uniform_real_distribution<float> randCoord(-size * 0.5f, size * 0.5f);
    std::uniform_real_distribution<float> randVertical(-verticalSpread, verticalSpread);

    // Generate random 3D points within the bounding volume
    for (int i = 0; i < numPoints; ++i) {
        glm::vec3 point(randCoord(mt), randVertical(mt), randCoord(mt));
        vertices.push_back(point);
    }

    // Create connections between points to form triangles in 3D space
    for (int i = 0; i < numPoints - 2; ++i) {
        for (int j = i + 1; j < numPoints - 1; ++j) {
            for (int k = j + 1; k < numPoints; ++k) {
                // Add triangles using the three selected points
                indices.push_back(i);
                indices.push_back(j);
                indices.push_back(k);
            }
        }
    }

    // Create the mesh
    mesh_builder meshb;

    // Add the vertices to the mesh_builder
    for (const auto& vertex : vertices) {
        meshb.push_vertex(cgra::mesh_vertex{ vertex, glm::vec3(0.0f), glm::vec2(0.0f) });
    }

    // Add the indices to the mesh_builder
    for (const auto& index : indices) {
        meshb.push_index(index);
    }

    // Build the mesh
    mesh = meshb.build();
}

void CloudModel::createMesh_v3(float threshold, float gradualFactor) {
    int vertexCount = resolution * resolution;
    std::vector<vec3> vertices(vertexCount);
    std::vector<vec3> normals(vertexCount); // Store normals
    std::vector<vec2> uvs(vertexCount); // Store UVs
    std::vector<unsigned int> triangles(6 * (resolution - 1) * (resolution - 1));

    int triIndex = 0;

    for (int y = 0; y < resolution; y++) {
        for (int x = 0; x < resolution; x++) {
            int i = x + y * resolution;
            vec2 percent = vec2(static_cast<float>(x), static_cast<float>(y));
            percent.x /= static_cast<float>(resolution - 1);
            percent.y /= static_cast<float>(resolution - 1);

            // Use Perlin noise to displace the vertices and create a random shape
            float noiseValue = glm::perlin(glm::vec2(percent.x * 5.0f, percent.y * 5.0f)) * 0.5f;
            float cloudShape = 1.0f + noiseValue; // Randomize the cloud shape

            vec3 pointOnUnitCube = localUp + (percent.x - 0.5f) * 2 * axisA + (percent.y - 0.5f) * 2 * axisB;

            // Scale the point based on the cloud shape
            pointOnUnitCube *= cloudShape;

            vertices[i] = pointOnUnitCube;

            // Compute the normal by taking the cross product of the axisA and axisB vectors
            // and use it as the normal for each vertex
            normals[i] = normalize(cross(axisA, axisB));

            // Compute UVs based on percent
            uvs[i] = percent;

            if (x != resolution - 1 && y != resolution - 1) {
                triangles[triIndex] = i;
                triangles[triIndex + 1] = i + resolution + 1;
                triangles[triIndex + 2] = i + resolution;

                triangles[triIndex + 3] = i;
                triangles[triIndex + 4] = i + 1;
                triangles[triIndex + 5] = i + resolution + 1;
                triIndex += 6;
            }
        }
    }
    meshb.clear();

    perlinMap_v2();

    float maxHeight = 50;

    // Adjust the Y value based on noiseMap and threshold
    for (int i = 0; i < vertices.size(); i++) {
        float noiseValue = noiseMap[i];
        if (noiseValue >= threshold) {
            vertices[i].y = maxHeight;
        }
        else if (noiseValue >= threshold - gradualFactor) {
            // Gradually increase Y value for vertices close to the threshold
            float blendFactor = (threshold - noiseValue) / gradualFactor;
            vertices[i].y = maxHeight * (1.0 - blendFactor);
        }
    }

    // Add the vertices, normals, and UVs to the mesh_builder
    for (int i = 0; i < vertexCount; i++) {
        meshb.push_vertex(cgra::mesh_vertex{ vertices[i], normals[i], uvs[i] });
    }

    for (int i = 0; i < triangles.size(); i += 3) {
        meshb.push_indices({ triangles[i], triangles[i + 1], triangles[i + 2] });
    }
    mesh = meshb.build();
}

void CloudModel::perlinMap_v2() {
    // Create a Perlin noise map
    noiseMap = std::vector<float>(resolution * resolution);

    // Calculate the center of the noise map
    float centerX = static_cast<float>(resolution - 1) / 2.0f;
    float centerY = static_cast<float>(resolution - 1) / 2.0f;

    // Seed the random number generator based on the current time
    unsigned int seed = static_cast<unsigned int>(std::time(0));
    std::srand(seed);

    for (int y = 0; y < resolution; y++) {
        for (int x = 0; x < resolution; x++) {
            int i = x + y * resolution;

            // Calculate the distance from the center
            float distance = glm::distance(vec2(x, y), vec2(centerX, centerY));

            // Define amplitude and frequency based on distance
            float amplitude = 10.0f * (1.0f - distance / (centerX > centerY ? centerX : centerY)); // Adjust amplitude based on distance
            float frequency = 0.01f + 0.1f * (distance / (centerX > centerY ? centerX : centerY)); // Adjust frequency based on distance

            float noiseHeight = 0;

            for (int o = 0; o < 6; o++) {  // Increase the number of octaves
                float sampleX = x * frequency;
                float sampleY = y * frequency;

                float perlinValue = glm::perlin(vec2(sampleX, sampleY)) * 2 - 1;

                // Invert the Perlin value and apply a custom function to add randomness
                float invertedPerlinValue = 1.0f - (perlinValue + 1.0f) / 2.0f;
                float normalizedDistance = distance / (centerX > centerY ? centerX : centerY);
                float smoothFactor = 1.0f - glm::smoothstep(0.0f, 1.0f, normalizedDistance);

                // Introduce randomness to the noiseHeight
                float randomness = (static_cast<float>(std::rand()) / RAND_MAX) * 0.7f; // Adjust the value to control randomness
                noiseHeight += (invertedPerlinValue + randomness) * amplitude * smoothFactor;

                amplitude *= 0.5f;  // You can use a factor smaller than 1 to decrease amplitude
                frequency *= 2.0f;  // You can use a factor larger than 1 to increase frequency
            }

            noiseMap[i] = noiseHeight;
        }
    }
}

// Function to calculate Worley noise value at a given 2D position
float CloudModel::calculateWorleyNoise(const glm::vec2& position) {
    // Define the number of feature points (cell centers)
    int numPoints = 10; // You can adjust this value

    // Initialize variables to keep track of the closest and second closest distances
    float closestDist = 1000.0f; // Set to a large initial value
    float secondClosestDist = 1000.0f; // Set to a large initial value

    // Loop through each feature point
    for (int i = 0; i < numPoints; ++i) {
        // Generate a random feature point (cell center)
        glm::vec2 featurePoint = glm::vec2(glm::linearRand(-1.0f, 1.0f), glm::linearRand(-1.0f, 1.0f));

        // Calculate the distance between the current position and the feature point
        float distance = glm::length(position - featurePoint);

        // Update the closest and second closest distances
        if (distance < closestDist) {
            secondClosestDist = closestDist;
            closestDist = distance;
        }
        else if (distance < secondClosestDist) {
            secondClosestDist = distance;
        }
    }

    // Calculate the Worley noise value based on the closest and second closest distances
    return closestDist - secondClosestDist;
}



//vec3 computeNormal_v2(const std::vector<float>& noiseMap, const std::vector<int>& isAboveThreshold, int x, int y, int resolution) {
//    // Ensure x and y are within valid bounds
//    x = std::max(0, std::min(x, resolution - 1));
//    y = std::max(0, std::min(y, resolution - 1));
//
//    // Check if all three vertices of the triangle are above the threshold
//    if (isAboveThreshold[x + y * resolution] == 1 &&
//        isAboveThreshold[(x + 1) + y * resolution] == 1 &&
//        isAboveThreshold[x + (y + 1) * resolution] == 1) {
//        // Calculate the normal as you did previously
//        float hL = noiseMap[x + y * resolution];
//        float hR = noiseMap[(x + 1) + y * resolution];
//        float hU = noiseMap[x + (y + 1) * resolution];
//        vec3 tangentRight = vec3(1, hR - hL, 0);
//        vec3 tangentUp = vec3(0, hU - hL, 1);
//        vec3 normal = normalize(cross(tangentUp, tangentRight));
//        return normal;
//    }
//    else {
//        // Return a default normal (e.g., pointing up) for triangles below the threshold
//        return vec3(0.0, 1.0, 0.0);
//    }
//}

//int computeIsAboveThreshold(const vec3& position, float threshold) {
//    // Compare the y-coordinate of the vertex position with the threshold
//    return (position.y >= threshold) ? 1 : 0;
//}

//float CloudModel::calculateWorleyNoise3(const glm::vec3& position) {
//    // Define the number of feature points (cell centers)
//    int numPoints = 10; // You can adjust this value
//
//    // Initialize variables to keep track of the closest and second closest distances
//    float closestDist = 1000.0f; // Set to a large initial value
//    float secondClosestDist = 1000.0f; // Set to a large initial value
//
//    // Loop through each feature point
//    for (int i = 0; i < numPoints; ++i) {
//        // Generate a random feature point (cell center)
//        glm::vec3 featurePoint = glm::linearRand(vec3(-1.0f), vec3(1.0f));
//
//        // Calculate the distance between the current position and the feature point
//        float distance = glm::length(position - featurePoint);
//
//        // Update the closest and second closest distances
//        if (distance < closestDist) {
//            secondClosestDist = closestDist;
//            closestDist = distance;
//        }
//        else if (distance < secondClosestDist) {
//            secondClosestDist = distance;
//        }
//    }
//
//    // Calculate the Worley noise value based on the closest and second closest distances
//    return closestDist - secondClosestDist;
//}

//void CloudModel::generateCloudGeometry(int resolution, float size, int seed) {
//    // Create a mesh_builder to build the mesh
//    cgra::mesh_builder builder;
//
//    // Set random seed
//    std::default_random_engine rng(seed);
//
//    // Loop to create vertices with positions on a sphere
//    for (int i = 0; i < resolution; ++i) {
//        for (int j = 0; j < resolution; ++j) {
//            // Calculate the spherical coordinates
//            float theta = glm::mix(0.0f, glm::two_pi<float>(), static_cast<float>(i) / (resolution - 1));
//            float phi = glm::mix(0.0f, glm::pi<float>(), static_cast<float>(j) / (resolution - 1));
//
//            // Calculate the position on the sphere
//            float x = size * sin(phi) * cos(theta);
//            float y = size * cos(phi);
//            float z = size * sin(phi) * sin(theta);
//
//            // Randomly displace the position within a bounding box
//            glm::vec3 displacement = glm::linearRand(glm::vec3(-size * 0.1f), glm::vec3(size * 0.1f));
//            glm::vec3 vertexPos = glm::vec3(x, y, z) + displacement;
//
//            // Compute a random normal vector (simulating the cloud's surface)
//            glm::vec3 normal = glm::normalize(glm::sphericalRand(1.0f));
//
//            // Generate random UV coordinates (optional)
//            glm::vec2 uv = glm::linearRand(glm::vec2(0.0f), glm::vec2(1.0f));
//
//            // Add the vertex to the builder
//            builder.push_vertex(cgra::mesh_vertex{ vertexPos, normal, uv });
//        }
//    }
//
//    // Generate triangles based on vertices (e.g., create a grid)
//    for (int i = 0; i < resolution - 1; ++i) {
//        for (int j = 0; j < resolution - 1; ++j) {
//            // Calculate indices for two triangles forming a quad
//            unsigned int idx0 = i * resolution + j;
//            unsigned int idx1 = idx0 + 1;
//            unsigned int idx2 = (i + 1) * resolution + j;
//            unsigned int idx3 = idx2 + 1;
//
//            initializer_list<unsigned int> inds = { idx0, idx1, idx2 };
//
//            // Define triangles (counter-clockwise winding order)
//            builder.push_indices(inds);
//            inds = { idx2, idx1, idx3 };
//            builder.push_indices(inds);
//        }
//    }
//
//    // Build the mesh
//    mesh = builder.build();
//}
//
//
//
//cgra::gl_mesh CloudModel::createSphere(float radius, int slices, int stacks) {
//    cgra::mesh_builder builder;
//
//    for (int i = 0; i <= stacks; ++i) {
//        float V = i / (float)stacks;
//        float phi = V * glm::pi<float>();
//
//        for (int j = 0; j <= slices; ++j) {
//            float U = j / (float)slices;
//            float theta = U * (glm::pi<float>() * 2);
//
//            float X = cos(theta) * sin(phi);
//            float Y = cos(phi);
//            float Z = sin(theta) * sin(phi);
//
//            // Create a mesh_vertex and set its position
//            cgra::mesh_vertex vertex;
//            vertex.pos = glm::vec3(X, Y, Z) * radius;
//
//            // Add the vertex to the builder
//            builder.push_vertex(vertex);
//        }
//    }
//
//    for (int i = 0; i < stacks; ++i) {
//        for (int j = 0; j < slices; ++j) {
//            int first = (i * (slices + 1)) + j;
//            int second = first + slices + 1;
//
//            builder.push_index(first);
//            builder.push_index(second);
//            builder.push_index(first + 1);
//
//            builder.push_index(second);
//            builder.push_index(second + 1);
//            builder.push_index(first + 1);
//        }
//    }
//
//    return builder.build();
//}
//
//cgra::gl_mesh CloudModel::createPlane(float width, float depth, int segmentsX, int segmentsZ) {
//    cgra::mesh_builder builder;
//
//    float segmentWidth = width / segmentsX;
//    float segmentDepth = depth / segmentsZ;
//
//    for (int i = 0; i <= segmentsZ; ++i) {
//        for (int j = 0; j <= segmentsX; ++j) {
//            float x = j * segmentWidth - width / 2.0f;
//            float z = i * segmentDepth - depth / 2.0f;
//
//            // Create a mesh_vertex and set its position
//            cgra::mesh_vertex vertex;
//            vertex.pos = glm::vec3(x, 0.0f, z);
//
//            // Add the vertex to the builder
//            builder.push_vertex(vertex);
//        }
//    }
//
//    for (int i = 0; i < segmentsZ; ++i) {
//        for (int j = 0; j < segmentsX; ++j) {
//            int first = (i * (segmentsX + 1)) + j;
//            int second = first + segmentsX + 1;
//
//            builder.push_index(first);
//            builder.push_index(second);
//            builder.push_index(first + 1);
//
//            builder.push_index(second);
//            builder.push_index(second + 1);
//            builder.push_index(first + 1);
//        }
//    }
//
//    return builder.build();
//}
//
//
//void CloudModel::generateVolumetricClouds(int resolution, float size, int seed) {
//    // Create a mesh_builder to build the volumetric cloud mesh
//    cgra::mesh_builder builder;
//
//    // Set random seed
//    std::default_random_engine rng(seed);
//
//    // Generate volumetric cloud geometry
//    // Modify this part to create a cloud volume
//    for (int i = 0; i < resolution; ++i) {
//        for (int j = 0; j < resolution; ++j) {
//            for (int k = 0; k < resolution; ++k) {
//                // Generate a random position within the cloud volume
//                glm::vec3 vertexPos = glm::linearRand(glm::vec3(-size), glm::vec3(size));
//
//                // Compute a random normal vector (if needed)
//                glm::vec3 normal = glm::normalize(glm::sphericalRand(1.0f));
//
//                // Add the vertex to the builder
//                builder.push_vertex(cgra::mesh_vertex{ vertexPos, normal, /* add UV here */ });
//            }
//        }
//    }
//
//    // Generate triangles based on vertices
//    // Modify this part to create triangles that form the cloud volume
//    for (int i = 0; i < resolution - 1; ++i) {
//        for (int j = 0; j < resolution - 1; ++j) {
//            for (int k = 0; k < resolution - 1; ++k) {
//                // Calculate indices for a voxel or other cloud representation
//                unsigned int idx0 = i * resolution * resolution + j * resolution + k;
//                unsigned int idx1 = idx0 + 1;
//                unsigned int idx2 = (i + 1) * resolution * resolution + j * resolution + k;
//                unsigned int idx3 = idx2 + 1;
//                unsigned int idx4 = i * resolution * resolution + (j + 1) * resolution + k;
//                unsigned int idx5 = idx4 + 1;
//                unsigned int idx6 = (i + 1) * resolution * resolution + (j + 1) * resolution + k;
//                unsigned int idx7 = idx6 + 1;
//
//                // Define triangles (counter-clockwise winding order)
//                builder.push_indices({ idx0, idx1, idx2 });
//                builder.push_indices({ idx2, idx1, idx3 });
//                builder.push_indices({ idx2, idx3, idx6 });
//                builder.push_indices({ idx6, idx3, idx7 });
//                // Add more triangles as needed for your volumetric representation
//            }
//        }
//    }
//
//    // Build the mesh
//    mesh = builder.build();
//}
