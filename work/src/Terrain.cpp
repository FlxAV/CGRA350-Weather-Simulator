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

#include "Terrain.hpp"
#include <glm/gtc/noise.hpp>

#include <chrono>

using namespace std;
using namespace cgra;
using namespace glm;

vec3 computeNormal(const std::vector<float>& noiseMap, int x, int y, int resolution) {
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

Terrain::Terrain(mesh_builder meshb, int resolution, glm::vec3 localUp) {
    this->meshb = meshb;
    this->resolution = resolution;
    this->localUp = localUp;
    this->axisA = vec3(localUp.y, localUp.z, localUp.x);
    this->axisB = cross(localUp, axisA);
    this->color = vec3(0.2, 0.2, 0.8);
}

void Terrain::createMesh() {
    int vertexCount = resolution * resolution;
    perlinMap();
    std::vector<vec3> vertices(vertexCount);
    std::vector<vec3> normals(vertexCount, vec3(0.0f)); // Store normals
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
            vertices[i].y += noiseMap[i];
            normals[i] = vec3(0, 1, 0.0f);// computeNormal(noiseMap, x, y, resolution);

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

    // Add the new vertices, normals, and UVs to the mesh_builder
    for (int i = 0; i < vertexCount; i++) {
        meshb.push_vertex(cgra::mesh_vertex{ vertices[i], normals[i], uvs[i], belowThreshold[i] });
    }

    for (int i = 0; i < triangles.size(); i += 3) {
        meshb.push_indices({ triangles[i], triangles[i + 1], triangles[i + 2] });
    }
    mesh = meshb.build();
}

void Terrain::perlinMap() {
    // Create a Perlin noise map
    noiseMap = std::vector<float>(resolution * resolution);
    belowThreshold = std::vector<float>(resolution * resolution);

    for (int y = 0; y < resolution; y++) {
        for (int x = 0; x < resolution; x++) {
            int i = x + y * resolution;
            float amplitude = 10;
            float frequency = 0.005;
            float noiseHeight = 0;

            for (int o = 0; o < 6; o++) {
                float sampleX = x * frequency;
                float sampleY = y * frequency;
                float perlinValue = glm::perlin(vec2(sampleX, sampleY)) * 2 - 1;
                noiseHeight += perlinValue * amplitude;
                amplitude *= 0.5f;
                frequency *= 2.0f;
            }

            if (noiseHeight < -25) {
                noiseHeight = -25;
                belowThreshold[i] = 1.0f;
            }

            noiseMap[i] = noiseHeight;
        }
    }
}

void Terrain::draw(const glm::mat4& view, const glm::mat4 proj) {
    mat4 modelview = view * modelTransform;
    glUseProgram(shader);
    glUniformMatrix4fv(glGetUniformLocation(shader, "uProjectionMatrix"), 1, false, value_ptr(proj));
    glUniformMatrix4fv(glGetUniformLocation(shader, "uModelViewMatrix"), 1, false, value_ptr(modelview));
    glUniform3fv(glGetUniformLocation(shader, "uColor"), 1, value_ptr(color));

    float currentTime = std::chrono::duration<float, std::chrono::seconds::period>(std::chrono::steady_clock::now().time_since_epoch()).count();
    glUniform1f(glGetUniformLocation(shader, "uTime"), currentTime);

    mesh.draw();
}