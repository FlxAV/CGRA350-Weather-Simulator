
#pragma once

#include <glm/glm.hpp>


using namespace std;
using namespace cgra;
using namespace glm;

class CloudModel {
public:
	// CloudModel(); // Constructor
	// ~CloudModel(); // Destructor

	 //// Initialize the cloud model with necessary parameters
	// void initialize();

	CloudModel(mesh_builder meshb, int resolution, glm::vec3 localUp);  //, GLuint cloudDensityTextureID
	void printNoiseMap();
	void createMesh_v1(int numPoints, float size, int seed, float verticalSpread);
	void createMesh_v2(float threshold, float gradualFactor);
	void createMesh_v4(float threshold, float gradualFactor, float amp, float freq, float nHeight);
	void perlinMap_v3(float amplitude, float frequency, float noiseHeight);
	void createMesh_v3(float threshold, float gradualFactor);
	float customCloudShapeFunction(const vec2& uv);
	//void createMesh(int numSpheres, float cloudRadius, float spacing, int seed);
	//void createMesh();
   // void generateVolumetricClouds(int resolution, float size, int seed);
	//void createMesh();
    void perlinMap();
	void perlinMap_v2();
	void NoiseMap(int resolution, float scale, int octaves, float persistance, float lacunarity, float seed);

	void draw(const glm::mat4& view, const glm::mat4& proj);


	// Function to calculate Worley noise value at a given 2D position
	float calculateWorleyNoise(const glm::vec2& position);

	//float calculateWorleyNoise3(const glm::vec3& position);



	// gl_mesh createSphere(float radius, int slices, int stacks);

	// gl_mesh createPlane(float width, float depth, int segmentsX, int segmentsZ);

	// void generateCloudGeometry(int resolution, float size, int seed);

	mesh_builder meshb;
	int resolution;
	glm::vec3 localUp;
	glm::vec3 axisA;
	glm::vec3 axisB;
	GLuint shader = 0;
	cgra::gl_mesh mesh;
	glm::vec3 color{ 230, 130, 230 };
	glm::mat4 modelTransform{ 1.0 };
	std::vector<float> noiseMap;

	GLuint normalMapTextureID = 0;







private:


	//glm::vec3 m_color{ 1.0f, 1.0f, 1.0f }; // Cloud color
	//float m_density = 0.5f; // Cloud density (adjust as needed)

	// You can add more parameters and data members as required
};
