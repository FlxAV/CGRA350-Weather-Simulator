#pragma once

using namespace std;
using namespace cgra;
using namespace glm;

class Terrain {
	
	
public:
	Terrain(mesh_builder meshb, int resolution, glm::vec3 localUp);
	void createMesh();
	void perlinMap();
	void draw();
	void NoiseMap(int resolution, float scale, int octaves, float persistance, float lacunarity, float seed);
	
	void draw(const glm::mat4& view, const glm::mat4 proj);


	GLuint normalMapTextureID = 0;
	mesh_builder meshb;
	int resolution;
	glm::vec3 localUp;
	glm::vec3 axisA;
	glm::vec3 axisB;
	GLuint m_texture;
	GLuint shader = 0;
	cgra::gl_mesh mesh;
	glm::vec3 color{ 230, 130, 230 };
	glm::mat4 modelTransform{ 1.0 };
	std::vector<float> noiseMap;
	std::vector<float> belowThreshold;


	std::vector<vec3> vertices;
};
