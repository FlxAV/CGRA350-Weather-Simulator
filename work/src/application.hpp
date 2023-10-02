
#pragma once

// glm
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

// project
#include "opengl.hpp"
#include "cgra/cgra_mesh.hpp"
#include "skeleton_model.hpp"


// Basic model that holds the shader, mesh and transform for drawing.
// Can be copied and modified for adding in extra information for drawing
// including textures for texture mapping etc.

//////////////////////// SCENE MODELS ////////////////////////

struct basic_model {
	GLuint shader = 0;
	cgra::gl_mesh mesh;
	glm::vec3 color{0.7};
	glm::mat4 modelTransform{1.0};
	GLuint texture;

	void draw(const glm::mat4 &view, const glm::mat4 proj);
};

struct Material {
	glm::vec3 albedo;
	glm::vec3 specular;
	glm::vec3 emission;
	float emissionStrength;
	float roughness;
	float specularHighlight;
	float specularExponent;
	Material();
	Material(const glm::vec3& nalbedo, const glm::vec3& nspecular, const glm::vec3& nemission, float emissionStrength, float roughness, float specularHighlight, float specularExponent);
};


struct basic_rayplane {
	cgra::gl_mesh mesh;
	glm::vec3 color{ 0.7 };
	glm::mat4 modelTransform{ 1.0 };
	GLuint texture;
	Material material;

	void draw(const glm::mat4& view, const glm::mat4 proj);
};

struct basic_raymodel {
	cgra::gl_mesh mesh;
	glm::vec3 color{ 0.7 };
	glm::mat4 modelTransform{ 1.0 };
	GLuint texture;
	Material material;

	void draw(const glm::mat4& view, const glm::mat4 proj);
};

// PLACEHOLDER FOR MORE COMPLEX MODEL STRUCTS
struct Object {
	int type; // Type 0 = none (invisible), Type 1 = sphere, Type 2 = box
	glm::vec3 position;
	glm::vec3 scale; // For spheres, only the x value will be used as the radius
	Material material;
	Object();
	Object(int ntype, const glm::vec3& nposition, const glm::vec3& nscale, Material & nmaterial);

};



//////////////////////// LIGHTING MODELS ////////////////////////

// Point light
struct PointLight {
	glm::vec3 position;
	float radius;
	glm::vec3 color;
	float power;
	float reach; // Only points within this distance of the light will be affected
	PointLight();
	PointLight(const glm::vec3& nposition, float nradius, const glm::vec3& ncolor, float power, float reach);
};

// Sun light
struct SunLight {
	glm::vec3 position;
	glm::vec3 color;
	float power;
	float reach; // Only points within this distance of the light will be affected

	SunLight(const std::initializer_list<float>& position, float radius, const std::initializer_list<float>& color, float power, float reach);
	SunLight();
};

//////////////////////// Main application class ////////////////////////

class Application {
private:
	// window
	glm::vec2 m_windowsize;
	GLFWwindow *m_window;

	// oribital camera
	float m_pitch = .86;
	float m_yaw = -.86;
	float m_distance = 20;
	glm::mat4 proj, view;

	// last input
	bool m_leftMouseDown = false;
	glm::vec2 m_mousePosition;

	// drawing flags
	bool m_show_axis = false;
	bool m_show_grid = false;
	bool m_showWireframe = false;

	//Asssignment
	int numLatitudeLines = 30;
	int numLongitudeLines = 30;
	float radius = 1.0f;


	// geometry
	basic_model m_model;
	basic_model m_plane;;
	basic_raymodel m_rayplane;
	PointLight pointLight;


	bool refreshRequired = false;
	int accumulatedPasses = 0;
	glm::mat4 preView = glm::mat4(1);
	//glm::vec3 testpos = vec3(0);
	glm::vec3 camPos;

public:
	// setup
	Application(GLFWwindow *);

	// ray - shader cache
	GLuint directOutPassUniformLocation, accumulatedPassesUniformLocation, timeUniformLocation, camPosUniformLocation, rotationMatrixUniformLocation, aspectRatioUniformLocation, debugKeyUniformLocation;
	GLuint rayshader = 0;
	GLuint boundShader;
	std::vector<Object> objects;
	std::vector<PointLight> lights;
	Material planeMaterial;

	glm::vec3 cameraPosition = glm::vec3(0, 1, 2);
	float cameraYaw = 0.0f, cameraPitch = 0.0f;
	int shadowResolution = 20;
	int lightBounces = 5;
	int framePasses = 10;

	float blur = 0.002f; // Slight blur (les than a pixel) = anti-aliasing
	float bloomRadius = 0.02f;
	float bloomIntensity = 0.5f;
	bool planeVisible = true;

	int selectedObjectIndex = -1;
	int width, height;


	// disable copy constructors (for safety)
	Application(const Application&) = delete;
	Application& operator=(const Application&) = delete;

	//CORE
	cgra::mesh_builder drawUVSphere();
	cgra::mesh_builder drawPlane();

	// draw the rays in the shaders
	void buildRayShader();
	void buildRayBasicShader();
	void drawRayShader();
	void recompileShader();

	//object data
	void sendObjectData(int objectIndex);
	void bind(GLuint shaderProgram);
	void unbind();
	void drawScene(const glm::mat4& view, const glm::mat4 proj);
	bool sphereIntersection(glm::vec3 position, float radius, glm::vec3 rayOrigin, glm::vec3 rayDirection, float* hitDistance);
	bool boxIntersection(glm::vec3 position, glm::vec3 size, glm::vec3 rayOrigin, glm::vec3 rayDirection, float* hitDistance);
	bool planeIntersection(glm::vec3 planeNormal, glm::vec3 planePoint, glm::vec3 rayOrigin, glm::vec3 rayDirection, float* hitDistance);
	void placeBasicScene();
	void drawBasicScene(const glm::mat4& view, const glm::mat4 proj, double time);
	// rendering callbacks (every frame)
	void render();
	void renderGUI();

	// input callbacks
	void cursorPosCallback(double xpos, double ypos);
	void mouseButtonCallback(int button, int action, int mods);
	void scrollCallback(double xoffset, double yoffset);
	void keyCallback(int key, int scancode, int action, int mods);
	void charCallback(unsigned int c);
};