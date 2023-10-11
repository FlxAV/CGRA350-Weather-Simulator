
#pragma once

// glm
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

// project
#include "opengl.hpp"
#include "cgra/cgra_mesh.hpp"
#include "skeleton_model.hpp"
#include "Terrain.hpp"
#include "cloudModel.hpp"

// Basic model that holds the shader, mesh and transform for drawing.
// Can be copied and modified for adding in extra information for drawing
// including textures for texture mapping etc.

//////////////////////// SCENE MODELS ////////////////////////
struct splineSeg {
	glm::vec3 a;
	glm::vec3 b;
	glm::vec3 c;
	glm::vec3 d;
};

struct basic_model {
	GLuint shader = 0;
	cgra::gl_mesh mesh;
	glm::vec3 color{ 0.7 };
	glm::mat4 modelTransform{ 1.0 };
	GLuint texture;
	std::vector<glm::vec3> interpolatedPts;

	void Teapot(const glm::mat4& view, const glm::mat4 proj, const std::vector<glm::vec3>& sliderPts, const float& start);
	void CoordPts(const glm::mat4& view, const glm::mat4 proj);
	// void draw2(const glm::mat4 &view, const glm::mat4 proj);
	void drawSplineLine(const glm::mat4& view, const glm::mat4 proj);

	void placePoints();
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
	Object(int ntype, const glm::vec3& nposition, const glm::vec3& nscale, Material& nmaterial);

};



//////////////////////// LIGHTING MODELS ////////////////////////

// Point light
struct PointLight {
	glm::vec3 position;
	glm::vec3 direction;
	float radius;
	glm::vec3 color;
	float power;
	float reach; // Only points within this distance of the light will be affected
	PointLight();
	PointLight(const glm::vec3& nposition, const glm::vec3& ndirection, float nradius, const glm::vec3& ncolor, float power, float reach);
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
	GLFWwindow* m_window;

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

	// Terrain
	glm::vec3 bottom = glm::vec3(0, 1, 0);
	mesh_builder floorMesh;
	int floorRes = 1000;
	Terrain plane = Terrain(floorMesh, floorRes, bottom);

	// Clouds
	float uLateralShift = 0.0;
	float uTransparency = 0.5;
	glm::vec3 top = glm::vec3(0, -1, 0);
	mesh_builder skymesh;
	int skyres = 1000;
	CloudModel clouds = CloudModel(skymesh, skyres, top);
	float m_threshold = 50.0;
	float m_gradualFactor = 40.0;

	float threshold_check = 0.0;
	float gradual_check = 0.0;

	float m_amp = 6.3;
	float m_freq = 0.01;
	float m_nHeight = -3.0;

	float amp_check = 0.0;
	float freq_check = 0.0;
	float nHeight_check = 0.0;

	bool redraw = false;

	// Raytracing fields
	glm::vec3 lightTranslate = glm::vec3(0, 100, 0);
	bool refreshRequired = false;
	int accumulatedPasses = 0;
	glm::vec3 camPos;
	int shadowResolution = 20;
	int lightBounces = 5;
	int framePasses = 10;
	float skyboxStrength = 1.0F;
	float skyboxGamma = 2.2F;
	float skyboxCeiling = 10.0F;
	

	float brightness = 1.0;

	GLuint skyboxTexture;
	GLuint m_skymap_shader;


	// Camera fields
	basic_model m_spline;       ///Teapot Spline
	basic_model m_cam_spline;   ///Camera Spline
	std::vector<glm::vec3> controlPts;
	std::vector<glm::vec3> interpolatedPoints;
	std::vector<glm::vec3> camControlPts;
	std::vector<glm::vec3> camInterpolatedPoints;
	std::vector<glm::vec3> sliderPts;
	std::vector<glm::vec3> camSliderPts;

public:
	// setup
	Application(GLFWwindow*);

	// ray - shader cache
	GLuint rayshader = 0;
	GLuint boundShader;
	Material planeMaterial;
	Material waterMaterial;

	// camera
	float cameraYaw = 0.0f, cameraPitch = 0.0f;

	int width, height;
	float vx = 0;
	float vy = 0;


	// disable copy constructors (for safety)
	Application(const Application&) = delete;
	Application& operator=(const Application&) = delete;

	//CORE
	cgra::mesh_builder drawUVSphere();
	cgra::mesh_builder drawPlane();

	// methods for the raytracing section
	void buildRayAdvancedShader();
	void buildRayBasicShader();
	void drawBasicScene(const glm::mat4& view, const glm::mat4 proj, double time);
	void drawSkybox(const glm::mat4& view, const glm::mat4 proj);


	// rendering callbacks (every frame)
	void render();
	void renderGUI();


	// input callbacks
	void cursorPosCallback(double xpos, double ypos);
	void mouseButtonCallback(int button, int action, int mods);
	void scrollCallback(double xoffset, double yoffset);
	void keyCallback(int key, int scancode, int action, int mods);
	void charCallback(unsigned int c);


	// Camera methods
	void calculateSpline(glm::vec3& p0, glm::vec3& p1, glm::vec3& p2, glm::vec3& p3, std::vector<glm::vec3>& intPoints, std::vector<glm::vec3>& sliderPoints);
	void animateTeapot();
	void plotPoints();
	int start = 0;
	int camStart = 0;
	int syncStart = 0;
	bool animate = false;
	bool camTrack = false;

};