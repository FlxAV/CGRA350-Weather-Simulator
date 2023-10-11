
// std
#include <iostream>
#include <string>
#include <chrono>

#define GLM_ENABLE_EXPERIMENTAL

// glm
#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/spline.hpp>

#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>


// project
#include "application.hpp"
#include "cgra/cgra_geometry.hpp"
#include "cgra/cgra_gui.hpp"
#include "cgra/cgra_image.hpp"
#include "cgra/cgra_shader.hpp"
#include "cgra/cgra_wavefront.hpp"


using namespace std;
using namespace cgra;
using namespace glm;


#define M_PI 3.1415926535897932384626433832795


void basic_model::CoordPts(const glm::mat4& view, const glm::mat4 proj) {
	mat4 modelview = view * modelTransform;
	modelview = glm::scale(modelview, vec3(0.1));

	glUseProgram(shader); // load shader and variables
	glUniformMatrix4fv(glGetUniformLocation(shader, "uProjectionMatrix"), 1, false, value_ptr(proj));
	glUniformMatrix4fv(glGetUniformLocation(shader, "uModelViewMatrix"), 1, false, value_ptr(modelview));
	glUniform3fv(glGetUniformLocation(shader, "uColor"), 1, value_ptr(color));

	///Draw the 5 Sphere interpolation points
	//drawSphere();
}


void basic_model::drawSplineLine(const glm::mat4& view, const glm::mat4 proj) {
	mat4 modelview = view * modelTransform;

	glUseProgram(shader); // load shader and variables
	glUniformMatrix4fv(glGetUniformLocation(shader, "uProjectionMatrix"), 1, false, value_ptr(proj));
	glUniformMatrix4fv(glGetUniformLocation(shader, "uModelViewMatrix"), 1, false, value_ptr(modelview));
	glUniform3fv(glGetUniformLocation(shader, "uColor"), 1, value_ptr(color));

	mesh.mode = GL_LINES;
	mesh.draw(); // draw
}

/// <----------------------------- >
void basic_model::placePoints() {
	mesh_builder spline;


	for (int k = 0; k < interpolatedPts.size(); k++) {

		spline.push_vertex(mesh_vertex{ interpolatedPts[k] });
		//spline.push_vertex(mesh_vertex{catmullRom(splinePointsRaw[0+k],splinePointsRaw[1+k],splinePointsRaw[2+k],splinePointsRaw[3+k],i)});
	}
	for (int k = 0; k < interpolatedPts.size() - 1; k++) {

		spline.push_index(k);
		spline.push_index(k + 1);
		//spline.push_vertex(mesh_vertex{catmullRom(splinePointsRaw[0+k],splinePointsRaw[1+k],splinePointsRaw[2+k],splinePointsRaw[3+k],i)});
	}

	mesh = spline.build();

}


mesh_builder Application::drawUVSphere() {
	// One vertex at every latitude-longitude intersection,
	mesh_builder mb;
	mesh_vertex mv;


	float latStep = M_PI / numLatitudeLines;
	float longStep = 2 * M_PI / numLongitudeLines;

	float x, y, z;
	float nx, ny, nz;
	float u, v;

	float uStep = 1.0f / numLongitudeLines;
	float vStep = 1.0f / numLatitudeLines;
	float uCoord = 0;
	float vCoord = 0;

	for (int i = 0; i < numLatitudeLines + 1; i++) {
		for (int j = 0; j < numLongitudeLines + 1; j++) {

			// Positions
			x = radius * sin(latStep * i) * cos(longStep * j);
			y = radius * sin(latStep * i) * sin(longStep * j);
			z = radius * cos(latStep * i);
			vec3 positions = vec3(x, y, z);

			//	UVs
			u = uCoord;
			v = vCoord;
			vec2 uv = vec2(u, v);


			// Add the vertex to the sphere mesh
			mb.push_vertex({ positions, positions, uv });

			uCoord += uStep;
		}

		uCoord = 0.0f;
		vCoord += vStep;
	}

	// Add indices for the sphere faces
	for (int i = 0; i < numLatitudeLines + 1; i++) {
		for (int j = 0; j < numLongitudeLines + 1; j++) {

			int p1 = i * numLongitudeLines + j;
			int p2 = i * numLongitudeLines + j + 1;

			int p3 = (i + 1) * numLongitudeLines + j + 1;
			int p4 = (i + 1) * numLongitudeLines + j;

			// Triangle 1
			mb.push_index(p1);
			mb.push_index(p2);
			mb.push_index(p3);

			// Triangle 2
			mb.push_index(p1);
			mb.push_index(p3);
			mb.push_index(p4);
		}
	}
	return mb;
}
mesh_builder Application::drawPlane() {
	mesh_builder mb;
	mesh_vertex mv;
	std::vector<vec3> positions;
	std::vector<vec3> normals;
	std::vector<vec2> uvs;

	std::vector<unsigned int> indices;
	int subdivisions = 5; // You can adjust this value to control the plane's subdivisions

	int size = 100;
	int pointsPerRow = (int)pow(2, subdivisions) + 1;

	float u, v;

	// Generate vertices and indices for the plane
	for (int i = 0; i <= subdivisions; ++i) {
		float y = static_cast<float>(i) / static_cast<float>(subdivisions) - 0.5f;
		//v = (float)i / (pointsPerRow - 1);
		for (int j = 0; j <= subdivisions; ++j) {
			float x = static_cast<float>(j) / static_cast<float>(subdivisions) - 0.5f;

			// Calculate vertex position
			positions.push_back(vec3(x, 0.0f, y) * vec3(size) + vec3(0, -2, 0));

			// Calculate vertex normal (for a simple plane, this will be the same for all vertices)
			normals.push_back(vec3(0, 1, 0.0f));

			// Calculate vertex texture coordinates
			float u = (float)i / subdivisions;
			float v = (float)j / subdivisions;
			//u = (float)j / (pointsPerRow - 1);
			//std::cout << u << "\n";
			//std::cout << v << "\n";
			uvs.push_back(vec2(u, v));
		}
	}

	// Generate indices for the plane
	for (int i = 0; i < subdivisions; ++i) {

		for (int j = 0; j < subdivisions; ++j) {
			int v0 = i * (subdivisions + 1) + j;
			int v1 = v0 + 1;
			int v2 = (i + 1) * (subdivisions + 1) + j;
			int v3 = v2 + 1;

			// First triangle
			indices.push_back(v0);
			indices.push_back(v2);
			indices.push_back(v1);

			// Second triangle
			indices.push_back(v1);
			indices.push_back(v2);
			indices.push_back(v3);
		}
	}

	// Create mesh from indices
	for (unsigned int i = 0; i < indices.size(); ++i) {
		mb.push_index(i);
		mb.push_vertex(mesh_vertex{
			positions[indices[i]],
			normals[indices[i]],
			uvs[indices[i]]
			});
	}

	//lighting 

	return mb;
}

//////////////////////// SCENE MODELS ////////////////////////
Object::Object() = default;
Object::Object(int ntype, const glm::vec3 & nposition, const glm::vec3 & nscale, Material & nmaterial) {
	type = ntype;
	position = nposition;
	scale = nscale;
	material = nmaterial;

}




// Material
Material::Material() = default;
Material::Material(const glm::vec3 & nalbedo, const glm::vec3 & nspecular, const glm::vec3 & nemission, float emissionStrength, float roughness, float specularHighlight, float specularExponent) {
	albedo = nalbedo;
	specular = nspecular;
	emission = nemission;
	this->emissionStrength = emissionStrength;
	this->roughness = roughness;
	this->specularHighlight = specularHighlight;
	this->specularExponent = specularExponent;
}
//////////////////////// LIGHTING MODELS ////////////////////////
PointLight::PointLight() = default;
PointLight::PointLight(const glm::vec3 & nposition, const glm::vec3 & ndirection, float nradius, const glm::vec3 & ncolor, float power, float reach) {
	position = nposition;
	direction = ndirection;
	this->radius = nradius;
	color = ncolor;

	this->power = power;
	this->reach = reach;
}


//////////////////////// Main application class ////////////////////////

Application::Application(GLFWwindow * window) : m_window(window) {

	//shader_builder sb;
	//sb.set_shader(GL_VERTEX_SHADER, CGRA_SRCDIR + std::string("//res//shaders//terrain_vert.glsl"));
	//sb.set_shader(GL_FRAGMENT_SHADER, CGRA_SRCDIR + std::string("//res//shaders//terrain_frag.glsl"));
	//GLuint shader = sb.build();

	////m_model.shader = shader;
	////m_model.mesh = load_wavefront_data(CGRA_SRCDIR + std::string("/res//assets//teapot.obj")).build();
	////m_model.color = vec3(1, 0, 0);
	shader_builder cb;
	cb.set_shader(GL_VERTEX_SHADER, CGRA_SRCDIR + std::string("//res//shaders//color_vert_cloud.glsl"));
	cb.set_shader(GL_FRAGMENT_SHADER, CGRA_SRCDIR + std::string("//res//shaders//color_frag_cloud.glsl"));
	GLuint cloudshader = cb.build();

	shader_builder sb;
	sb.set_shader(GL_VERTEX_SHADER, CGRA_SRCDIR + std::string("//res//shaders//color_vert.glsl"));
	sb.set_shader(GL_FRAGMENT_SHADER, CGRA_SRCDIR + std::string("//res//shaders//color_frag.glsl"));
	GLuint shader = sb.build();


	//m_model.shader = shader;
	//m_model.mesh = drawUVSphere().build();

	//m_rayplane.mesh = drawPlane().build();
	planeMaterial = Material(vec3(0.6F, 0.2F, 0.2F), vec3(0.0F, 0.0F, 0.0F), vec3(0.0F, 0.0F, 0.0F), 0.0F, 0.0F, 0.0F, 0.0F);
	waterMaterial = Material(vec3(0.1F, 0.1F, 0.4F), vec3(0.25F, 0.25F, 0.75F), vec3(0.0F, 0.0F, 0.5F), 0.0F, 0.0F, 0.0F, 0.0F);
	//m_rayplane.material = planeMaterial;
	//placeBasicScene();
	//recompileShader();

	// Position , radius , color , power , reach
	pointLight = PointLight(vec3(1.0, 5.0, 1.0), vec3(-1, -1, 0), 6, vec3(0.7), 500, 500);
	buildRayBasicShader();

	//glUseProgram(rayshader);
	//glUniform1i(glGetUniformLocation(rayshader, "u_screenTexture"), 0);

	//plane.shader = shader;
	plane.resolution = 850;
	plane.modelTransform = glm::scale(glm::mat4(1), glm::vec3(150, 1, 150));
	plane.createMesh();
	//plane.color = vec3(1);


	clouds.shader = cloudshader;
	clouds.resolution = 1500;
	clouds.modelTransform = glm::translate(glm::mat4(1), glm::vec3(0, 55, 0)); // Move up by 2 units in the Y-direction
	clouds.modelTransform *= glm::scale(glm::mat4(1), glm::vec3(1000, 1, 1000));
	clouds.createMesh_v4(m_threshold, m_gradualFactor, m_amp, m_freq, m_nHeight);
	clouds.color = vec3(1);



	

	// camera stuff
	//glBegin(GL_LINES);
	m_spline.mesh.mode = GL_LINES;
	m_cam_spline.shader = shader;
	m_cam_spline.color = vec3(1, 1, 1);
	m_cam_spline.mesh.mode = GL_LINES;

	// points
	plotPoints();


	for (int i = 0; i < controlPts.size() - 3; i++) {
		calculateSpline(controlPts[i], controlPts[i + 1], controlPts[i + 2], controlPts[i + 3],
			interpolatedPoints, sliderPts);
	}

	for (int i = 0; i < camControlPts.size() - 3; i++) {
		calculateSpline(camControlPts[i], camControlPts[i + 1], camControlPts[i + 2],
			camControlPts[i + 3], camInterpolatedPoints, camSliderPts);
	}

	// calculate spline
	m_spline.interpolatedPts = interpolatedPoints;
	m_spline.shader = shader;
	m_spline.placePoints();

}


void Application::render() {

	//TIME
	double preTime = glfwGetTime();

	// retrieve the window hieght
	glfwGetFramebufferSize(m_window, &width, &height);

	m_windowsize = vec2(width, height); // update window size
	glViewport(0, 0, width, height); // set the viewport to draw to the entire window

	// clear the back-buffer
	glClearColor(0.3f, 0.3f, 0.4f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// enable flags for normal/forward rendering
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);

	// projection matrix - CHANGED Z FOR MORE DISTANCE 
	proj = perspective(1.f, float(width) / height, 0.1f, 100000000.f);

	// view matrix
	view = translate(mat4(1), vec3(vx, vy, -m_distance))
		* rotate(mat4(1), m_pitch, vec3(1, 0, 0))
		* rotate(mat4(1), m_yaw, vec3(0, 1, 0));


	camPos = vec3(0, 0, 0);

	if (camTrack) {
		m_pitch = -0.03;
		m_yaw = -0.08;
		//view = translate(mat4(1), -camSliderPts[start]);
		view = glm::lookAt(camSliderPts[camStart], sliderPts[start], vec3(0, 1, 0));
	}


	// helpful draw options
	if (m_show_grid) drawGrid(view, proj);
	if (m_show_axis) drawAxis(view, proj);
	glPolygonMode(GL_FRONT_AND_BACK, (m_showWireframe) ? GL_LINE : GL_FILL);



	
	drawBasicScene(view, proj, preTime);

	if (redraw) {
		clouds.createMesh_v4(m_threshold, m_gradualFactor, m_amp, m_freq, m_nHeight);
		redraw = false;  // Reset the flag
	}


	clouds.draw(view, proj);

}


void Application::renderGUI() {

	// setup window
	ImGui::SetNextWindowPos(ImVec2(5, 5), ImGuiSetCond_Once);
	ImGui::SetNextWindowSize(ImVec2(300, 200), ImGuiSetCond_Once);
	ImGui::Begin("Options", 0);

	// display current camera parameters
	ImGui::Text("Application %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
	ImGui::SliderFloat("Pitch", &m_pitch, -pi<float>() / 2, pi<float>() / 2, "%.2f");
	ImGui::SliderFloat("Yaw", &m_yaw, -pi<float>(), pi<float>(), "%.2f");
	ImGui::SliderFloat("Distance", &m_distance, 0, 100, "%.2f", 2.0f);
	ImGui::Separator();

	ImGui::SliderInt("framePasses", &framePasses, 0, 100);
	//ImGui::SliderInt("shadowResolution", &shadowResolution, 0, 100);
	ImGui::SliderInt("lightBounces", &lightBounces, 0, 15);
	ImGui::SliderFloat("light power", &pointLight.power, 0, 5000);
	ImGui::SliderFloat("light reach", &pointLight.reach, 0, 5000);
	ImGui::SliderFloat3("light position", value_ptr(lightTranslate), -1, 1);
	ImGui::SliderFloat("light Height", &lightTranslate.y, 0, 250);

	ImGui::SliderFloat3("light color", value_ptr(pointLight.color), 0, 1);
	ImGui::Separator();
	ImGui::SliderFloat3("Object albedo", value_ptr(planeMaterial.albedo), 0, 1);
	ImGui::SliderFloat3("Object specular", value_ptr(planeMaterial.specular), 0, 1);
	ImGui::SliderFloat3("Object emission", value_ptr(planeMaterial.emission), 0, 1);
	ImGui::SliderFloat("Object emissionStrength", &planeMaterial.emissionStrength, 0, 1);
	ImGui::SliderFloat("Object roughness", &planeMaterial.roughness, 0, 1);
	//ImGui::SliderFloat("Object spec Highlight", &planeMaterial.specularHighlight, 0, 1);
	//ImGui::SliderFloat("Object spec EXPO", &planeMaterial.specularExponent, 0, 11);
	// helpful drawing options
	ImGui::Checkbox("Show axis", &m_show_axis);
	ImGui::SameLine();
	ImGui::Checkbox("Show grid", &m_show_grid);
	ImGui::Checkbox("Wireframe", &m_showWireframe);
	ImGui::SameLine();
	if (ImGui::Button("Screenshot")) rgba_image::screenshot(true);


	ImGui::Separator();

	ImGui::SliderFloat("Threshold", &m_threshold, 0, 100, "%.2f", 2.0f);
	ImGui::SliderFloat("GradualFactor", &m_gradualFactor, 20, 120, "%.2f", 2.0f);
	ImGui::SliderFloat("Amp", &m_amp, -10, 15, "%.2f", 2.0f);
	ImGui::SliderFloat("Freq", &m_freq, -1, 1, "%.2f", 2.0f);
	ImGui::SliderFloat("nHeight", &m_nHeight, -20, 20, "%.2f", 2.0f);

	if (ImGui::Button("Redraw")) {
		redraw = true;  // Set the flag to trigger a redraw
	}

	ImGui::Separator();

	float animationSpeed = 1;

	if (ImGui::SliderInt("Camera", &syncStart, 0, camSliderPts.size() - 101)) {
		camStart = syncStart;
		start = syncStart;
	}
	if (ImGui::Button("Camera Track")) {
		int count = 0;
		if (camTrack && count == 0) {
			camTrack = false;
			count++;
		}
		else if (camTrack == false && count == 0) {
			camTrack = true;
			count++;
		}
	}

	// Update syncStart to simulate animation when in camTrack mode

	if (camTrack == true) {
		// Increment the slider value based on the animation speed
		syncStart += animationSpeed;
		start = syncStart;
		camStart = syncStart;

		// Ensure the slider value stays within the allowed range
		if (syncStart > camSliderPts.size() - 101) {
			//ImGui::SliderInt("Animate", &syncStart, 0, camSliderPts.size() - 101);
			syncStart = 0;  // Reset to the beginning if it exceeds the maximum value
		}
	}

	ImGui::Separator();

	// example of how to use input boxes
	static float exampleInput;
	if (ImGui::InputFloat("example input", &exampleInput)) {
		cout << "example input changed to " << exampleInput << endl;
	}

	// finish creating window
	ImGui::End();
}

// Builds ray shader program -> RAYTRACING SECTION
void Application::buildRayBasicShader() {

	shader_builder rayshaderbuilder;
	rayshaderbuilder.set_shader(GL_VERTEX_SHADER, CGRA_SRCDIR + std::string("//res//shaders//basic_raytracing_vert.glsl"));
	rayshaderbuilder.set_shader(GL_FRAGMENT_SHADER, CGRA_SRCDIR + std::string("//res//shaders//basic_raytracing_frag.glsl"));
	rayshader = rayshaderbuilder.build();

}

void Application::buildRayAdvancedShader() {

	shader_builder rayshaderbuilder;
	rayshaderbuilder.set_shader(GL_VERTEX_SHADER, CGRA_SRCDIR + std::string("//res//shaders//advanced_raytracing_vert.glsl"));
	rayshaderbuilder.set_shader(GL_FRAGMENT_SHADER, CGRA_SRCDIR + std::string("//res//shaders//advanced_raytracing_frag.glsl"));
	rayshader = rayshaderbuilder.build();

}

void Application::drawBasicScene(const glm::mat4 & view, const glm::mat4 proj, double time) {

	// DRAW PLANE
	mat4 modelview = view * plane.modelTransform;

	glUseProgram(rayshader); // load shader and variables
	accumulatedPasses += 1;

	glUniformMatrix4fv(glGetUniformLocation(rayshader, "uProjectionMatrix"), 1, false, value_ptr(proj));
	glUniformMatrix4fv(glGetUniformLocation(rayshader, "uModelViewMatrix"), 1, false, value_ptr(modelview));
	glUniform3fv(glGetUniformLocation(rayshader, "uColor"), 1, value_ptr(m_rayplane.color));
	glUniform1f(glGetUniformLocation(rayshader, "u_aspectRatio"), (float)width / height);
	glUniform1i(glGetUniformLocation(rayshader, "u_shadowResolution"), shadowResolution);
	glUniform1i(glGetUniformLocation(rayshader, "u_lightBounces"), lightBounces);
	glUniform1i(glGetUniformLocation(rayshader, "u_framePasses"), framePasses);
	glUniform1f(glGetUniformLocation(rayshader, "u_time"), (float)time);

	float currentTime = std::chrono::duration<float, std::chrono::seconds::period>(std::chrono::steady_clock::now().time_since_epoch()).count();
	glUniform1f(glGetUniformLocation(rayshader, "u_waveTime"), currentTime);

	glUniform3fv(glGetUniformLocation(rayshader, "u_cameraPosition"), 1, value_ptr(camPos));
	glUniform3fv(glGetUniformLocation(rayshader, "u_lightTranslation"), 1, value_ptr(lightTranslate));

	// light stuff
	glUniform3fv(glGetUniformLocation(rayshader, "u_light.position"), 1, value_ptr(pointLight.position));
	glUniform1f(glGetUniformLocation(rayshader, "u_light.radius"), pointLight.radius);
	glUniform3f(glGetUniformLocation(rayshader, "u_light.color"), pointLight.color.x, pointLight.color.y, pointLight.color.z);
	glUniform1f(glGetUniformLocation(rayshader, "u_light.power"), pointLight.power);
	glUniform1f(glGetUniformLocation(rayshader, "u_light.reach"), pointLight.reach);


	// plane material stuff
	glUniform3fv(glGetUniformLocation(rayshader, "u_planeMaterial.albedo"), 1, value_ptr(planeMaterial.albedo));
	glUniform3fv(glGetUniformLocation(rayshader, "u_planeMaterial.specular"), 1, value_ptr(planeMaterial.specular));
	glUniform3fv(glGetUniformLocation(rayshader, "u_planeMaterial.emission"), 1, value_ptr(planeMaterial.emission));
	glUniform1f(glGetUniformLocation(rayshader, "u_planeMaterial.emissionStrength"), planeMaterial.emissionStrength);
	glUniform1f(glGetUniformLocation(rayshader, "u_planeMaterial.roughness"), planeMaterial.roughness);
	glUniform1f(glGetUniformLocation(rayshader, "u_planeMaterial.specularHighlight"), planeMaterial.specularHighlight);
	glUniform1f(glGetUniformLocation(rayshader, "u_planeMaterial.specularExponent"), planeMaterial.specularExponent);

	// water material stuff
	glUniform3fv(glGetUniformLocation(rayshader, "u_waterMaterial.albedo"), 1, value_ptr(waterMaterial.albedo));
	glUniform3fv(glGetUniformLocation(rayshader, "u_waterMaterial.specular"), 1, value_ptr(waterMaterial.specular));
	glUniform3fv(glGetUniformLocation(rayshader, "u_waterMaterial.emission"), 1, value_ptr(waterMaterial.emission));
	glUniform1f(glGetUniformLocation(rayshader, "u_waterMaterial.emissionStrength"), waterMaterial.emissionStrength);
	glUniform1f(glGetUniformLocation(rayshader, "u_waterMaterial.roughness"), waterMaterial.roughness);
	glUniform1f(glGetUniformLocation(rayshader, "u_waterMaterial.specularHighlight"), waterMaterial.specularHighlight);
	glUniform1f(glGetUniformLocation(rayshader, "u_waterMaterial.specularExponent"), waterMaterial.specularExponent);

	// face buffer
	//GLuint ssbo;
	//glGenBuffers(1, &ssbo);
	//glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo);
	//glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(glm::vec3) * plane.vertices.size(), &plane.vertices[0], GL_STATIC_DRAW);
	//glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0); // Unbind

	plane.mesh.draw(); // draw
	// do the same thing for UVs but bind it to location=2 - the data is treated in lots of 2 (2 floats = vec2)
	glEnableVertexAttribArray(3);
	glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(vec3), value_ptr(pointLight.position));

	if (refreshRequired) {
		accumulatedPasses = 0;
		refreshRequired = false;
	}
	glUniform1i(glGetUniformLocation(rayshader, "u_accumulatedPasses"), accumulatedPasses);
	// DRAW OBJECTS

	//plane.mesh.draw(); // draw
}



/////////////////////////////////////// camera controls

void Application::cursorPosCallback(double xpos, double ypos) {
	if (m_leftMouseDown) {
		vec2 whsize = m_windowsize / 2.0f;

		// clamp the pitch to [-pi/2, pi/2]
		m_pitch += float(acos(glm::clamp((m_mousePosition.y - whsize.y) / whsize.y, -1.0f, 1.0f))
			- acos(glm::clamp((float(ypos) - whsize.y) / whsize.y, -1.0f, 1.0f)));
		m_pitch = float(glm::clamp(m_pitch, -pi<float>() / 2, pi<float>() / 2));

		// wrap the yaw to [-pi, pi]
		m_yaw += float(acos(glm::clamp((m_mousePosition.x - whsize.x) / whsize.x, -1.0f, 1.0f))
			- acos(glm::clamp((float(xpos) - whsize.x) / whsize.x, -1.0f, 1.0f)));
		if (m_yaw > pi<float>()) m_yaw -= float(2 * pi<float>());
		else if (m_yaw < -pi<float>()) m_yaw += float(2 * pi<float>());
	}

	// updated mouse position
	m_mousePosition = vec2(xpos, ypos);
}



void Application::mouseButtonCallback(int button, int action, int mods) {
	(void)mods; // currently un-used

	// capture is left-mouse down
	if (button == GLFW_MOUSE_BUTTON_LEFT)
		m_leftMouseDown = (action == GLFW_PRESS); // only other option is GLFW_RELEASE
}


void Application::scrollCallback(double xoffset, double yoffset) {
	(void)xoffset; // currently un-used
	m_distance *= pow(1.1f, -yoffset);
}


void Application::keyCallback(int key, int scancode, int action, int mods) {
	if (key == 265) {
		vy--;
	}
	else if (key == 264) {
		vy++;
	}
	else if (key == 263) {
		vx++;
	}
	else if (key == 262) {
		vx--;
	}
}


void Application::charCallback(unsigned int c) {
	(void)c; // currently un-used
}


// Camera stuff 

/// <----------------------------- >


void Application::calculateSpline(vec3& p0, vec3& p1, vec3& p2, vec3& p3, std::vector<glm::vec3>& intPoints, std::vector<glm::vec3>& sliderPoints) {
	float tension = 0;
	float alpha = 0.5;
	float t0 = 0.0f;
	float t1 = t0 + pow(distance(p0, p1), alpha);
	float t2 = t1 + pow(distance(p1, p2), alpha);
	float t3 = t2 + pow(distance(p2, p3), alpha);

	vec3 m1 = (1.0f - tension) * (t2 - t1) *
		((p1 - p0) / (t1 - t0) - (p2 - p0) / (t2 - t0) + (p2 - p1) / (t2 - t1));
	vec3 m2 = (1.0f - tension) * (t2 - t1) *
		((p2 - p1) / (t2 - t1) - (p3 - p1) / (t3 - t1) + (p3 - p2) / (t3 - t2));

	splineSeg segment;
	segment.a = 2.0f * (p1 - p2) + m1 + m2;
	segment.b = -3.0f * (p1 - p2) - m1 - m1 - m2;
	segment.c = m1;
	segment.d = p1;

	for (float t = 0; t < 1; t = t + 0.1) {
		vec3 point = segment.a * t * t * t +
			segment.b * t * t +
			segment.c * t +
			segment.d;

		intPoints.push_back(point);
	}

	for (float t = 0; t < 1; t = t + 0.01) {
		vec3 point = segment.a * t * t * t +
			segment.b * t * t +
			segment.c * t +
			segment.d;

		sliderPoints.push_back(point);
	}
}


void Application::plotPoints() {
	//<--------------------------------------- CAMERA ANIMATION AND SETTINGS--------------------------------------------------->

	//CONTROL POINTS

	//Start
	controlPts.push_back(vec3(120, 4, 100));
	controlPts.push_back(vec3(120, 4.5, 120));
	controlPts.push_back(vec3(90, 4.7, 80));
	controlPts.push_back(vec3(90, 4.3, 50));
	controlPts.push_back(vec3(90, 4.2, 0));

	controlPts.push_back(vec3(120, 4, -30));
	controlPts.push_back(vec3(120, 4.5, -60));
	controlPts.push_back(vec3(90, 4.7, -90));
	controlPts.push_back(vec3(90, 4.3, -120));
	controlPts.push_back(vec3(90, 4.2, -120));
	///Curve
	controlPts.push_back(vec3(70, 4, -120));
	controlPts.push_back(vec3(50, 4.5, -120));
	controlPts.push_back(vec3(30, 4.7, -120));
	controlPts.push_back(vec3(10, 4.3, -120));
	controlPts.push_back(vec3(10, 4.2, -120));

	controlPts.push_back(vec3(10, 4.3, 0));
	controlPts.push_back(vec3(10, 4.2, 0));

	controlPts.push_back(vec3(0, 4.3, 100));
	controlPts.push_back(vec3(0, 4.2, 110));

	controlPts.push_back(vec3(-70, 5, 120));
	controlPts.push_back(vec3(-70, 5.8, 120));

	controlPts.push_back(vec3(-110, 5, 120));
	controlPts.push_back(vec3(-110, 5.8, 120));
	//5
	controlPts.push_back(vec3(-110, 5, -120));
	controlPts.push_back(vec3(-110, 5.8, -120));

	controlPts.push_back(vec3(0, 5, -120));
	controlPts.push_back(vec3(0, 5.8, -120));

	controlPts.push_back(vec3(50, 5, -120));
	controlPts.push_back(vec3(50, 5.8, -120));

	controlPts.push_back(vec3(50, 5, -50));
	controlPts.push_back(vec3(50, 5.4, -45));

	controlPts.push_back(vec3(50, 5, 50));
	controlPts.push_back(vec3(50, 5.4, 45));
	//10
	controlPts.push_back(vec3(10, 0, 50));
	controlPts.push_back(vec3(-0, 0.4, 45));

	controlPts.push_back(vec3(-40, -3, 50));
	controlPts.push_back(vec3(-100, -3.4, 45));

	controlPts.push_back(vec3(-110, -6, 50));
	controlPts.push_back(vec3(-120, -6.4, 45));

	controlPts.push_back(vec3(-130, -8, 150));
	controlPts.push_back(vec3(-140, -8.4, 165));

	controlPts.push_back(vec3(-144, -10.9, 168));
	controlPts.push_back(vec3(-139, -10.4, 174));
	//15
	controlPts.push_back(vec3(-114, -5, 173));
	controlPts.push_back(vec3(-90, -5.9, 171.4));

	controlPts.push_back(vec3(-40, -2.7, 168));
	controlPts.push_back(vec3(0, -2.4, 164));

	controlPts.push_back(vec3(100, 0.7, 160));
	controlPts.push_back(vec3(130, 0.4, 158));

	controlPts.push_back(vec3(145, 1.7, 135));
	controlPts.push_back(vec3(150, 1.4, 120));

	controlPts.push_back(vec3(150, 3.7, 110));
	controlPts.push_back(vec3(150, 3.4, 100));
	//20
	controlPts.push_back(vec3(150, 3.7, 0));
	controlPts.push_back(vec3(150, 3.4, -60));

	controlPts.push_back(vec3(150, 3.7, -120));
	controlPts.push_back(vec3(150, 3.4, -140));

	controlPts.push_back(vec3(135, 3.7, -140));
	controlPts.push_back(vec3(120, 3.4, -140));

	controlPts.push_back(vec3(0, 3.7, -140));
	controlPts.push_back(vec3(-100, 3.4, -140));

	controlPts.push_back(vec3(-140, 0.7, -140));
	controlPts.push_back(vec3(-160, 0.4, -140));
	//25
	controlPts.push_back(vec3(-155, -1.7, -125));
	controlPts.push_back(vec3(-150, -1.4, -110));

	controlPts.push_back(vec3(-145, -2.5, -95));
	controlPts.push_back(vec3(-140, -5.7, -80));

	controlPts.push_back(vec3(-130, -10, -65));
	controlPts.push_back(vec3(-120, -17, -50));

	controlPts.push_back(vec3(0, -23, 0));
	controlPts.push_back(vec3(10, -23.4, 80));

	controlPts.push_back(vec3(20, -24.7, 100));
	controlPts.push_back(vec3(50, -25.5, 130));
	//30
	controlPts.push_back(vec3(80, -24.7, 160));
	controlPts.push_back(vec3(110, -25.5, 140));

	controlPts.push_back(vec3(105, -21.7, 100));
	controlPts.push_back(vec3(85, -21.5, 50));

	controlPts.push_back(vec3(65, -19.7, 0));
	controlPts.push_back(vec3(45, -18.5, -30));

	controlPts.push_back(vec3(25, -17.7, -60));
	controlPts.push_back(vec3(5, -16.5, -90));

	controlPts.push_back(vec3(-15, -17.7, -100));
	controlPts.push_back(vec3(-65, -16.5, -120));
	//35
	controlPts.push_back(vec3(-125, -17.7, -135));
	controlPts.push_back(vec3(-205, -16.5, -150));





	//<-------------------------------------------------------------------------------------------------------------------->

		//CAM-CONTROL POINTS


		//Start
	camControlPts.push_back(vec3(100, 0.2, 40));
	camControlPts.push_back(vec3(103, 0.5, 30));
	camControlPts.push_back(vec3(106, 0.6, 20));
	camControlPts.push_back(vec3(109, 0.77, 10));
	camControlPts.push_back(vec3(112, 0.9, 0));

	camControlPts.push_back(vec3(115, 1.4, -10));
	camControlPts.push_back(vec3(118, 1.5, -20));
	camControlPts.push_back(vec3(115, 1, -30));
	camControlPts.push_back(vec3(112, 1.6, -40));
	camControlPts.push_back(vec3(109, 1.5, -50));
	//Curve
	camControlPts.push_back(vec3(106, 2, -60));
	camControlPts.push_back(vec3(103, 3.5, -70));
	camControlPts.push_back(vec3(100, 4, -72));
	camControlPts.push_back(vec3(99.8, 4.6, -69));
	camControlPts.push_back(vec3(99.4, 5.5, -65));
	//1st
	camControlPts.push_back(vec3(98, 6.6, -64.5));
	camControlPts.push_back(vec3(97.8, 6.5, -63.9));

	camControlPts.push_back(vec3(98, 8.6, -63));
	camControlPts.push_back(vec3(98.4, 8.5, -62.8));

	camControlPts.push_back(vec3(99, 10.6, -62.35));
	camControlPts.push_back(vec3(99.7, 10.5, -61.2));

	camControlPts.push_back(vec3(100, 14.6, -60));
	camControlPts.push_back(vec3(105, 14.5, -55));

	camControlPts.push_back(vec3(110, 14.6, -50));
	camControlPts.push_back(vec3(114.7, 14.5, -45));
	//5
	camControlPts.push_back(vec3(117.95, 14.6, -40));
	camControlPts.push_back(vec3(121, 14.5, -35));

	camControlPts.push_back(vec3(126.4, 14.6, -30));
	camControlPts.push_back(vec3(130, 14.5, -25));

	camControlPts.push_back(vec3(145, 14.6, -20));
	camControlPts.push_back(vec3(155, 14.5, -15));

	camControlPts.push_back(vec3(140, 7.6, -10));
	camControlPts.push_back(vec3(125, 7.5, 0));

	camControlPts.push_back(vec3(107, 3.6, 10));
	camControlPts.push_back(vec3(90, 3.5, 20));
	//10
	camControlPts.push_back(vec3(65, 0.6, 20.6));
	camControlPts.push_back(vec3(57, 0.5, 23.9));

	camControlPts.push_back(vec3(44, -3.6, 26.4));
	camControlPts.push_back(vec3(22, -3.5, 29));

	camControlPts.push_back(vec3(0, -6.2, 33));
	camControlPts.push_back(vec3(-1.5, -6.2, 38.5));

	camControlPts.push_back(vec3(-3, -8.37, 43));
	camControlPts.push_back(vec3(-3.3, -8.4, 48));
	//15
	camControlPts.push_back(vec3(-4.7, -10.9, 66));
	camControlPts.push_back(vec3(-5.4, -10.97, 74));

	camControlPts.push_back(vec3(-6, -5, 66));
	camControlPts.push_back(vec3(-14, -5, 58));

	camControlPts.push_back(vec3(-22, -3.6, 50));
	camControlPts.push_back(vec3(-30, -3.67, 42));

	camControlPts.push_back(vec3(-38, 0, 30));
	camControlPts.push_back(vec3(-46, 0.5, 15));

	camControlPts.push_back(vec3(-54, 5, -25));
	camControlPts.push_back(vec3(-62, 4.5, -55));
	//20
	camControlPts.push_back(vec3(-70, 4.4, -59));
	camControlPts.push_back(vec3(-81, 4.35, -57));

	camControlPts.push_back(vec3(-84, 4.2, -53));
	camControlPts.push_back(vec3(-87, 3.15, -47));

	camControlPts.push_back(vec3(-89, 3.10, -43));
	camControlPts.push_back(vec3(-90, 3, -40));

	camControlPts.push_back(vec3(-93, 2.9, -39));
	camControlPts.push_back(vec3(-97, 2.67, -37));

	camControlPts.push_back(vec3(-100, 0.22, -33));
	camControlPts.push_back(vec3(-103, 0.03, -30));
	//25
	camControlPts.push_back(vec3(-107, -1.22, -28));
	camControlPts.push_back(vec3(-110, -1.03, -26));

	camControlPts.push_back(vec3(-113, 2.1, -16));
	camControlPts.push_back(vec3(-117, -5.5, -6));

	camControlPts.push_back(vec3(-120, -10.5, 0));
	camControlPts.push_back(vec3(-123, -17.5, 6));

	camControlPts.push_back(vec3(-127, -20.3, 16));
	camControlPts.push_back(vec3(-130, -20.1, 26));

	camControlPts.push_back(vec3(-133, -20.7, 36));
	camControlPts.push_back(vec3(-136, -20.5, 46));
	//30
	camControlPts.push_back(vec3(-137, -18.7, 56));
	camControlPts.push_back(vec3(-138, -18.5, 66));

	camControlPts.push_back(vec3(-139, -19.7, 76));
	camControlPts.push_back(vec3(-140, -19.5, 86));

	camControlPts.push_back(vec3(-141, -18.7, 96));
	camControlPts.push_back(vec3(-142, -16.5, 106));

	camControlPts.push_back(vec3(-143, -14.7, 116));
	camControlPts.push_back(vec3(-144, -12.5, 125));

	camControlPts.push_back(vec3(-145, -10.7, 135));
	camControlPts.push_back(vec3(-146, -8.5, 145));
	//35
	camControlPts.push_back(vec3(-147, -6.7, 155));
	camControlPts.push_back(vec3(-148, -4.5, 165));
	//<------------------------------------------------------------>






}