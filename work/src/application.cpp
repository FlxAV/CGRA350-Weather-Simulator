
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


using namespace std;
using namespace cgra;
using namespace glm;


#define M_PI 3.1415926535897932384626433832795


void basic_model::draw(const glm::mat4 &view, const glm::mat4 proj) {
	mat4 modelview = view * modelTransform;
	
	glUseProgram(shader); // load shader and variables
	glUniformMatrix4fv(glGetUniformLocation(shader, "uProjectionMatrix"), 1, false, value_ptr(proj));
	glUniformMatrix4fv(glGetUniformLocation(shader, "uModelViewMatrix"), 1, false, value_ptr(modelview));
	glUniform3fv(glGetUniformLocation(shader, "uColor"), 1, value_ptr(color));

	mesh.draw(); // draw
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

	int size = 50;
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
			std::cout << u << "\n";
			std::cout << v << "\n";
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
	m_plane.color = { 1,1,1 };

	return mb;
}
 
//////////////////////// SCENE MODELS ////////////////////////
Object::Object() = default;
Object::Object(int ntype, const glm::vec3& nposition, const glm::vec3 & nscale, Material &nmaterial) {
	type = ntype; 
	position = nposition; 
	scale = nscale; 
	material = nmaterial;

}




// Material
Material::Material() = default;
Material::Material(const glm::vec3&nalbedo, const glm::vec3&nspecular, const glm::vec3&nemission, float emissionStrength, float roughness, float specularHighlight, float specularExponent) {
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
PointLight::PointLight(const glm::vec3& nposition, const glm::vec3& ndirection, float nradius, const glm::vec3&ncolor, float power, float reach) {
	position = nposition;
	direction = ndirection;
	radius = nradius;
	color = ncolor;

	this->power = power;
	this->reach = reach;
}


//////////////////////// Main application class ////////////////////////

Application::Application(GLFWwindow *window) : m_window(window) {
	
	//shader_builder sb;
	//sb.set_shader(GL_VERTEX_SHADER, CGRA_SRCDIR + std::string("//res//shaders//color_vert.glsl"));
	//sb.set_shader(GL_FRAGMENT_SHADER, CGRA_SRCDIR + std::string("//res//shaders//color_frag.glsl"));
	//GLuint shader = sb.build();

	////m_model.shader = shader;
	////m_model.mesh = load_wavefront_data(CGRA_SRCDIR + std::string("/res//assets//teapot.obj")).build();
	////m_model.color = vec3(1, 0, 0);



	//m_model.shader = shader;
	//m_model.mesh = drawUVSphere().build();

	m_rayplane.mesh = drawPlane().build();
	planeMaterial = Material(vec3(0.4F, 0.4F, 0.4F), vec3(0.75F, 0.75F, 0.75F), vec3(0.0F, 0.0F, 0.0F), 0.0F, 0.0F, 0.0F, 0.0F);
	m_rayplane.material = planeMaterial;
	//placeBasicScene();
	//recompileShader();

	// Position , radius , color , power , reach
	pointLight = PointLight(vec3(1.0,5.0,1.0),vec3(-1,-1,0),6, vec3(0.7), 15, 15);
	buildRayBasicShader();

	//pointLight.position += vec3(5, 0, 0);
	glUseProgram(rayshader);
	glUniform1i(glGetUniformLocation(rayshader, "u_screenTexture"), 0);
	//placeBasicScene();
	//recompileShader();

	//glUniform1i(glGetUniformLocation(rayshader, "u_screenTexture"), 0);

	
	int charles = 5000000;
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

	// projection matrix
	proj = perspective(1.f, float(width) / height, 0.1f, 1000.f);

	// view matrix
	view = translate(mat4(1), vec3(vx, vy, -m_distance))
		* rotate(mat4(1), m_pitch, vec3(1, 0, 0))
		* rotate(mat4(1), m_yaw,   vec3(0, 1, 0));


	camPos = vec3(0,0,0);


	// helpful draw options
	if (m_show_grid) drawGrid(view, proj);
	if (m_show_axis) drawAxis(view, proj);
	glPolygonMode(GL_FRONT_AND_BACK, (m_showWireframe) ? GL_LINE : GL_FILL);



	// draw the model
	//m_model.draw(view, proj);

	//m_model.draw(view, proj);
	//m_plane.draw(view, proj);
	drawBasicScene(view, proj, preTime);

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

	ImGui::SliderInt("framePasses", &framePasses, 0, 100);
	ImGui::SliderInt("shadowResolution", &shadowResolution, 0, 100);
	ImGui::SliderInt("lightBounces", &lightBounces, 0, 15);
	ImGui::SliderFloat("light power", &pointLight.power, 0, 1500);
	ImGui::SliderFloat3("light position", value_ptr(lightTranslate), -10, 10);

	// helpful drawing options
	ImGui::Checkbox("Show axis", &m_show_axis);
	ImGui::SameLine();
	ImGui::Checkbox("Show grid", &m_show_grid);
	ImGui::Checkbox("Wireframe", &m_showWireframe);
	ImGui::SameLine();
	if (ImGui::Button("Screenshot")) rgba_image::screenshot(true);

	
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
void Application::drawBasicScene(const glm::mat4& view, const glm::mat4 proj,double time ) {

	// DRAW PLANE
	vec3 pos = pointLight.position;
	mat4 modelview = view * m_rayplane.modelTransform;
	
	cout << "old position: " << pointLight.position.x << "--" << pointLight.position.y << "--" << pointLight.position.z << "\n";

	mat4 transform = view * mat4(1);
	vec3 v(modelview * vec4(pointLight.position, 1));

	//pointLight.position = v;

	//cout << "new position: " << v.x << "--" << v.y << "--" << v.z << "--" << v.w << "\n";
	//pointLight.position = v;

	preView = modelview;

	//pointLight.position = proj * modelview * vec4(pointLight.position, 1);
	// plane vertex stuff
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
	glUniform3fv(glGetUniformLocation(rayshader, "u_cameraPosition"), 1, value_ptr(camPos));
	glUniform3fv(glGetUniformLocation(rayshader, "u_lightTranslation"), 1, value_ptr(lightTranslate));

	// light stuff
	glUniform3fv(glGetUniformLocation(rayshader, "u_light.position"), 1,value_ptr(pointLight.position));
	glUniform3fv(glGetUniformLocation(rayshader, "u_light.position"), 1, value_ptr(pointLight.direction));
	glUniform1f(glGetUniformLocation(rayshader, "u_light.radius"), pointLight.radius);
	glUniform3f(glGetUniformLocation(rayshader, "u_light.color"), pointLight.color.x, pointLight.color.y, pointLight.color.z);
	glUniform1f(glGetUniformLocation(rayshader, "u_light.power"), pointLight.power);
	glUniform1f(glGetUniformLocation(rayshader, "u_lights.reach"), pointLight.reach);


	// plane  material stuff
	glUniform3f(glGetUniformLocation(rayshader, "u_planeMaterial.albedo"), m_rayplane.material.albedo.x, m_rayplane.material.albedo.y, m_rayplane.material.albedo.z);
	glUniform3f(glGetUniformLocation(rayshader, "u_planeMaterial.specular"), m_rayplane.material.specular.x, m_rayplane.material.specular.y, m_rayplane.material.specular.z);
	glUniform3f(glGetUniformLocation(rayshader, "u_planeMaterial.emission"), m_rayplane.material.emission.x, m_rayplane.material.emission.y, m_rayplane.material.emission.z);
	glUniform1f(glGetUniformLocation(rayshader, "u_planeMaterial.emissionStrength"), m_rayplane.material.emissionStrength);
	glUniform1f(glGetUniformLocation(rayshader, "u_planeMaterial.roughness"), m_rayplane.material.roughness);
	glUniform1f(glGetUniformLocation(rayshader, "u_planeMaterial.specularHighlight"), m_rayplane.material.specularHighlight);
	glUniform1f(glGetUniformLocation(rayshader, "u_planeMaterial.specularExponent"), m_rayplane.material.specularExponent);

	m_rayplane.mesh.draw(); // draw


	// do the same thing for UVs but bind it to location=2 - the data is treated in lots of 2 (2 floats = vec2)
	glEnableVertexAttribArray(3);
	glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(vec3),value_ptr(pointLight.position));

	if (refreshRequired) {
		accumulatedPasses = 0;
		refreshRequired = false;
	}
	glUniform1i(glGetUniformLocation(rayshader, "u_accumulatedPasses"), accumulatedPasses);
	// DRAW OBJECTS
	//pointLight.position = pos;
}




// Draws ray tracer using the compiled shaders -> RAYTRACING SECTION
void Application::drawRayShader() {
	// TODO
}

// ?? draw  -> RAYTRACING SECTION







bool Application::sphereIntersection(glm::vec3 position, float radius, glm::vec3 rayOrigin, glm::vec3 rayDirection, float* hitDistance) {
	float t = glm::dot(position - rayOrigin, rayDirection);
	glm::vec3 p = rayOrigin + rayDirection * t;

	float y = glm::length(position - p);
	if (y < radius) {
		float x = sqrt(radius * radius - y * y);
		float t1 = t - x;
		if (t1 > 0) {
			*hitDistance = t1;
			return true;
		}

	}

	return false;
}

bool Application::boxIntersection(glm::vec3 position, glm::vec3 size, glm::vec3 rayOrigin, glm::vec3 rayDirection, float* hitDistance) {
	float t1 = -1000000000000.0;
	float t2 = 1000000000000.0;

	glm::vec3 boxMin = position - size / 2.0f;
	glm::vec3 boxMax = position + size / 2.0f;

	glm::vec3 t0s = (boxMin - rayOrigin) / rayDirection;
	glm::vec3 t1s = (boxMax - rayOrigin) / rayDirection;

	glm::vec3 tsmaller = min(t0s, t1s);
	glm::vec3 tbigger = max(t0s, t1s);

	t1 = std::max({ t1, tsmaller.x, tsmaller.y, tsmaller.z });
	t2 = std::min({ t2, tbigger.x, tbigger.y, tbigger.z });

	*hitDistance = t1;

	return t1 >= 0 && t1 <= t2;
}

bool Application::planeIntersection(glm::vec3 planeNormal, glm::vec3 planePoint, glm::vec3 rayOrigin, glm::vec3 rayDirection, float* hitDistance)
{
	float denom = glm::dot(planeNormal, rayDirection);
	if (abs(denom) > 0.0001) {
		glm::vec3 d = planePoint - rayOrigin;
		*hitDistance = glm::dot(d, planeNormal) / denom;
		return (*hitDistance >= 0.0001);
	}

	return false;
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


//
void Application::drawScene(const glm::mat4& view, const glm::mat4 proj) {

	double preTime = glfwGetTime();

	glUseProgram(rayshader);
	// scene uniform variables
	if (refreshRequired) {
		accumulatedPasses = 0;
		refreshRequired = false;
		glUniform1i(accumulatedPassesUniformLocation, accumulatedPasses); // If the shader receives a value of 0 for accumulatedPasses, it will discard the buffer and just output what it rendered on that frame.
	}

	glUniform1f(timeUniformLocation, (float)preTime);
	glUniform3f(camPosUniformLocation, 0, 0, 0);
	mat4 rotationMatrix = mat4(1);

	glUniformMatrix4fv(rotationMatrixUniformLocation, 1, GL_FALSE, glm::value_ptr(rotationMatrix));
	glUniform1f(aspectRatioUniformLocation, (float)width / height);


	glUniform1i(directOutPassUniformLocation, 0);
	accumulatedPasses += 1;

	glUniform1i(directOutPassUniformLocation, 1);
	glUniform1i(accumulatedPassesUniformLocation, accumulatedPasses);



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