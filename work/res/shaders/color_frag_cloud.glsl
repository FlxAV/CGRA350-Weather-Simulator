#version 330 core

// uniform data
uniform mat4 uProjectionMatrix;
uniform mat4 uModelViewMatrix;
uniform vec3 uColor;
uniform float u_Brightness;

uniform float uTransparency;

// Fog Constants and Uniforms
#define FOG_START 0
#define FOG_END 500
#define FOG_STRENGTH 1
#define FOG_COLOR vec3(0.7, 0.7, 0.9)

// viewspace data (this must match the output of the fragment shader)
in VertexData {
	vec3 position;
	vec3 normal;
	vec2 textureCoord;
    flat int isAboveThreshold; // Use the attribute from the vertex shader
} f_in;

// framebuffer output
out vec4 fb_color;

void main() {
	// calculate lighting (hack)
	vec3 eye = normalize(-f_in.position);
	float light = abs(dot(normalize(f_in.normal), eye));
	vec3 color;

	// Check the attribute to determine the color
	if (f_in.isAboveThreshold == 1) {
		color = mix(uColor / 4, uColor, light);
	} else {

		//color = vec3(1.0, 0.0, 0.0); // Red color for vertices below or on the threshold
		discard;
	}

 // Calculate the distance to the camera for the fog effect
    float distanceToCamera = length(f_in.position - eye);

    // Fog calculation
    if (distanceToCamera > FOG_START && distanceToCamera < FOG_END) {
        float normalizedDistance = (distanceToCamera - FOG_START) / (FOG_END - FOG_START);
        float fogFactor = normalizedDistance * FOG_STRENGTH;

        color = mix(color, FOG_COLOR * u_Brightness, fogFactor);
    } else if (distanceToCamera >= FOG_END) {
        color = FOG_COLOR ;
    }

    // Output to the framebuffer
    fb_color = vec4(color * u_Brightness , uTransparency);
}
