#version 330 core

// uniform data
uniform mat4 uProjectionMatrix;
uniform mat4 uModelViewMatrix;
uniform vec3 uColor;

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

	// output to the framebuffer
	fb_color = vec4(color, 1);
}