#version 330 core

// uniform data
uniform mat4 uProjectionMatrix;
uniform mat4 uModelViewMatrix;
uniform vec3 uColor;

// mesh data
layout(location = 0) in vec3 aPosition;
layout(location = 1) in vec3 aNormal;
layout(location = 2) in vec2 aTexCoord;

// model data (this must match the input of the vertex shader)
out VertexData {
	vec3 position;
	vec3 normal;
	vec2 textureCoord;
    flat int isAboveThreshold; // Use flat qualifier to ensure consistent value
} v_out;

void main() {
    // Apply threshold to vertices in object space (local coordinates)
    float threshold = -1.0; // Adjust the threshold value as needed
    if (aPosition.y >= threshold) {
        v_out.position = (uModelViewMatrix * vec4(aPosition, 1)).xyz;
        v_out.normal = normalize((uModelViewMatrix * vec4(aNormal, 0)).xyz);
        v_out.textureCoord = aTexCoord;
        v_out.isAboveThreshold = 0; // Mark vertices below threshold
    } else {
        v_out.position = (uModelViewMatrix * vec4(aPosition, 1)).xyz;
        v_out.normal = normalize((uModelViewMatrix * vec4(aNormal, 0)).xyz);
        v_out.textureCoord = aTexCoord;
        v_out.isAboveThreshold = 1; // Mark vertices above threshold
    }

    // Set the screenspace position (needed for converting to fragment data)
    gl_Position = uProjectionMatrix * uModelViewMatrix * vec4(aPosition, 1);
}