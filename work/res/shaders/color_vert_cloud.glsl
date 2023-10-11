#version 330 core

// uniform data
uniform mat4 uProjectionMatrix;
uniform mat4 uModelViewMatrix;
uniform vec3 uColor;

// Define a uniform for the lateral shift and height variation
uniform float uTime;
uniform float uHeightVariation = 5.0;  // Increasing this will create larger height variations
uniform float uLateralShift;



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
    float threshold = -1.0; // Adjust the threshold value as needed
    vec3 pos = aPosition;
    pos.x += sin(uTime* (0.1 + uLateralShift)) * 0.1 ;

    if (aPosition.y >= threshold) {
        v_out.isAboveThreshold = 0; // Mark vertices above threshold
    } else {
        v_out.isAboveThreshold = 1; // Mark vertices below threshold

    }

 pos.y += (sin(uTime * 0.01 + (pos.x + pos.z) * 1.1) + cos(uTime * 0.015 + (pos.x - pos.z) * 1.3)) * (uHeightVariation * 1.5);




    v_out.position = vec3(uModelViewMatrix * vec4(pos, 1.0));
    v_out.normal = normalize((uModelViewMatrix * vec4(aNormal, 0)).xyz);
    v_out.textureCoord = aTexCoord;

    // Set the screenspace position (needed for converting to fragment data)
    gl_Position = uProjectionMatrix * uModelViewMatrix * vec4(pos, 1);
}
