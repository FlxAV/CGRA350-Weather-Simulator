#version 330 core

// uniform data
uniform mat4 uProjectionMatrix;
uniform mat4 uModelViewMatrix;
uniform vec3 uColor;
uniform float uTime;


// mesh data
layout(location = 0) in vec3 aPosition;
layout(location = 1) in vec3 aNormal;
layout(location = 2) in vec2 aTexCoord;
layout(location = 4) in float aBelowThreshold;

// model data (this must match the input of the vertex shader)
out VertexData {
    vec3 position;
    vec3 normal;
    vec2 textureCoord;
} v_out;

// A basic noise function to generate some randomness
float noise(vec2 st) {
    return fract(sin(dot(st.xy, vec2(12.9898,78.233))) * 43758.5453);
}

float waveFunction(vec3 position, vec2 direction, float waveLength, float waveSpeed, float waveHeight) {
    float dist = dot(position.xz, direction);
    return waveHeight * sin(uTime * waveSpeed + dist * 2.0 * 3.14159 / waveLength);
}

float dWaveFunctiondx(vec3 position, vec2 direction, float waveLength, float waveSpeed, float waveHeight) {
    float distDerivativeX = direction.x;
    return waveHeight * cos(uTime * waveSpeed + dot(position.xz, direction) * 2.0 * 3.14159 / waveLength) * distDerivativeX;
}

float dWaveFunctiondz(vec3 position, vec2 direction, float waveLength, float waveSpeed, float waveHeight) {
    float distDerivativeZ = direction.y;
    return waveHeight * cos(uTime * waveSpeed + dot(position.xz, direction) * 2.0 * 3.14159 / waveLength) * distDerivativeZ;
}

void main() {
    vec3 pos = aPosition;

    if (aBelowThreshold > 0.5) {
        vec2 st = aPosition.xz * 0.1;
        vec2 waveDirection = normalize(vec2(noise(st), noise(st + 10.0)));
        float waveLength = 0.1;
        float waveSpeed = 3;
        float waveHeight = 0.1;

        pos.y += waveFunction(pos, waveDirection, waveLength, waveSpeed, waveHeight);
        
        // Calculate the tangent vectors
        vec3 tangentX = vec3(1.0, dWaveFunctiondx(pos, waveDirection, waveLength, waveSpeed, waveHeight), 0.0);
        vec3 tangentZ = vec3(0.0, dWaveFunctiondz(pos, waveDirection, waveLength, waveSpeed, waveHeight), 1.0);
        
        // Compute the normal by crossing the two tangents
        v_out.normal = normalize(cross(tangentZ, tangentX));
    } else {
        v_out.normal = aNormal;
    }
    
    v_out.position = vec3(uModelViewMatrix * vec4(pos, 1.0));
    v_out.textureCoord = aTexCoord;

    gl_Position = uProjectionMatrix * uModelViewMatrix * vec4(pos, 1.0);
}