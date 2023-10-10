#version 330 core

// uniform data
uniform mat4 uProjectionMatrix;
uniform mat4 uModelViewMatrix;
uniform vec3 uColor;
uniform vec3 u_lightTranslation;
uniform float u_wavetime;
uniform float u_time;

// mesh data
layout(location = 0) in vec3 aPosition;
layout(location = 1) in vec3 aNormal;
layout(location = 2) in vec2 aTexCoord;

layout(location = 3) in vec3 lightPosition;

layout(location = 4) in float aBelowThreshold;

 //layout(std430, binding = 0) buffer VertexBuffer {
  //     vec3 vertices[];
  // };


// model data (this must match the input of the vertex shader)
out VertexData {
	vec3 position;
	vec3 normal;
	vec2 textureCoord;
} v_out;

out LightData{
	vec3 position;
}vl_out;

out float threshold;

///////////////// Terrain Methods - KAHU

// A basic noise function to generate some randomness
float noise(vec2 st) {
    return fract(sin(dot(st.xy, vec2(12.9898,78.233))) * 43758.5453);
}

float waveFunction(vec3 position, vec2 direction, float waveLength, float waveSpeed, float waveHeight) {
    float dist = dot(position.xz, direction);
    return waveHeight * sin(u_time * waveSpeed + dist * 2.0 * 3.14159 / waveLength);
}

float dWaveFunctiondx(vec3 position, vec2 direction, float waveLength, float waveSpeed, float waveHeight) {
    float distDerivativeX = direction.x;
    return waveHeight * cos(u_time * waveSpeed + dot(position.xz, direction) * 2.0 * 3.14159 / waveLength) * distDerivativeX;
}

float dWaveFunctiondz(vec3 position, vec2 direction, float waveLength, float waveSpeed, float waveHeight) {
    float distDerivativeZ = direction.y;
    return waveHeight * cos(u_time * waveSpeed + dot(position.xz, direction) * 2.0 * 3.14159 / waveLength) * distDerivativeZ;
}




///////////////////

void main() {
	// transform vertex data to viewspace
	vec3 pos = aPosition;
	
    if (aBelowThreshold > 0.5) {
    	threshold = 1;
        vec2 st = aPosition.xz * 0.1;
        vec2 waveDirection = normalize(vec2(noise(st), noise(st + 10.0)));
        float waveLength = 0.1;
        float waveSpeed = 3;
        float waveHeight = 0.1;

        pos.y += waveFunction(pos, waveDirection, waveLength, waveSpeed, waveHeight);
        
        // Calculate the tangent vectors
        vec3 tangentX = vec3(1.0, dWaveFunctiondx(pos, waveDirection, waveLength, waveSpeed, waveHeight), 0.0);
        vec3 tangentZ = vec3(0.0, dWaveFunctiondz(pos, waveDirection, waveLength, waveSpeed, waveHeight), 1.0);
        
        //tangentX = (uModelViewMatrix * vec4(tangentX, 0)).xyz;
        //tangentZ = (uModelViewMatrix * vec4(tangentZ, 0)).xyz;

         // Compute the normal by crossing the two tangents
        //v_out.normal = normalize(uModelViewMatrix * vec4(aNormal, 0)).xyz;
        v_out.normal = normalize(cross(tangentZ, tangentX));
        
        // Transform the normal by the uModelViewMatrix
       
    } else {
    	threshold = 0;
        v_out.normal = normalize((uModelViewMatrix * vec4(aNormal, 0)).xyz);
    }
    
	
    v_out.position =  (uModelViewMatrix * vec4(pos, 1)).xyz;
	v_out.textureCoord = aTexCoord;

	// light position 
	vl_out.position = (uModelViewMatrix * vec4(lightPosition+u_lightTranslation, 1)).xyz;


	// set the screenspace position (needed for converting to fragment data)
	gl_Position = uProjectionMatrix * uModelViewMatrix * vec4(pos, 1);
}