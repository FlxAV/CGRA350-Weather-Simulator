#version 330 core

// uniform data
uniform vec3 uColor;

// viewspace data
in VertexData {
    vec3 position;
    vec3 normal;
    vec2 textureCoord;
} f_in;

// framebuffer output
out vec4 fb_color;

// Light properties
const vec3 lightDir = normalize(vec3(1.0, -1.0, 1.0)); // Direction of the light

void main() {
    vec3 norm = normalize(f_in.normal);

    // Calculate diffuse lighting
    float diff = max(dot(norm, -lightDir), 0.0);
    
    // Ambient lighting
    const float ambientStrength = 0.1f;
    vec3 ambient = ambientStrength * uColor;

    // Multiply our color by the diffuse lighting
    vec3 diffuse = uColor * diff;

    // Combine ambient and diffuse lighting
    vec3 litColor = ambient + diffuse;

    // Set the computed color to the framebuffer
    fb_color = vec4(litColor, 1.0);
}