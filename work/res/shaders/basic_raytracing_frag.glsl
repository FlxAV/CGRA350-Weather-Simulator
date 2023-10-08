#version 430 core

#define RENDER_DISTANCE 1000000
#define EPSILON 0.0001
#define PI 3.1415926538
#define OUTLINE_WIDTH 0.004
#define OUTLINE_COLOR vec4(1.0, 0.0, 1.0, 1.0)


// uniform data
uniform mat4 uProjectionMatrix;
uniform mat4 uModelViewMatrix;
uniform vec3 uColor;
uniform float u_aspectRatio;
uniform sampler2D u_screenTexture;
uniform int u_framePasses;
uniform float u_time;

uniform int u_lightBounces;
uniform int u_shadowResolution;
uniform int u_accumulatedPasses;
uniform vec3 u_cameraPosition;

//layout(std430, binding = 0) buffer VertexBuffer {
   // vec3 vertices[];
//};

// viewspace data (this must match the output of the fragment shader) -> plane data
in VertexData {
	vec3 position;
	vec3 normal;
	vec2 textureCoord;
} f_in;

in LightData{
	vec3 position;
}fl_in;

// framebuffer output
out vec4 fb_color;

////////////////////////////// STRUCTS


struct PointLight {
	vec3 position;
	vec3 direction;
	float radius;
	vec3 color;
	float power;
	float reach; // Only points within this distance of the light will be affected
};
uniform PointLight u_light;


struct Ray {
	vec3 origin;
	vec3 direction;
};

struct Material {
	vec3 albedo;
	vec3 specular;
	vec3 emission;
	float emissionStrength;
	float roughness;
	float specularHighlight;
	float specularExponent;
};
uniform Material u_planeMaterial;

struct SurfacePoint {
	vec3 position;
	vec3 normal;
	Material material;
};


/////////////////////////////// METHODS

float rand(vec2 co){
    return fract(sin(dot(co, vec2(12.9898, 78.233))) * 43758.5453);
}


bool planeIntersection(Ray ray, out float hitDistance) 
{ 
    
	vec3 norm = f_in.normal;
    float denom = dot(norm, ray.direction); 
    
    if (abs(denom) > EPSILON) { 
        vec3 d = f_in.position - ray.origin; 
        hitDistance = dot(d, norm) / denom; 
        return (hitDistance >= EPSILON); 
    } 
 
    return false; 
} 



bool terrainIntersection(Ray ray) {
    bool hit = false;
    float t, u, v;

   // for (int i = 0; i < vertices.length(); i += 3) {
        
        vec3 v0 = vec3(0);//vertices[i];
        vec3 v1 = vec3(0);//vertices[i+1];
        vec3 v2 = vec3(0);//vertices[i+2];


        vec3 edge1 = v1 - v0;
        vec3 edge2 = v2 - v0;
        vec3 h = cross(ray.direction, edge2);
        float a = dot(edge1, h);

        //if (a > -EPSILON && a < EPSILON) {
           // continue; // Ray is parallel to the triangle.
       // }

        float f = 1.0 / a;
        vec3 s = ray.origin - v0;
        u = f * dot(s, h);

        //if (u < 0.0 || u > 1.0) {
            //continue;
       // }

        vec3 q = cross(s, edge1);
        v = f * dot(ray.direction, q);

        //if (v < 0.0 || u + v > 1.0) {
           // continue;
        //}

        t = f * dot(edge2, q);

        if (t > EPSILON) {
         
            hit = true;
            // If you only want the first intersection, uncomment the following line:
            // break;
        }
    //}

    return hit;
}


// TODO: if coordinate blocks plane -> MOUNTAINS 
// Iterate through each triangle(face) -> check whether ray collides with face

bool terrainLightBlock(Ray lightRay){
	
    bool hit = false;
   // for (int i = 0; i < vertices.length(); i += 3) {
        vec3 v0 = vec3(0);//vertices[i];
        vec3 v1 = vec3(0);//vertices[i+1];
        vec3 v2 = vec3(0);//vertices[i+2];

        vec3 e1 = v1 - v0;
    	vec3 e2 = v2 - v0;
    	vec3 h = cross(lightRay.direction, e2);
    	float a = dot(e1, h);
    
    	if (a > -EPSILON && a < EPSILON)
        	return false;
    
    	float f = 1.0 / a;
    	vec3 s = lightRay.origin - v0;
    	float u = f * dot(s, h);
    
    	if (u < 0.0 || u > 1.0)
        	return false;
    
    	vec3 q = cross(s, e1);
    	float v = f * dot(lightRay.direction, q);
    
    	if (v < 0.0 || u + v > 1.0)
        	return false;


    	;
    //}

    return false;

}


bool raycast(Ray ray, out SurfacePoint hitPoint) {
	bool didHit = false;
	float minHitDist = RENDER_DISTANCE;

	float hitDist = distance(f_in.position, ray.origin);

	if (planeIntersection(ray, hitDist)) {
		didHit = true;
		if (hitDist < minHitDist) {
			minHitDist = hitDist;
			hitPoint.position = f_in.position;//ray.origin + ray.direction * minHitDist;
			hitPoint.normal = f_in.normal;
			hitPoint.material = u_planeMaterial;
		}
	}

	return didHit;
}

mat3x3 getTangentSpace(vec3 normal)
{
    // Choose a helper vector for the cross product
    vec3 helper = vec3(1, 0, 0);
    if (abs(normal.x) > 0.99)
        helper = vec3(0, 0, 1);

    // Generate vectors
    vec3 tangent = normalize(cross(normal, helper));
    vec3 binormal = normalize(cross(normal, tangent));
    return mat3x3(tangent, binormal, normal);
}

vec3 sampleHemisphere(vec3 normal, float alpha, vec2 seed)
{
    // Sample the hemisphere, where alpha determines the kind of the sampling
    float cosTheta = pow(rand(seed), 1.0 / (alpha + 1.0));
    float sinTheta = sqrt(1.0 - cosTheta * cosTheta);
    float phi = 2 * PI * rand(seed.yx);
    vec3 tangentSpaceDir = vec3(cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta);

    // Transform direction to world space
    return getTangentSpace(normal) * tangentSpaceDir;
}

//////////////PART 2: DIRECT ILLUMINATION

// Adds up the total light received directly from all light sources
vec3 computeDirectIllumination(SurfacePoint point, vec3 observerPos, float seed, PointLight light) {
	
	//light.position = (uModelViewMatrix * vec4(light.position, 1)).xyz;

	vec3 directIllumination = vec3(0);
	float lightDistance = length(light.position - point.position);
	
	if (lightDistance > light.reach) {return directIllumination;}

	float diffuse = clamp(dot(point.normal, normalize(light.position-point.position)), 0.0, 1.0);

	if (diffuse > EPSILON || point.material.roughness < 1.0) {
		// Shadow raycasting
		int shadowRays = int(u_shadowResolution/(lightDistance*lightDistance)+1); 
		// There must be a better way to find the right amount of shadow rays
		int shadowRayHits = 0;
		for (int i = 0; i<shadowRays; i++) {
			
			// Sample a point on the light sphere
			vec3 lightSurfacePoint = light.position + normalize(vec3(rand(vec2(i+seed, 1)+point.position.xy), rand(vec2(i+seed, 2)+point.position.yz), rand(vec2(i+seed, 3)+point.position.xz)));
			vec3 lightDir = normalize(lightSurfacePoint - point.position);
			vec3 rayOrigin = point.position + lightDir * EPSILON * 2.0;

			float maxRayLength = length(lightSurfacePoint - rayOrigin);
			Ray shadowRay = Ray(rayOrigin, lightDir);
			SurfacePoint SR_hit;
			if (raycast(shadowRay, SR_hit)) {
				if (length(SR_hit.position-rayOrigin) < maxRayLength) {
					shadowRayHits += 1;
				}
			}

		}
		// Diffuse 
		float attenuation = lightDistance * lightDistance;
		float shadowFactor = (1.0 - float(shadowRayHits) / shadowRays);
		shadowFactor = pow(shadowFactor, 2.0); 

		directIllumination += light.color * light.power* diffuse * point.material.albedo * shadowFactor / attenuation;
	
		// Specular highlight
		vec3 lightDir = normalize(f_in.position - point.position);

		vec3 reflectedLightDir = reflect(lightDir, point.normal);
		vec3 cameraDir = normalize(observerPos - point.position);
		directIllumination += point.material.specularHighlight * light.color * (light.power/(lightDistance*lightDistance)) * pow(max(dot(cameraDir, reflectedLightDir), 0.0), 1.0/max(point.material.specularExponent, EPSILON));

	}
	

	return directIllumination;
}


// SunlightIllumination
void sunlightIllumination(){
	
	

}


vec3 computeSceneColor(Ray cameraRay, float seed) {
	vec3 totalIllumination = vec3(0);
	vec3 rayOrigin = cameraRay.origin;
	vec3 rayDirection = cameraRay.direction;
	vec3 energy = vec3(1.0);

	for (int depth = 0; depth < u_lightBounces; depth++) {
		
		SurfacePoint hitPoint;
		if (raycast(Ray(rayOrigin, rayDirection), hitPoint)) {

			// Part one: Hit object's emission
			totalIllumination += energy * hitPoint.material.emission * hitPoint.material.emissionStrength;

			// Part two: Direct light (received directly from light sources)
			// IF statement to chceck whether some of the terrain is blocking
			PointLight light = u_light;
			light.position = fl_in.position;


			vec3 lightDir = normalize(f_in.position - light.position);
			Ray lightRay = Ray(light.position, lightDir);


			//if (!terrainIntersection(lightRay))

			totalIllumination += energy * computeDirectIllumination(hitPoint, rayOrigin, seed,light);


			// Part three: Indirect light (other objects)
			float specChance = dot(hitPoint.material.specular, vec3(1.0/3.0));
			float diffChance = dot(hitPoint.material.albedo, vec3(1.0/3.0));

			float sum = specChance + diffChance;
			specChance /= sum;
			diffChance /= sum;

			//specChance = 0;
			//diffChance = 0;

			// Roulette-select the ray's path
			float roulette = rand(hitPoint.position.zx+vec2(hitPoint.position.y)+vec2(seed, depth));

			if (roulette < specChance)
			{
				// Specular reflection
				float smoothness = 1.0-hitPoint.material.roughness;
				float alpha = pow(1000.0, smoothness*smoothness);
				if (smoothness == 1.0) {
					rayDirection = reflect(rayDirection, hitPoint.normal);
				} else {
					rayDirection = sampleHemisphere(reflect(rayDirection, hitPoint.normal), alpha, hitPoint.position.zx+vec2(hitPoint.position.y)+vec2(seed, depth));
				}
				rayOrigin = hitPoint.position + rayDirection * EPSILON;
				float f = (alpha + 2) / (alpha + 1);
				energy *= hitPoint.material.specular * clamp(dot(hitPoint.normal, rayDirection) * f, 0.0, 1.0);
			}
			else if (diffChance > 0 && roulette < specChance + diffChance)
			{
				// Diffuse reflection
				rayOrigin = hitPoint.position + hitPoint.normal * EPSILON;
				rayDirection = sampleHemisphere(hitPoint.normal, 1.0, hitPoint.position.zx+vec2(hitPoint.position.y)+vec2(seed, depth));
				energy *= hitPoint.material.albedo * clamp(dot(hitPoint.normal, rayDirection), 0.0, 1.0);
			} else {
				// This means both the hit material's albedo and specular are totally black, so there won't be anymore light. We can stop here.
				break;
			}
		} else {
			// The ray didn't hit anything :(
			totalIllumination += energy;
			break;
		}
	}

	return totalIllumination;
}

////////////////////////////////////////////// Main

void main() {

	/////////////////////////// 



	///////////////////////////

	vec3 rayDir = f_in.normal - u_cameraPosition;
	vec3 eye = normalize(-u_cameraPosition);
	Ray cameraRay = Ray(u_cameraPosition, f_in.position);

	// Camera raycasting
	vec3 colorSum = computeSceneColor(cameraRay, u_time); // calculate the method
	for (int i = 0; i<u_framePasses-1; i++) colorSum += computeSceneColor(cameraRay, u_time+i); 
	
	// get u_framepasses and u_time
	fb_color = vec4(colorSum / u_framePasses, 1.0);

	if (u_accumulatedPasses > 0) {

		SurfacePoint hitPoint;
		vec3 offsetDirection = cameraRay.direction + vec3(rand(vec2(1,u_time)+f_in.textureCoord), rand(vec2(2, u_time)+f_in.textureCoord), rand(vec2(3, u_time)+f_in.textureCoord));
			
		if (raycast(Ray(cameraRay.origin, offsetDirection), hitPoint)) {
			fb_color += vec4(hitPoint.material.emission*hitPoint.material.emissionStrength, 1.0);
		}

		// Add last frame back (progressive sampling)
		fb_color += texture(u_screenTexture, f_in.textureCoord);
	}
}