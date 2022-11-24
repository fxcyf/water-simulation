#version 330 core
out vec4 col;
in vec3 Normal;
in vec3 FragPos;
uniform vec3 color;
uniform vec3 lightColor;
uniform vec3 lightPos;
uniform vec3 viewPos;
uniform sampler2D floorTexture;
void main()
{
	float ambientStrength = 0.1;
	vec3 ambient = ambientStrength * lightColor;

	float diffuseStrength = 0.8;
	vec3 norm = normalize(Normal);
	vec3 lightDir = normalize(lightPos - FragPos);
	float diff = max(dot(norm, lightDir), 0.0);
	vec3 diffuse = diffuseStrength * diff * lightColor;

	float specularStrength = 0.5;
	vec3 viewDir = normalize(viewPos - FragPos);
	vec3 reflectDir = reflect(-lightDir, norm);
	float spec = pow(max(dot(viewDir, reflectDir), 0.0), 32);
	vec3 specular = specularStrength * spec * lightColor;

//	float floor_kd = diffuseStrength * diff;

//	vec3 incident = normalize(viewPos - FragPos);
//	vec3 refraction = refract(incident, Normal, 0.7);

//	vec3 floor_hit_point = FragPos + refraction * ((3 - FragPos.z) / refraction.z);
//	vec2 floor_texture_coord = ((floor_hit_point + vec3(3.0, 3.0, 0.0)) / 6.0).xy;

	vec3 result = (ambient + diffuse + specular) * color;
   col = vec4(result, 0.8f);
}