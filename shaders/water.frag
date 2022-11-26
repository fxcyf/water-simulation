#version 330 core
out vec4 col;
in vec3 Normal;
in vec3 FragPos;
in vec2 FloorTexCoord;
in vec2 SkyTexCoord;

uniform vec3 color;
uniform vec3 lightColor;
uniform vec3 lightPos;
uniform vec3 viewPos;
uniform sampler2D floorTexture;
uniform sampler2D skyTexture;
uniform float ambientStrength;
uniform float diffuseStrength;
uniform float floorTextureStrength;
uniform float skyTextureStrength;
uniform float colorStrength;
uniform float waterAlpha;
void main()
{
//	float ambientStrength = 0.2;
	vec3 ambient = ambientStrength * lightColor;

//	float diffuseStrength = 0.6;
	vec3 norm = normalize(Normal);
	vec3 lightDir = normalize(lightPos - FragPos);
	float diff = max(dot(norm, lightDir), 0.0);
	vec3 diffuse = diffuseStrength * diff * lightColor;

	float specularStrength = 0.5;
	vec3 viewDir = normalize(viewPos - FragPos);
	vec3 reflectDir = reflect(-lightDir, norm);
	float spec = pow(max(dot(viewDir, reflectDir), 0.0), 32);
	vec3 specular = specularStrength * spec * lightColor;

	float kd = diffuseStrength * diff;

//	float textureStrength = 0.8;
//	float colorStrength = 0.6;
//	float waterAlpha = 0.5;
    vec3 floor_texel = texture(floorTexture, FloorTexCoord).xyz;
	vec3 sky_texel = texture(skyTexture, SkyTexCoord).xyz;

	vec3 result = colorStrength * ((ambient + diffuse + specular) * color + kd * (floorTextureStrength * floor_texel + skyTextureStrength * sky_texel));
//	vec3 result = floor_kd * 0.5 * floor_texel;

   col = vec4(result, waterAlpha);
}