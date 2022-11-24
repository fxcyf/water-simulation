#version 330 core
out vec4 col;
in vec3 Normal;
in vec3 FragPos;
in vec2 TexCoord;

uniform vec3 color;
uniform vec3 lightColor;
uniform vec3 lightPos;
uniform sampler2D floorTexture;
void main()
{
	float ambientStrength = 0.3;
	vec3 ambient = ambientStrength * lightColor;

	float diffuseStrength = 0.5;
	vec3 norm = normalize(Normal);
	vec3 lightDir = normalize(lightPos - FragPos);
	float diff = max(dot(norm, lightDir), 0.0);
	vec3 diffuse = diffuseStrength * diff * lightColor;

	vec3 texel = texture(floorTexture, TexCoord).xyz;
	vec3 result = ( ambient + diffuse ) * texel;
   col = vec4(result, 0.8f);
}