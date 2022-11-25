#version 330 core
layout (location = 0) in vec4 iPos;
layout (location = 1) in vec3 aNormal;
layout (location = 2) in vec2 aTexCoord;

out vec3 FragPos;
out vec3 Normal;
out vec2 TexCoord;

uniform mat4 modelview;
uniform mat4 model;
uniform mat3 normalMat;
void main()
{
   vec4 oPos = modelview * iPos;
   gl_Position = vec4(oPos.x, oPos.y, oPos.z, oPos.w);
   FragPos = vec3(model * iPos);
   Normal = normalMat * normalize(aNormal);
   TexCoord = vec2(aTexCoord.x, aTexCoord.y);

}