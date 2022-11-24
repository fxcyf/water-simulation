#version 330 core
layout (location = 0) in vec3 waterPos;
layout (location = 1) in vec3 causticPos;
out vec3 oldPos;
out vec3 newPos;
uniform mat4 modelview;
uniform mat4 model;
void main()
{
	oldPos = waterPos;
	newPos = causticPos;
   vec4 oPos = modelview * vec4(causticPos, 1.f);
   gl_Position = oPos;
}