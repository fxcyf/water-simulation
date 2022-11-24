#version 330 core
layout (location = 0) in vec4 iPos;
uniform mat4 modelview;
void main()
{
   vec4 oPos = modelview*iPos;
   gl_Position = vec4(oPos.x, oPos.y, oPos.z, oPos.w);
}