#version 330 core
#extension GL_OES_standard_derivatives : enable
out vec4 col;
in vec3 oldPos;
in vec3 newPos;
uniform vec3 color;
void main()
{
	float causticsFactor = 0.8;
	float causticsIntensity = 0.0;
	float oldArea = length(dFdx(oldPos)) * length(dFdy(oldPos));
	float newArea = length(dFdx(newPos)) * length(dFdy(newPos));
	float ratio = oldArea / newArea;
	if (newArea == 0.) {
		ratio = 2.0e+20;
	} else {
		ratio = oldArea / newArea;
	}
	causticsIntensity = causticsFactor * ratio;
	col = vec4(color * causticsIntensity, 0.6f);
}