#version 330 core
#extension GL_OES_standard_derivatives : enable
out vec4 col;
in vec3 oldPos;
in vec3 newPos;
//in bool gl_FrontFacing;
uniform vec3 color;
uniform float causticsStrength;
uniform float causticsAlpha;
void main()
{
//	float causticsStrength = 0.8;
//	float causticsAlpha = 0.6;
	float causticsIntensity = 0.0;
	float oldArea = length(dFdx(oldPos)) * length(dFdy(oldPos));
	float newArea = length(dFdx(newPos)) * length(dFdy(newPos));
	float ratio = oldArea / newArea;
	if (newArea == 0.) {
		ratio = 2.0e+20;
	} else {
		ratio = oldArea / newArea;
	}
	causticsIntensity = causticsStrength * ratio;
	if (gl_FrontFacing) {
		col = vec4(color * causticsIntensity, causticsAlpha);
	} else {
		col = vec4(color * causticsIntensity * 5, causticsAlpha / 5);
	}
}