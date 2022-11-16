/*
 Ruled surface                
(C) Bedrich Benes 2022
Purdue University
bbenes@purdue.edu
*/

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"


#include <iostream>
#include "glad/glad.h"
#include "GLFW/glfw3.h"

#include "glm/glm.hpp"
#include "glm/gtc/type_ptr.hpp"
#include "glm/gtc/matrix_transform.hpp"

#include <stdio.h>
#include <iostream>
#include <string.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <time.h>
#include <string>
#include <vector>			//Standard template library
#include <array>

#include "triangle.h"  //triangles
#include "helper.h"         
#include "objGen.h"         //to save OBJ file format for 3D printing
#include "trackball.h"
#include <cstdlib>

#pragma warning(disable : 4996)
#pragma comment(lib, "glfw3.lib")


using namespace std;
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
//void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void processInput(GLFWwindow* window);
TrackBallC trackball;
static void KbdCallback(GLFWwindow* window, int key, int scancode, int action, int mods);
bool mouseLeft, mouseMid, mouseRight;

glm::vec3 cameraPos = glm::vec3(3.f, 0.f, 3.f);
glm::vec3 cameraFront = glm::vec3(-5.0f, -0.0f, -7.0f);
glm::vec3 cameraUp = glm::vec3(0.0f, 0.0f, 1.0f);
float deltaTime = 0.0f;	// time between current frame and last frame
float lastFrame = 0.0f;
int staticBall = 0;
float PI = 3.14159;
//glm::vec3 SPositions[3];
int bulletBall = 0;
int moveBall = 0;
float playerspeed = 2.5f;
std::vector<glm::vec3> SPositions;
std::vector<glm::vec3> BPositions;
std::vector<glm::vec3> BVelocity;
std::vector<glm::vec3> MPositions;
std::vector<glm::vec3> MVelocity;
vector <TriangleC> tri;   //all the triangles will be stored here
std::string filename = "geometry.obj";

GLuint points = 0; //number of points to display the object
GLuint points2 = 0;
int steps = 4;//# of subdivisions
int pointSize = 1;
int lineWidth = 1;
float r = .1f;
bool firstMouse = true;
float yaw = -90.0f;	// yaw is initialized to -90.0 degrees since a yaw of 0.0 results in a direction vector pointing to the right so we initially rotate a bit to the left.
float pitch = 0.0f;
float lastX = 800 / 2.0;
float lastY = 600 / 2.0;
float fov = 90.0f;

//Vertex array object and vertex buffer object indices 
GLuint VAO, VBO;
GLuint VAO2, VBO2;

std::vector<int> indices;
std::vector<int> lineIndices;
int k1, k2;

class water_surface {
private:
	GLuint points_vbo, normals_vbo;

public:
	static constexpr int N = 61;
	static constexpr int width = 60, height = 60;

	float c = 16.0;
	float u[width][height];
	float u_new[width][height];
	float v[width][height];
	float control_point_heights[width][height];

	GLuint vao, elements_vbo;

	glm::vec3 points[N * N];
	glm::vec3 normals[N * N];
	glm::vec3 elements[(N - 1) * (N - 1) * 2];

	float points_buffer[N * N * 3];
	float normals_buffer[N * N * 3];
	int elements_buffer[(N - 1) * (N - 1) * 2 * 3];

	water_surface() {
		glGenVertexArrays(1, &vao);
		glGenBuffers(1, &elements_vbo);
		glGenBuffers(1, &normals_vbo);
		glGenBuffers(1, &points_vbo);

		// Initialize water heights
		for (int i = 0; i < this->width; i++) {
			for (int j = 0; j < this->height; j++) {
				float x = -3 + 6 * (1 - (i / (float)this->width));
				float y = -3 + 6 * (1 - (j / (float)this->height));

				/*if (i == 30 && j == 30) {
					this->u[i][j] = 1.2;
				}
				else if ((i==29 && j ==29) || (i==29 && j == 31) || (i == 31 && j == 29)||(i==31 && j==31)) {
					this->u[i][j] = 0.4;
				}
				else if (i > 28 && i < 32 && j > 28 && j < 32) {
					this->u[i][j] = 0.7;
				}
				else if((i == 28 && j == 30) || (i == 32 && j == 30) || (i == 30 && j == 28) || (i == 30 && j == 32)) {
					this->u[i][j] = 0.2;
				}
				else {
					this->u[i][j] = 0;
				}*/
				this->u[i][j] = 0;

				this->v[i][j] = 0;
			}
		}
	}

	void update(float dt) {
		for (int i = 0; i < this->width; i++) {
			for (int j = 0; j < this->height; j++) {
				float v1, v2, v3, v4;

				if (i == 0) {
					v1 = this->u[i][j];
				}
				else {
					v1 = this->u[i - 1][j];
				}

				if (i == this->width - 1) {
					v2 = this->u[i][j];
				}
				else {
					v2 = this->u[i + 1][j];
				}

				if (j == 0) {
					v3 = this->u[i][j];
				}
				else {
					v3 = this->u[i][j - 1];
				}

				if (j == this->height - 1) {
					v4 = this->u[i][j];
				}
				else {
					v4 = this->u[i][j + 1];
				}

				float f = c * c * ((v1 + v2 + v3 + v4) - 4 * this->u[i][j]);
				this->v[i][j] += f * dt;
				this->v[i][j] *= 0.997;
				this->u_new[i][j] = u[i][j] + v[i][j] * dt;
			}
		}

		for (int i = 0; i < this->width; i++) {
			for (int j = 0; j < this->height; j++) {
				this->u[i][j] = this->u_new[i][j];
				this->control_point_heights[i][j] = this->u[i][j];
			}
		}

		static int x_delta[9] = { 0, -1, -1, 0, 1, 1, 1, 0, -1 };
		static int y_delta[9] = { 0, 0, -1, -1, -1, 0, 1, 1, 1 };

		for (int i = 3; i < this->width; i += 3) {
			for (int j = 3; j < this->height; j += 3) {
				glm::vec3 points[9];

				for (int k = 0; k < 9; k++) {
					int x_index = i + x_delta[k];
					int y_index = j + y_delta[k];

					points[k].x = -3 + 6 * (1 - (x_index / (float)this->width));
					points[k].y = -3 + 6 * (1 - (y_index / (float)this->height));
					points[k].z = this->control_point_heights[x_index][y_index];
				}

				float sum_xx = 0.0;
				float sum_yy = 0.0;
				float sum_xy = 0.0;
				float sum_yz = 0.0;
				float sum_xz = 0.0;

				for (int k = 0; k < 9; k++) {
					sum_xx += points[k].x * points[k].x;
					sum_yy += points[k].y * points[k].y;
					sum_xy += points[k].x * points[k].y;
					sum_yz += points[k].y * points[k].z;
					sum_xz += points[k].x * points[k].z;
				}

				float D = sum_xx * sum_yy - sum_xy * sum_xy;
				float a = (sum_yz * sum_xy - sum_xz * sum_yy) / D;
				float b = (sum_xy * sum_xz - sum_xx * sum_yz) / D;

				glm::vec3 n(a, b, 1);
				glm::vec3 p = points[0];

				for (int k = 1; k < 9; k++) {
					glm::vec3 p0 = points[k];

					float z = (n.x * (p.x - p0.x) + n.y * (p.y - p0.y)) / n.z + p.z;

					int x_index = i + x_delta[k];
					int y_index = j + y_delta[k];
					this->control_point_heights[x_index][y_index] = z;
				}
			}
		}

		// Calculate points
		int point_index = 0;
		for (float x = -3; x < 3; x += 0.1) {
			for (float y = -3; y < 3; y += 0.1) {
				int i = ((x + 3.0) / 6.0) * this->width;
				int j = ((y + 3.0) / 6.0) * this->height;

				// Take the average height of the blocks around us, weighted by their distances
				float sum_of_weights = 0.0;
				float avg_height = 0.0;
				for (int k = i - 1; k <= i + 1; k++) {
					if (k < 0 || k >= this->width) {
						continue;
					}

					for (int l = j - 1; l <= j + 1; l++) {
						if (l < 0 || l >= this->height) {
							continue;
						}

						float u = -3 + 6 * (k / (float)this->width);
						float v = -3 + 6 * (l / (float)this->height);

						// Make sure the weight is larger for smaller distances
						float weight = 100 - (u - x) * (u - x) + (y - v) * (y - v);
						avg_height += this->u[k][l] * weight;
						sum_of_weights += weight;
					}
				}
				avg_height /= sum_of_weights;

				points[point_index++] = glm::vec3(x, y, avg_height);
			}
		}

		// Calculate normals
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				// Average the normals for each triangle around us
				int p1_i = i + N * j;
				int p2_i = i + 1 + N * j;
				int p3_i = i - 1 + N * j;
				int p4_i = i + N * (j - 1);
				int p5_i = i + N * (j + 1);

				glm::vec3 normal;

				if (i == 0 && j == 0) {
					glm::vec3 p1 = points[p1_i];
					glm::vec3 p2 = points[p2_i];
					glm::vec3 p5 = points[p5_i];

					normal = glm::cross(p5 - p1, p2 - p1);
				}
				else if (i == 0 && j == N - 1) {
					glm::vec3 p1 = points[p1_i];
					glm::vec3 p2 = points[p2_i];
					glm::vec3 p4 = points[p4_i];

					normal = glm::cross(p2 - p1, p4 - p1);
				}
				else if (i == N - 1 && j == 0) {
					glm::vec3 p1 = points[p1_i];
					glm::vec3 p3 = points[p3_i];
					glm::vec3 p5 = points[p5_i];

					normal = glm::cross(p3 - p1, p5 - p1);
				}
				else if (i == N - 1 && j == N - 1) {
					glm::vec3 p1 = points[p1_i];
					glm::vec3 p3 = points[p3_i];
					glm::vec3 p4 = points[p4_i];

					normal = glm::cross(p4 - p1, p3 - p1);
				}
				else if (i == 0) {
					glm::vec3 p1 = points[p1_i];
					glm::vec3 p2 = points[p2_i];
					glm::vec3 p4 = points[p4_i];
					glm::vec3 p5 = points[p5_i];

					glm::vec3 n1 = glm::cross(p2 - p1, p4 - p1);
					n1 = glm::normalize(n1);
					//n1.make_unit_length();
					glm::vec3 n2 = glm::cross(p5 - p1, p2 - p1);
					n2 = glm::normalize(n2);
					//n2.make_unit_length();

					normal = 0.5f * (n1 + n2);
				}
				else if (j == 0) {
					glm::vec3 p1 = points[p1_i];
					glm::vec3 p2 = points[p2_i];
					glm::vec3 p3 = points[p3_i];
					glm::vec3 p5 = points[p5_i];

					glm::vec3 n1 = glm::cross(p3 - p1, p5 - p1);
					n1 = glm::normalize(n1);
					glm::vec3 n2 = glm::cross(p5 - p1, p2 - p1);
					n2 = glm::normalize(n2);

					normal = 0.5f * (n1 + n2);
				}
				else if (i == N - 1) {
					glm::vec3 p1 = points[p1_i];
					glm::vec3 p3 = points[p3_i];
					glm::vec3 p4 = points[p4_i];
					glm::vec3 p5 = points[p5_i];

					glm::vec3 n1 = glm::cross(p4 - p1, p3 - p1);
					n1 = glm::normalize(n1);
					glm::vec3 n2 = glm::cross(p3 - p1, p5 - p1);
					n2 = glm::normalize(n2);

					normal = 0.5f * (n1 + n2);
				}
				else if (j == N - 1) {
					glm::vec3 p1 = points[p1_i];
					glm::vec3 p2 = points[p2_i];
					glm::vec3 p3 = points[p3_i];
					glm::vec3 p4 = points[p4_i];

					glm::vec3 n1 = glm::cross(p2 - p1, p4 - p1);
					n1 = glm::normalize(n1);
					glm::vec3 n2 = glm::cross(p4 - p1, p3 - p1);
					n2 = glm::normalize(n2);

					normal = 0.5f * (n1 + n2);
				}
				else {
					glm::vec3 p1 = points[p1_i];
					glm::vec3 p2 = points[p2_i];
					glm::vec3 p3 = points[p3_i];
					glm::vec3 p4 = points[p4_i];
					glm::vec3 p5 = points[p5_i];

					glm::vec3 n1 = glm::cross(p4 - p1, p3 - p1);
					n1 = glm::normalize(n1);
					glm::vec3 n2 = glm::cross(p2 - p1, p4 - p1);
					n2 = glm::normalize(n2);
					glm::vec3 n3 = glm::cross(p5 - p1, p2 - p1);
					n3 = glm::normalize(n3);
					glm::vec3 n4 = glm::cross(p3 - p1, p5 - p1);
					n4 = glm::normalize(n4);

					normal = 0.25f * (n1 + n2 + n3 + n4);
				}

				normal = glm::normalize(normal);
				normals[p1_i] = normal;
			}
		}

		// Calculate the elements for each triangle
		int e_i = 0;
		for (int i = 0; i < N - 1; i++) {
			for (int j = 0; j < N - 1; j++) {
				// First triangle
				int p1 = i + N * j;
				int p2 = i + 1 + N * j;
				int p3 = i + N * (j + 1);

				elements[e_i++] = glm::vec3(p1, p2, p3);

				// Second triangle
				int p4 = i + 1 + N * j;
				int p5 = i + 1 + N * (j + 1);
				int p6 = i + N * (j + 1);

				elements[e_i++] = glm::vec3(p4, p5, p6);
			}
		}

		glBindVertexArray(vao);

		// Set up the points vbo
		for (int i = 0; i < N * N; i++) {
			points_buffer[i * 3] = points[i].x;
			points_buffer[i * 3 + 1] = points[i].y;
			points_buffer[i * 3 + 2] = points[i].z;
		}
		glBindBuffer(GL_ARRAY_BUFFER, points_vbo);
		glBufferData(GL_ARRAY_BUFFER, sizeof(float) * N * N * 3, points_buffer, GL_STATIC_DRAW);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);
		glEnableVertexAttribArray(0);

		// Set up the normals vbo
		for (int i = 0; i < N * N; i++) {
			normals_buffer[i * 3] = normals[i].x;
			normals_buffer[i * 3 + 1] = normals[i].y;
			normals_buffer[i * 3 + 2] = normals[i].z;
		}
		glBindBuffer(GL_ARRAY_BUFFER, normals_vbo);
		glBufferData(GL_ARRAY_BUFFER, sizeof(float) * N * N * 3, normals_buffer, GL_STATIC_DRAW);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, NULL);
		glEnableVertexAttribArray(1);


		// Set up elements vbo
		for (int i = 0; i < (N - 1) * (N - 1) * 2; i++) {
			elements_buffer[i * 3] = elements[i].x;
			elements_buffer[i * 3 + 1] = elements[i].y;
			elements_buffer[i * 3 + 2] = elements[i].z;
		}
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, elements_vbo);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(int) * (N - 1) * (N - 1) * 2 * 3, elements_buffer, GL_STATIC_DRAW);

		glBindVertexArray(0);
	}
};

inline void AddVertex(vector <GLfloat>* a, glm::vec3 A) //check this!
{
	a->push_back(A[0]); a->push_back(A[1]); a->push_back(A[2]);
}


glm::vec3 P(GLfloat t)
{
	return glm::vec3(0.3 * cos(2 * M_PI * t + M_PI / 2), 0, 0.6 * sin(2 * M_PI * t + M_PI / 2));
}

inline glm::vec3 Q(GLfloat t)
{
	return glm::vec3(0.6 * cos(2 * M_PI * t), 1, 0.6 * sin(2 * M_PI * t));
}

inline glm::vec3 S(GLfloat u, GLfloat t)
{
	return glm::vec3(u * P(t) + (1 - u) * Q(t));
}
void CreatePlane(vector <GLfloat>* vv, int n, float t, float sharp) {
	GLfloat slices = 60/2;
	GLfloat stacks = 60/2;
	GLfloat deltaTheta = 1.0f/n;
	GLfloat deltaPhi = 1.0f/n;
	for (int i = 0; i < stacks*n; i++)
	{
		GLfloat phi = i * deltaPhi -30/2;
		for (int j = 0; j < slices*n; j++)
		{
			GLfloat theta = j * deltaTheta -30/2 ;
			float k = 1.0f;
			float f = k * (theta + phi - t);
			//float f2 = k * (phi - t);
			float a = sharp / k;
			float A = 1.0f;

			glm::vec3 ori = glm::vec3(0.0,0.0,0.0);
			glm::vec3 V_LD = glm::vec3(theta + a * cos(f), a * sin(f), phi + a * cos(f));
			glm::vec3 V_RD = glm::vec3(theta + deltaTheta + a * cos(f + deltaTheta), a * sin(f + deltaTheta), phi + a * cos(f + deltaTheta));
			glm::vec3 V_LU = glm::vec3(theta + a * cos(f + deltaPhi), a * sin(f + deltaPhi), phi + deltaPhi + a * cos(f + deltaPhi));
			glm::vec3 V_RU = glm::vec3(theta + deltaTheta + a * cos(f + deltaTheta + deltaPhi), a * sin(f + deltaTheta + deltaPhi), phi + deltaPhi + a * cos(f + deltaTheta + deltaPhi));

			float D_LD = glm::distance(V_LD, ori);
			float D_RD = glm::distance(V_RD, ori);
			float D_LU = glm::distance(V_LU, ori);
			float D_RU = glm::distance(V_RU, ori);
			float D_LD_ = 1;
			float D_RD_ = 1;
			float D_LU_ = 1; 
			float D_RU_ = 1;
			if (D_LD > t) {
				D_LD_ = 0;
			}
			if (D_RD > t) {
				D_RD_ = 0;
			}
			if (D_LU > t) {
				D_LU_ = 0;
			}
			if (D_RU > t) {
				D_RU_ = 0;
			}
			if (t > 20) {
				A = 0;
			}
			//lower triangle
			AddVertex(vv, V_LD + glm::vec3(0.0, (1-t*0.05)*A * (1 / ((D_LD+1))) * sin(-PI * D_LD_*(t-D_LD)), 0.0));
			AddVertex(vv, V_RD + glm::vec3(0.0, (1 - t * 0.05) * A * (1 / ((D_RD + 1))) * sin(-PI * D_RD_ * (t-D_RD)), 0.0));
			AddVertex(vv, V_LU + glm::vec3(0.0, (1 - t * 0.05) * A * (1 / ((D_LU + 1))) * sin(-PI * D_LU_ * (t-D_LU) ),0.0));
			//upper triangle
			AddVertex(vv, V_RD + glm::vec3(0.0, (1 - t * 0.05) * A * (1 / ((D_RD + 1))) * sin(-PI * D_RD_ * (t - D_RD)),0.0));
			AddVertex(vv, V_LU + glm::vec3(0.0, (1 - t * 0.05) * A * (1 / ((D_LU + 1))) * sin(-PI * D_LU_ * (t - D_LU)),0.0));
			AddVertex(vv, V_RU + glm::vec3(0.0, (1 - t * 0.05) * A * (1 / ((D_RU + 1))) * sin(-PI * D_RU_ * (t - D_RU)),0.0));

		}
	}
}

void CreateSphere(vector <GLfloat>* vv, int n) {
	GLfloat step = 1.f / n;
	GLfloat slices = 30;
	GLfloat stacks = 20;
	GLfloat deltaTheta = 2 * M_PI / (GLfloat)slices;
	GLfloat deltaPhi = M_PI / (GLfloat)stacks;
	for (int i = 0; i < stacks; i++)
	{
		GLfloat phi = i * deltaPhi;
		for (int j = 0; j < slices; j++)
		{
			GLfloat theta = j * deltaTheta;

			//lower triangle
			AddVertex(vv, glm::vec3(r * cos(theta) * sin(phi),
				r * sin(theta) * sin(phi),
				r * cos(phi)));
			AddVertex(vv, glm::vec3(r * cos(theta + deltaTheta) * sin(phi),
				r * sin(theta + deltaTheta) * sin(phi),
				r * cos(phi)));
			AddVertex(vv, glm::vec3(r * cos(theta) * sin(phi + deltaPhi),
				r * sin(theta) * sin(phi + deltaPhi),
				r * cos(phi + deltaPhi)));
			//upper triangle
			AddVertex(vv, glm::vec3(r * cos(theta + deltaTheta) * sin(phi),
				r * sin(theta + deltaTheta) * sin(phi),
				r * cos(phi)));
			AddVertex(vv, glm::vec3(r * cos(theta) * sin(phi + deltaPhi),
				r * sin(theta) * sin(phi + deltaPhi),
				r * cos(phi + deltaPhi)));
			AddVertex(vv, glm::vec3(r * cos(theta + deltaTheta) * sin(phi + deltaPhi),
				r * sin(theta + deltaTheta) * sin(phi + deltaPhi),
				r * cos(phi + deltaPhi)));

		}
	}
}
void CreateRuled(vector <GLfloat>* vv, int n)
{
	GLfloat step = 1.f / n;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			//lower triangle
			AddVertex(vv, S(i * step, j * step));
			AddVertex(vv, S((i + 1) * step, j * step));
			AddVertex(vv, S((i + 1) * step, (j + 1) * step));
			//upper triangle
			AddVertex(vv, S(i * step, j * step));
			AddVertex(vv, S((i + 1) * step, (j + 1) * step));
			AddVertex(vv, S(i * step, (j + 1) * step));
		}
	}
}

int CompileShaders() {
	//Vertex Shader
	const char* vsSrc= "#version 330 core\n"
		"layout (location = 0) in vec4 iPos;\n"
		"layout (location = 1) in vec3 aNormal;\n"
		"uniform mat4 modelview;\n"
		"void main()\n"
		"{\n"
		"   vec4 oPos=modelview*iPos;\n"
		"   gl_Position = vec4(oPos.x, oPos.y, oPos.z, oPos.w);\n"
		"}\0";

	//Fragment Shader
	const char* fsSrc = "#version 330 core\n"
		"out vec4 col;\n"
		"uniform vec4 color;\n"
		"void main()\n"
		"{\n"
		"   col = color;\n"
		"}\n\0";

	//Create VS object
	GLuint vs = glCreateShader(GL_VERTEX_SHADER);
	//Attach VS src to the Vertex Shader Object
	glShaderSource(vs, 1, &vsSrc, NULL);
	//Compile the vs
	glCompileShader(vs);

	//The same for FS
	GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(fs, 1, &fsSrc, NULL);
	glCompileShader(fs);

	//Get shader program object
	GLuint shaderProg = glCreateProgram();
	//Attach both vs and fs
	glAttachShader(shaderProg, vs);
	glAttachShader(shaderProg, fs);
	//Link all
	glLinkProgram(shaderProg);

	//Clear the VS and FS objects
	glDeleteShader(vs);
	glDeleteShader(fs);
	return shaderProg;
}

void BuildScene2(GLuint& VBO, GLuint& VAO, int n, float t, float sharp) { //return VBO and VAO values n is the subdivision
	vector<GLfloat> v;
	CreatePlane(&v, n, t, sharp);
	//now get it ready for saving as OBJ
	

	//make VAO
	glGenVertexArrays(1, &VAO2);
	glGenBuffers(1, &VBO2);

	//bind it
	glBindVertexArray(VAO2);

	//bind the VBO
	glBindBuffer(GL_ARRAY_BUFFER, VBO2);
	//send the data to the GPU
	points2 = v.size();
	glBufferData(GL_ARRAY_BUFFER, points2 * sizeof(GLfloat), &v[0], GL_STATIC_DRAW);

	//Configure the attributes
//	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
	glVertexAttribPointer((GLuint)0, 3, GL_FLOAT, GL_FALSE, 0, 0);
	//Make it valid
	glEnableVertexAttribArray(0);

	v.clear(); //no need for the data, it is on the GPU now

}
void BuildScene(GLuint& VBO, GLuint& VAO, int n) { //return VBO and VAO values n is the subdivision
	vector<GLfloat> v;
	CreateSphere(&v, n);
	//now get it ready for saving as OBJ
	tri.clear();
	for (unsigned int i = 0; i < v.size(); i += 9) { //stride 3 - 3 vertices per triangle
		TriangleC tmp;
		glm::vec3 a, b, c;
		a=glm::vec3(v[i], v[i + 1], v[i + 2]);
		b=glm::vec3(v[i + 3], v[i + 4], v[i + 5]);
		c = glm::vec3(v[i + 6], v[i + 7], v[i + 8]);
		tmp.Set(a, b, c); //store them for 3D export
		tri.push_back(tmp);
	}

	//make VAO
	glGenVertexArrays(1, &VAO);
	glGenBuffers(1, &VBO);

	//bind it
	glBindVertexArray(VAO);

	//bind the VBO
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	//send the data to the GPU
	points = v.size();
	glBufferData(GL_ARRAY_BUFFER, points * sizeof(GLfloat), &v[0], GL_STATIC_DRAW);

	//Configure the attributes
//	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
	glVertexAttribPointer((GLuint)0, 3, GL_FLOAT, GL_FALSE, 0, 0);
	//Make it valid
	glEnableVertexAttribArray(0);

	v.clear(); //no need for the data, it is on the GPU now

}

//Quit when ESC is released
static void KbdCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if (key == GLFW_KEY_SPACE && action == GLFW_PRESS) {
		if (bulletBall < 3) {
			BPositions.push_back(cameraPos);
			BVelocity.push_back(20.0f * cameraFront);
			bulletBall++;
		}
	}
}


int main()
{
	srand(time(NULL));
	glfwInit();

	//negotiate with the OpenGL
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	//make OpenGL window
	GLFWwindow* window = glfwCreateWindow(1000,1000, "Simple", NULL, NULL);
	//is all OK?
	lastX = 500;
	lastY = 500;
	if (window == NULL)
	{
		std::cout << "Cannot open GLFW window" << std::endl;
		glfwTerminate();
		return -1;
	}
	//Paste the window to the current context
	glfwMakeContextCurrent(window);

	//Load GLAD to configure OpenGL
	gladLoadGL();
	//Set the viewport
	//glViewport(0, 0, 800, 800);
	float X = 800;
	float Y = 800;
	//once the OpenGL context is done, build the scene and compile shaders
	//BuildScene2(VBO2, VAO2, steps, 0, 0.3);
	water_surface water_surface;
	//BuildScene(VBO, VAO, steps);
	int shaderProg = CompileShaders();
	GLint modelviewParameter = glGetUniformLocation(shaderProg, "modelview");
	//Background color
	glClearColor(0.2f, 0.2f, 0.2f, 1.0f);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	//Use shader
	glUseProgram(shaderProg);
	glPointSize(pointSize);

	// Initialize ImGUI
	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImGuiIO& io = ImGui::GetIO(); (void)io;
	ImGui::StyleColorsDark();
	ImGui_ImplGlfw_InitForOpenGL(window, true);
	ImGui_ImplOpenGL3_Init("#version 330");

	bool drawScene = true;
	float color[4] = { 0.8f, 0.8f, 0.2f, 1.0f };
	//send the color to the fragment shader
	glUniform4f(glGetUniformLocation(shaderProg, "color"), color[0], color[1], color[2], color[3]);


	glfwSetKeyCallback(window, KbdCallback); //set keyboard callback to quit
	//glfwSetCursorPosCallback(window, mouse_callback);
	//glfwSetMouseButtonCallback(window, MouseButtonCallback);;
	//glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
	//glEnable(GL_DEPTH_TEST);

	GLfloat ver[] = { 500,500 };
	// Main while loop
	bool is_drawing_continous = true;
	
	while (!glfwWindowShouldClose(window))
	{
		
		float currentFrame = static_cast<float>(glfwGetTime());
		deltaTime = currentFrame - lastFrame;
		lastFrame = currentFrame;
		

		//Clean the window
		processInput(window);
		//processInput(window);
		//glClear(GL_COLOR_BUFFER_BIT);
		glClear(GL_COLOR_BUFFER_BIT );
		
		
		
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();
		
		//BuildScene2(VBO, VAO, steps, currentFrame, 0.3);
		
		
		
		//if (ImGui::SliderInt("point Size", &pointSize, 1, 10, "%d", 0)) {
		glPointSize(pointSize); //set the new point size if it has been changed			
		//}
		//if (ImGui::SliderInt("line width", &lineWidth, 1, 10, "%d", 0)) {
		glLineWidth(lineWidth); //set the new point size if it has been changed			
		//}
		//if (ImGui::ColorEdit4("Color", color)) { //set the new color only if it has changed
		glUniform4f(glGetUniformLocation(shaderProg, "color"), color[0], color[1], color[2], color[3]);
		//}

		//set the projection matrix
		glm::mat4 proj = glm::perspective(fov, 1.0f, 0.1f, 100.0f);
		//set the viewing matrix (looking from [0,0,5] to [0,0,0])
		glm::mat4 view = glm::lookAt(cameraPos, cameraPos + cameraFront, cameraUp);
		
		{
			static bool is_mouse_down = false;

			double xpos, ypos;
			glfwGetCursorPos(window, &xpos, &ypos);
			float x = (2.0f * xpos) / 1000 - 1.0f;
			float y = 1.0f - (2.0f * ypos) / 1000;
			float z = 1.0f;
			glm::vec3 ray_nds = glm::vec3(x, y, z);
			glm::vec4 ray_clip = glm::vec4(ray_nds.x, ray_nds.y, -1.0, 1.0);
			glm::vec4 ray_eye = glm::inverse(proj) * ray_clip;
			ray_eye = glm::vec4(ray_eye.x, ray_eye.y, -1.0, 0.0);
			glm::vec4 ray_wor = (glm::inverse(view) * ray_eye);
			//glm::vec4 mouse_world1 = (glm::inverse(view) * glm::inverse(proj) * glm::vec4((2.0 * xpos) / 1000 - 1.0, 1.0 - (2.0 * ypos) /1000, -1.0, 1.0));
			//glm::vec3 mouse_world = (mouse_world1.x, mouse_world1.y, mouse_world1.z);
			glm::vec3 mouse_world = glm::normalize(glm::vec3(ray_wor.x, ray_wor.y, ray_wor.z));
			glm::vec3 mouse_intersection = cameraPos + (-cameraPos.z / mouse_world.z) * mouse_world;

			int mouse_state = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
			if (mouse_state == GLFW_PRESS && !is_mouse_down) {
				is_mouse_down = true;

				if (mouse_intersection.x > -3.0 &&
					mouse_intersection.x < 3.0 &&
					mouse_intersection.y > -3.0 &&
					mouse_intersection.y < 3.0) {
					int i = (mouse_intersection.x + 3.0) / 6.0 * water_surface.width;
					int j = (mouse_intersection.y + 3.0) / 6.0 * water_surface.height;

					if (i > 0 && j > 0 && i < water_surface.width - 1 && j < water_surface.height - 1) {
						water_surface.u[i][j] = 1.2;
						water_surface.u[i - 1][j - 1] = 0.7;
						water_surface.u[i - 1][j] = 0.7;
						water_surface.u[i - 1][j + 1] = 0.7;
						water_surface.u[i + 1][j - 1] = 0.7;
						water_surface.u[i + 1][j] = 0.7;
						water_surface.u[i + 1][j + 1] = 0.7;
						water_surface.u[i][j + 1] = 0.7;
						water_surface.u[i][j - 1] = 0.7;
					}
				}
			}
			else if (mouse_state == GLFW_RELEASE && is_mouse_down) {
				is_mouse_down = false;
			}
		}

		
		water_surface.update(deltaTime);

		glUniform4f(glGetUniformLocation(shaderProg, "color"), 1.0f, 1.0f, 1.0f, 0.0f);
		glm::mat4 model2 = glm::mat4(1.0f);
		model2 = glm::translate(model2, glm::vec3(0.0f, 0.0f, 0.0f));
		glm::mat4 modelView2 = proj * view * model2;
		glUniformMatrix4fv(modelviewParameter, 1, GL_FALSE, glm::value_ptr(modelView2));
		glBindVertexArray(water_surface.vao);
		//glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, water_surface.elements_vbo);
		glDrawElements(GL_TRIANGLES, (water_surface.N - 1) * (water_surface.N - 1) * 2 * 3, GL_UNSIGNED_INT, 0);

		// Renders the ImGUI elements
		ImGui::Render();
		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
		//Swap the back buffer with the front buffer
		glfwSwapBuffers(window);
		//make sure events are served
		glfwPollEvents();
	}
	//Cleanup
	glDeleteVertexArrays(1, &VAO2);
	glDeleteBuffers(1, &VBO2);
	glDeleteProgram(shaderProg);
	glfwDestroyWindow(window);
	glfwTerminate();
	return 0;
}

void mouse_callback(GLFWwindow* window, double xposIn, double yposIn)
{
	float xpos = static_cast<float>(xposIn);
	float ypos = static_cast<float>(yposIn);

	if (firstMouse)
	{
		lastX = xpos;
		lastY = ypos;
		firstMouse = false;
	}

	float xoffset = xpos - lastX;
	float yoffset = lastY - ypos; // reversed since y-coordinates go from bottom to top
	lastX = xpos;
	lastY = ypos;

	float sensitivity = 0.01f; // change this value to your liking
	xoffset *= sensitivity;
	yoffset *= sensitivity;

	yaw += xoffset;
	pitch += yoffset;

	// make sure that when pitch is out of bounds, screen doesn't get flipped
	if (pitch > 89.0f)
		pitch = 89.0f;
	if (pitch < -89.0f)
		pitch = -89.0f;

	glm::vec3 front;
	front.x = cos(glm::radians(yaw)) * cos(glm::radians(pitch));
	front.z = sin(glm::radians(pitch));
	front.y = sin(glm::radians(yaw)) * cos(glm::radians(pitch));
	cameraFront = glm::normalize(front);
}

void processInput(GLFWwindow* window)
{
	if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
		glfwSetWindowShouldClose(window, true);

	float cameraSpeed = playerspeed * deltaTime;
	if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
		//cout << "  w" << endl;
		cameraPos += cameraSpeed * cameraFront /4.0f;
	if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
		cameraPos -= cameraSpeed * cameraFront/4.0f;
	if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
		cameraPos -= glm::normalize(glm::cross(cameraFront, cameraUp)) * cameraSpeed;
	if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
		cameraPos += glm::normalize(glm::cross(cameraFront, cameraUp)) * cameraSpeed;
	if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS)
		cameraPos += glm::normalize(glm::normalize(cameraUp)) * cameraSpeed;
	if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS)
		cameraPos -= glm::normalize(glm::normalize(cameraUp)) * cameraSpeed;
	//cameraPos.y = 1.0f;
}