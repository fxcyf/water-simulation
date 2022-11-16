/*
Water Simulation
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
void processInput(GLFWwindow* window);
static void KbdCallback(GLFWwindow* window, int key, int scancode, int action, int mods);

glm::vec3 cameraPos = glm::vec3(3.f, 0.f, 3.f);
glm::vec3 cameraFront = glm::vec3(-5.0f, -0.0f, -7.0f);
glm::vec3 cameraUp = glm::vec3(0.0f, 0.0f, 1.0f);
float deltaTime = 0.0f;	// time between current frame and last frame
float lastFrame = 0.0f;
float PI = 3.14159;

float playerspeed = 2.5f;
int pointSize = 1;
int lineWidth = 1;
bool firstMouse = true;
float yaw = -90.0f;	// yaw is initialized to -90.0 degrees since a yaw of 0.0 results in a direction vector pointing to the right so we initially rotate a bit to the left.
float pitch = 0.0f;
float lastX = 800 / 2.0;
float lastY = 600 / 2.0;
float fov = 90.0f;

glm::vec3 waterColor = glm::vec3(0.30f, 0.51f, 0.66);
glm::vec3 lightColor = glm::vec3(1.f, 1.f, 1.f);
glm::vec3 lightPos = glm::vec3(0.f, 0.f, 5.f);

//Vertex array object and vertex buffer object indices 
GLuint lightVAO, lightVBO;


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
				this->v[i][j] *= 0.99;
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
					//printf("p1: (%.2f, %.2f, %.2f)\n", p1.x, p1.y, p1.z);
					//printf("p2: (%.2f, %.2f, %.2f)\n", p2.x, p2.y, p2.z);
					//printf("p3: (%.2f, %.2f, %.2f)\n", p3.x, p3.y, p3.z);
					//printf("p4: (%.2f, %.2f, %.2f)\n", p4.x, p4.y, p4.z);
					//printf("p5: (%.2f, %.2f, %.2f)\n", p5.x, p5.y, p5.z);
					glm::vec3 n1 = glm::cross(p4 - p1, p3 - p1);
					n1 = glm::normalize(n1);
					glm::vec3 n2 = glm::cross(p2 - p1, p4 - p1);
					n2 = glm::normalize(n2);
					glm::vec3 n3 = glm::cross(p5 - p1, p2 - p1);
					n3 = glm::normalize(n3);
					glm::vec3 n4 = glm::cross(p3 - p1, p5 - p1);
					n4 = glm::normalize(n4);

					normal = 0.25f * (n1 + n2 + n3 + n4);
					//printf("normal: (%.2f, %.2f, %.2f)\n", normal.x, normal.y, normal.z);

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


int CompileShaders() {
	//Vertex Shader
	const char* vsSrc = "#version 330 core\n"
		"layout (location = 0) in vec4 iPos;\n"
		"layout (location = 1) in vec3 aNormal;\n"
		"out vec3 FragPos;\n"
		"out vec3 Normal;\n"
		"uniform mat4 modelview;\n"
		"uniform mat4 model;\n"
		"uniform mat3 normalMat;\n"
		"void main()\n"
		"{\n"
		"   vec4 oPos=modelview * iPos;\n"
		"   gl_Position = vec4(oPos.x, oPos.y, oPos.z, oPos.w);\n"
		"   FragPos = vec3(model * iPos);\n"
		"   Normal = normalMat * normalize(aNormal);\n"
		"}\0";

	//Fragment Shader
	const char* fsSrc = "#version 330 core\n"
		"out vec4 col;\n"
		"in vec3 Normal;\n"
		"in vec3 FragPos;\n"
		"uniform vec3 color;\n"
		"uniform vec3 lightColor;\n"
		"uniform vec3 lightPos;\n"
		"uniform vec3 viewPos;\n"
		"void main()\n"
		"{\n"
		"	float ambientStrength = 0.1;\n"
		"	vec3 ambient = ambientStrength * lightColor;\n"

		"	float diffuseStrength = 0.8;\n"
		"	vec3 norm = normalize(Normal);\n"
		"	vec3 lightDir = normalize(lightPos - FragPos);\n"
		"	float diff = max(dot(norm, lightDir), 0.0);\n"
		"	vec3 diffuse = diffuseStrength * diff * lightColor;\n"

		"	float specularStrength = 0.5;\n"
		"	vec3 viewDir = normalize(viewPos - FragPos);\n"
		"	vec3 reflectDir = reflect(-lightDir, norm);\n"
		"	float spec = pow(max(dot(viewDir, reflectDir), 0.0), 32);\n"
		"	vec3 specular = specularStrength * spec * lightColor;\n"

		"	vec3 result = (ambient + diffuse + specular) * color;\n"
		"   col = vec4(result, 1.0);\n"
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

int CompileLightShaders() {
	//Vertex Shader
	const char* vsSrc = "#version 330 core\n"
		"layout (location = 0) in vec4 iPos;\n"
		"uniform mat4 modelview;\n"
		"void main()\n"
		"{\n"
		"   vec4 oPos=modelview*iPos;\n"
		"   gl_Position = vec4(oPos.x, oPos.y, oPos.z, oPos.w);\n"
		"}\0";

	//Fragment Shader
	const char* fsSrc = "#version 330 core\n"
		"out vec4 col;\n"
		"void main()\n"
		"{\n"
		"   col = vec4(1.0);\n"
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

void BuildLightScene(GLuint& VBO, GLuint& VAO) { //return VBO and VAO values n is the subdivision
	float vertices[] = {
		-0.5f, -0.5f, -0.5f,
		 0.5f, -0.5f, -0.5f,
		 0.5f,  0.5f, -0.5f,
		 0.5f,  0.5f, -0.5f,
		-0.5f,  0.5f, -0.5f,
		-0.5f, -0.5f, -0.5f,

		-0.5f, -0.5f,  0.5f,
		 0.5f, -0.5f,  0.5f,
		 0.5f,  0.5f,  0.5f,
		 0.5f,  0.5f,  0.5f,
		-0.5f,  0.5f,  0.5f,
		-0.5f, -0.5f,  0.5f,

		-0.5f,  0.5f,  0.5f,
		-0.5f,  0.5f, -0.5f,
		-0.5f, -0.5f, -0.5f,
		-0.5f, -0.5f, -0.5f,
		-0.5f, -0.5f,  0.5f,
		-0.5f,  0.5f,  0.5f,

		 0.5f,  0.5f,  0.5f,
		 0.5f,  0.5f, -0.5f,
		 0.5f, -0.5f, -0.5f,
		 0.5f, -0.5f, -0.5f,
		 0.5f, -0.5f,  0.5f,
		 0.5f,  0.5f,  0.5f,

		-0.5f, -0.5f, -0.5f,
		 0.5f, -0.5f, -0.5f,
		 0.5f, -0.5f,  0.5f,
		 0.5f, -0.5f,  0.5f,
		-0.5f, -0.5f,  0.5f,
		-0.5f, -0.5f, -0.5f,

		-0.5f,  0.5f, -0.5f,
		 0.5f,  0.5f, -0.5f,
		 0.5f,  0.5f,  0.5f,
		 0.5f,  0.5f,  0.5f,
		-0.5f,  0.5f,  0.5f,
		-0.5f,  0.5f, -0.5f,
	};

	//make VAO
	glGenVertexArrays(1, &VAO);
	glGenBuffers(1, &VBO);

	//bind it
	glBindVertexArray(VAO);

	//bind the VBO
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	//send the data to the GPU
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

	//Configure the attributes
	glVertexAttribPointer((GLuint)0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), 0);
	//Make it valid
	glEnableVertexAttribArray(0);
}


//Quit when ESC is released
static void KbdCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if (key == GLFW_KEY_ESCAPE && action == GLFW_RELEASE) glfwSetWindowShouldClose(window, GLFW_TRUE);
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
	if (window == NULL)
	{
		std::cout << "Cannot open GLFW window" << std::endl;
		glfwTerminate();
		return -1;
	}
	lastX = 500;
	lastY = 500;
	//Paste the window to the current context
	glfwMakeContextCurrent(window);

	//Load GLAD to configure OpenGL
	gladLoadGL();

	
	water_surface waterSurface;
	int shaderProg = CompileShaders();
	glUseProgram(shaderProg);
	GLint modelviewParameter = glGetUniformLocation(shaderProg, "modelview");
	GLint modelParameter = glGetUniformLocation(shaderProg, "model");
	GLint normalMatParameter = glGetUniformLocation(shaderProg, "normalMat");
	GLint lightPosParameter = glGetUniformLocation(shaderProg, "lightPos");
	GLint viewPosParameter = glGetUniformLocation(shaderProg, "viewPos");

	//Background color
	glClearColor(0.2f, 0.2f, 0.2f, 1.0f);
	//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glPointSize(pointSize);
	glLineWidth(lineWidth);

	BuildLightScene(lightVBO, lightVAO);
	int shaderProgLight = CompileLightShaders();


	// Initialize ImGUI
	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImGuiIO& io = ImGui::GetIO(); (void)io;
	ImGui::StyleColorsDark();
	ImGui_ImplGlfw_InitForOpenGL(window, true);
	ImGui_ImplOpenGL3_Init("#version 330");

	glfwSetKeyCallback(window, KbdCallback); //set keyboard callback to quit
	//glfwSetCursorPosCallback(window, mouse_callback);
	//glfwSetMouseButtonCallback(window, MouseButtonCallback);;
	//glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
	//glEnable(GL_DEPTH_TEST);

	// Main while loop
	while (!glfwWindowShouldClose(window))
	{
		
		float currentFrame = static_cast<float>(glfwGetTime());
		deltaTime = currentFrame - lastFrame;
		lastFrame = currentFrame;
		
		//Clean the window
		processInput(window);
		glClear(GL_COLOR_BUFFER_BIT);

		
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();

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
			}
			if (mouse_state == GLFW_PRESS || mouse_state == GLFW_REPEAT) {
				if (mouse_intersection.x > -3.0 &&
					mouse_intersection.x < 3.0 &&
					mouse_intersection.y > -3.0 &&
					mouse_intersection.y < 3.0) {
					int i = (mouse_intersection.x + 3.0) / 6.0 * waterSurface.width;
					int j = (mouse_intersection.y + 3.0) / 6.0 * waterSurface.height;

					if (i > 0 && j > 0 && i < waterSurface.width - 1 && j < waterSurface.height - 1) {
						waterSurface.u[i][j] = 0.7;
						waterSurface.u[i - 1][j - 1] = 0.4;
						waterSurface.u[i - 1][j] = 0.6;
						waterSurface.u[i - 1][j + 1] = 0.4;
						waterSurface.u[i + 1][j - 1] = 0.4;
						waterSurface.u[i + 1][j] = 0.6;
						waterSurface.u[i + 1][j + 1] = 0.4;
						waterSurface.u[i][j + 1] = 0.6;
						waterSurface.u[i][j - 1] = 0.6;
					}
				}
			}
			else if (mouse_state == GLFW_RELEASE && is_mouse_down) {
				is_mouse_down = false;
			}
		}

		
		waterSurface.update(deltaTime);

		glUseProgram(shaderProg);
		glBindVertexArray(waterSurface.vao);

		//send the color to the fragment shader
		glUniform3f(glGetUniformLocation(shaderProg, "color"), waterColor[0], waterColor[1], waterColor[2]);
		glUniform3f(glGetUniformLocation(shaderProg, "lightColor"), lightColor[0], lightColor[1], lightColor[2]);

		glm::mat4 model = glm::mat4(1.0f);
		glm::mat4 modelView = proj * view * model;
		glm::mat3 normalMat = glm::mat3(glm::transpose(glm::inverse(model)));
		glUniformMatrix4fv(modelviewParameter, 1, GL_FALSE, glm::value_ptr(modelView));
		glUniformMatrix4fv(modelParameter, 1, GL_FALSE, glm::value_ptr(model));
		glUniformMatrix3fv(normalMatParameter, 1, GL_FALSE, glm::value_ptr(normalMat));
		glUniform3f(lightPosParameter, lightPos[0], lightPos[1], lightPos[2]);
		glUniform3f(viewPosParameter, cameraPos[0], cameraPos[1], cameraPos[2]);

		glDrawElements(GL_TRIANGLES, (waterSurface.N - 1) * (waterSurface.N - 1) * 2 * 3, GL_UNSIGNED_INT, 0);

		glUseProgram(shaderProgLight);
		glBindVertexArray(lightVAO);

		model = glm::mat4();
		model = glm::translate(model, lightPos);
		model = glm::scale(model, glm::vec3(0.2f));
		modelView = proj * view * model;
		glUniformMatrix4fv(glGetUniformLocation(shaderProgLight, "modelview"), 1, GL_FALSE, glm::value_ptr(modelView));

		glDrawArrays(GL_TRIANGLES, 0, 36);

		// Renders the ImGUI elements
		ImGui::Render();
		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
		//Swap the back buffer with the front buffer
		glfwSwapBuffers(window);
		//make sure events are served
		glfwPollEvents();
	}
	//Cleanup
	glDeleteVertexArrays(1, &lightVAO);
	glDeleteBuffers(1, &lightVBO);
	glDeleteVertexArrays(1, &waterSurface.vao);
	glDeleteBuffers(1, &waterSurface.elements_vbo);
	glDeleteProgram(shaderProg);
	glDeleteProgram(shaderProgLight);
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