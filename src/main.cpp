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
#include <algorithm>

#include "triangle.h"  //triangles
#include "helper.h"         
#include "objGen.h"         //to save OBJ file format for 3D printing
#include "trackball.h"
#include "shaders.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include <cstdlib>


#pragma warning(disable : 4996)
#pragma comment(lib, "glfw3.lib")


using namespace std;
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void processInput(GLFWwindow* window);
static void KbdCallback(GLFWwindow* window, int key, int scancode, int action, int mods);

glm::vec3 camera_pos = glm::vec3(5.f, 2.f, 1.f);
glm::vec3 camera_front = glm::vec3(-5.0f, -2.0f, -1.f);
glm::vec3 camera_up = glm::vec3(0.0f, 0.0f, 1.0f);
float delta_time = 0.0f;	// time between current frame and last frame
float last_frame = 0.0f;
float PI = 3.14159;

float player_speed = 2.5f;
int point_size = 1;
int line_width = 1;
bool is_first_mouse = true;
float yaw = -90.0f;	// yaw is initialized to -90.0 degrees since a yaw of 0.0 results in a direction vector pointing to the right so we initially rotate a bit to the left.
float pitch = 0.0f;
float last_x = 800 / 2.0;
float last_y = 600 / 2.0;
float fov = 90.0f;

float water_color[] = { 0.34, 0.75, 0.87 };
glm::vec3 light_color = glm::vec3(1.f, 1.f, 1.f);
float light_pos[] = { 0.f, 0.f, 10.f };

//Vertex array object and vertex buffer object indices 
GLuint light_vao, light_vbo;
GLuint floor_vao, floor_vbo;
GLuint floor_texture;


class WaterSurface {
private:
	GLuint points_vbo, normals_vbo, elements_vbo, caustic_vbo, tex_coords_vbo;

	void pointsCalculate() {
		int point_index = 0;
		for (int x = 0; x < N; x++) {
			//printf("u: (%.2f)\n", x);
			for (int y = 0; y < M; y++) {
				float i = -this->water_width / 2 + this->water_width * ((x / (float)this->width));
				float j = -this->water_length / 2 + this->water_length * ((y / (float)this->length));

				// Take the average height of the blocks around us, weighted by their distances
				float sum_of_weights = 0.0;
				float avg_height = 0.0;
				for (int k = x - 1; k <= x + 1; k++) {
					if (k < 0 || k >= this->width) {
						continue;
					}

					for (int l = y - 1; l <= y + 1; l++) {
						if (l < 0 || l >= this->length) {
							continue;
						}

						float u = -this->water_width / 2 + this->water_width * (k / (float)this->width);
						float v = -this->water_length / 2 + this->water_length * (l / (float)this->length);

						// Make sure the weight is larger for smaller distances
						float weight = 100 - (u - i) * (u - i) + (j - v) * (j - v);
						avg_height += this->u[k][l] * weight;
						sum_of_weights += weight;
					}
				}
				avg_height /= sum_of_weights;

				points[point_index++] = glm::vec3(i, j, avg_height);
			}
		}
	}

	void normalsCalculate() {
		for (int i = 0; i < N; i++) {
			//printf("u: (%i)\n", i);
			for (int j = 0; j < M; j++) {
				//printf("u: (%i)\n", j);
				// Average the normals for each triangle around us
				int p1_i = i * M + j;
				int p2_i = (i + 1) * M + j;
				int p3_i = (i - 1) * M + j;
				int p4_i = i * M + (j - 1);
				int p5_i = i * M + (j + 1);

				glm::vec3 normal;

				if (i == 0 && j == 0) {
					glm::vec3 p1 = points[p1_i];
					glm::vec3 p2 = points[p2_i];
					glm::vec3 p5 = points[p5_i];

					normal = glm::cross(p5 - p1, p2 - p1);
				}
				else if (i == 0 && j == M - 1) {
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
				else if (i == N - 1 && j == M - 1) {
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
				else if (j == M - 1) {
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

				normal = -glm::normalize(normal);
				normals[p1_i] = normal;

			}
		}
	}

	void elementsCalculate() {
		int e_i = 0;
		for (int i = 0; i < N - 1; i++) {
			for (int j = 0; j < M - 1; j++) {
				// First triangle
				int p1 = i * M + j;
				int p2 = (i + 1) * M + j;
				int p3 = i * M + (j + 1);

				elements[e_i++] = glm::vec3(p1, p2, p3);

				// Second triangle
				int p4 = (i + 1) * M + j;
				int p5 = (i + 1) * M + (j + 1);
				int p6 = i * M + (j + 1);

				elements[e_i++] = glm::vec3(p4, p5, p6);
			}
		}

	}

	void causticCalculate() {
		int caustic_index = 0;
		glm::vec3 light_pos_vec = glm::vec3(light_pos[0], light_pos[1], light_pos[2]);
		//waterSurface.points;
		for (int i = 0; i < N * M; i++) {
			glm::vec3 water_pos = points[i];
			glm::vec3 incident = water_pos - light_pos_vec;
			glm::vec3 water_normal = normals[i];

			glm::vec3 refraction = glm::normalize(glm::refract(incident, water_normal, 0.7));
			float constant = abs((water_pos.z + caustic_depth) / refraction.z);
			glm::vec3 caustic_pos = water_pos + refraction * constant;
			if (caustic_pos.x < -water_width / 2 || caustic_pos.x > water_width / 2 || caustic_pos.y < -water_length / 2 || caustic_pos.y > water_length / 2) {
				float minc = 1;
				if (caustic_pos.x < 0 && caustic_pos.y < 0) {
					minc = max(abs(caustic_pos.x - water_pos.x) / abs(-water_width / 2 - water_pos.x), abs(caustic_pos.y - water_pos.y) / abs(-water_length / 2 - water_pos.y));
				}
				else if (caustic_pos.x < 0 && caustic_pos.y >= 0) {
					minc = max(abs(caustic_pos.x - water_pos.x) / abs(-water_width / 2 - water_pos.x), abs(caustic_pos.y - water_pos.y) / abs(water_length / 2 - water_pos.y));
				}
				else if (caustic_pos.x >= 0 && caustic_pos.y >= 0) {
					minc = max(abs(caustic_pos.x - water_pos.x) / abs(water_width / 2 - water_pos.x), abs(caustic_pos.y - water_pos.y) / abs(water_length / 2 - water_pos.y));
				}
				else {
					minc = max(abs(caustic_pos.x - water_pos.x) / abs(water_width / 2 - water_pos.x), abs(caustic_pos.y - water_pos.y) / abs(-water_length / 2 - water_pos.y));
				}
				caustic_pos = water_pos + refraction * constant / minc;
			}
			caustic[caustic_index++] = caustic_pos;
		}

	}

	void texCoordsCalculate() {
		int tex_coord_index = 0;
		//waterSurface.points;
		for (int i = 0; i < N * M; i++) {
			glm::vec3 water_pos = points[i];
			float ratio = 0.7;
			if (water_pos.z > camera_pos.z) {
				tex_coords[tex_coord_index++] = glm::vec2(0.f, 0.f);
				continue;
			}
			glm::vec3 incident = water_pos - camera_pos;
			glm::vec3 water_normal = normals[i];

			glm::vec3 refraction = glm::normalize(glm::refract(incident, water_normal, ratio));
			float constant = abs((water_pos.z + caustic_depth) / refraction.z);
			glm::vec3 floor_hit_pos = water_pos + refraction * constant;
			glm::vec2 tex_coord = glm::vec2(
				(floor_hit_pos.x + water_width / 2) / water_width,
				(floor_hit_pos.y + water_length / 2) / water_length
			);
			if (floor_hit_pos.x <= -water_width / 2 || floor_hit_pos.x >= water_width / 2 || floor_hit_pos.y <= -water_length / 2 || floor_hit_pos.y >= water_length / 2) {
				float minc = 1;
				if (floor_hit_pos.x < 0 && floor_hit_pos.y < 0) {
					float x_const = abs(floor_hit_pos.x - water_pos.x) / abs(-water_width / 2 - water_pos.x);
					float y_const = abs(floor_hit_pos.y - water_pos.y) / abs(-water_length / 2 - water_pos.y);
					minc = max(x_const, y_const);
					floor_hit_pos = water_pos + refraction * constant / minc;
					if (x_const > y_const) {
						tex_coord = glm::vec2(
							(caustic_depth * (depth_overflow - 1) - floor_hit_pos.z) / (caustic_depth * depth_overflow),
							(floor_hit_pos.y + water_length / 2) / water_length
						);
					}
					else {
						tex_coord = glm::vec2(
							(floor_hit_pos.x + water_width / 2) / water_width,
							(floor_hit_pos.z + caustic_depth) / caustic_depth / depth_overflow
						);
					}
				}
				else if (floor_hit_pos.x < 0 && floor_hit_pos.y >= 0) {
					minc = max(abs(floor_hit_pos.x - water_pos.x) / abs(-water_width / 2 - water_pos.x), abs(floor_hit_pos.y - water_pos.y) / abs(water_length / 2 - water_pos.y));
					float x_const = abs(floor_hit_pos.x - water_pos.x) / abs(-water_width / 2 - water_pos.x);
					float y_const = abs(floor_hit_pos.y - water_pos.y) / abs(water_length / 2 - water_pos.y);
					minc = max(x_const, y_const);
					floor_hit_pos = water_pos + refraction * constant / minc;
					if (x_const > y_const) {
						tex_coord = glm::vec2(
							//(floor_hit_pos.z + caustic_depth) / caustic_depth / depth_overflow,
							(caustic_depth * (depth_overflow - 1) - floor_hit_pos.z) / (caustic_depth * depth_overflow),
							(floor_hit_pos.y + water_length / 2) / water_length
						);
						
					}
					else {
						tex_coord = glm::vec2(
							(floor_hit_pos.x + water_width / 2) / water_width,
							(floor_hit_pos.z + caustic_depth) / caustic_depth / depth_overflow
						);
					}
					
				}
				else if (floor_hit_pos.x >= 0 && floor_hit_pos.y >= 0) {
					minc = max(abs(floor_hit_pos.x - water_pos.x) / abs(water_width / 2 - water_pos.x), abs(floor_hit_pos.y - water_pos.y) / abs(water_length / 2 - water_pos.y));
					float x_const = abs(floor_hit_pos.x - water_pos.x) / abs(water_width / 2 - water_pos.x);
					float y_const = abs(floor_hit_pos.y - water_pos.y) / abs(water_length / 2 - water_pos.y);
					minc = max(x_const, y_const);
					floor_hit_pos = water_pos + refraction * constant / minc;
					if (x_const > y_const) {
						tex_coord = glm::vec2(
							//(floor_hit_pos.z + caustic_depth) / caustic_depth / depth_overflow,
							(caustic_depth * (depth_overflow - 1) - floor_hit_pos.z) / (caustic_depth * depth_overflow),
							(floor_hit_pos.y + water_length / 2) / water_length
						);
					}
					else {
						tex_coord = glm::vec2(
							(floor_hit_pos.x + water_width / 2) / water_width,
							(floor_hit_pos.z + caustic_depth) / caustic_depth / depth_overflow
						);
					}
				}
				else {
					minc = max(abs(floor_hit_pos.x - water_pos.x) / abs(water_width / 2 - water_pos.x), abs(floor_hit_pos.y - water_pos.y) / abs(-water_length / 2 - water_pos.y));
					float x_const = abs(floor_hit_pos.x - water_pos.x) / abs(water_width / 2 - water_pos.x);
					float y_const = abs(floor_hit_pos.y - water_pos.y) / abs(-water_length / 2 - water_pos.y);
					minc = max(x_const, y_const);
					floor_hit_pos = water_pos + refraction * constant / minc;
					if (x_const > y_const) {
						tex_coord = glm::vec2(
							//(floor_hit_pos.z + caustic_depth) / caustic_depth / depth_overflow,
							(caustic_depth * (depth_overflow - 1) - floor_hit_pos.z) / (caustic_depth * depth_overflow),
							(floor_hit_pos.y + water_length / 2) / water_length
						);
					}
					else {
						tex_coord = glm::vec2(
							(floor_hit_pos.x + water_width / 2) / water_width,
							(floor_hit_pos.z + caustic_depth) / caustic_depth / depth_overflow
						);
					}
				}
			}

			//if (abs(floor_hit_pos.z + caustic_depth) < 0.007) {
			//	printf("caustic: (%.2f, %.2f, %.2f), tex_coord: (%.2f, %.2f)\n", floor_hit_pos[0], floor_hit_pos[1], floor_hit_pos[2], tex_coord[0], tex_coord[1]);
			//}
			tex_coords[tex_coord_index++] = tex_coord;
		}

	}

	void buildScene() {
		glBindVertexArray(vao);

		// Set up the points vbo
		for (int i = 0; i < N * M; i++) {
			points_buffer[i * 3] = points[i].x;
			points_buffer[i * 3 + 1] = points[i].y;
			points_buffer[i * 3 + 2] = points[i].z;
		}
		glBindBuffer(GL_ARRAY_BUFFER, points_vbo);
		glBufferData(GL_ARRAY_BUFFER, sizeof(float) * N * M * 3, points_buffer, GL_STATIC_DRAW);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);
		glEnableVertexAttribArray(0);

		// Set up the normals vbo
		for (int i = 0; i < N * M; i++) {
			normals_buffer[i * 3] = normals[i].x;
			normals_buffer[i * 3 + 1] = normals[i].y;
			normals_buffer[i * 3 + 2] = normals[i].z;
		}
		glBindBuffer(GL_ARRAY_BUFFER, normals_vbo);
		glBufferData(GL_ARRAY_BUFFER, sizeof(float) * N * M * 3, normals_buffer, GL_STATIC_DRAW);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, NULL);
		glEnableVertexAttribArray(1);

		// Set up the texture coordinates vbo
		for (int i = 0; i < N * M; i++) {
			tex_coords_buffer[i * 2] = tex_coords[i].x;
			tex_coords_buffer[i * 2 + 1] = tex_coords[i].y;
		}
		glBindBuffer(GL_ARRAY_BUFFER, tex_coords_vbo);
		glBufferData(GL_ARRAY_BUFFER, sizeof(float) * N * M * 2, tex_coords_buffer, GL_STATIC_DRAW);
		glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 0, NULL);
		glEnableVertexAttribArray(2);

		// Set up elements vbo
		for (int i = 0; i < (N - 1) * (M - 1) * 2; i++) {
			elements_buffer[i * 3] = elements[i].x;
			elements_buffer[i * 3 + 1] = elements[i].y;
			elements_buffer[i * 3 + 2] = elements[i].z;
		}
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, elements_vbo);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(int) * (N - 1) * (M - 1) * 2 * 3, elements_buffer, GL_STATIC_DRAW);
		glBindVertexArray(0);

		glBindVertexArray(caustic_vao);
		// Set up the caustic vbo
		for (int i = 0; i < N * M; i++) {
			caustic_buffer[i * 3] = caustic[i].x;
			caustic_buffer[i * 3 + 1] = caustic[i].y;
			caustic_buffer[i * 3 + 2] = caustic[i].z;
		}

		glBindBuffer(GL_ARRAY_BUFFER, points_vbo);
		glBufferData(GL_ARRAY_BUFFER, sizeof(float) * N * M * 3, points_buffer, GL_STATIC_DRAW);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);
		glEnableVertexAttribArray(0);

		glBindBuffer(GL_ARRAY_BUFFER, caustic_vbo);
		glBufferData(GL_ARRAY_BUFFER, sizeof(float) * N * M * 3, caustic_buffer, GL_STATIC_DRAW);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, NULL);
		glEnableVertexAttribArray(1);


		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, elements_vbo);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(int) * (N - 1) * (M - 1) * 2 * 3, elements_buffer, GL_STATIC_DRAW);

		glBindVertexArray(0);
	}

public:
	static constexpr int width = 60, length = 60;
	static constexpr int N = width + 1, M = length + 1;

	float water_width = (float)width / 10, water_length = (float)length / 10;
	float caustic_depth = 3.f;
	float depth_overflow = 1.5;

	float c = 16.0;
	float damp = 0.99;
	float u[width][length];
	float v[width][length];
	//float u_new[width][height];
	//float control_point_heights[width][height];

	GLuint vao, caustic_vao;

	//glm::vec3 points[N * M];
	//glm::vec3 normals[N * M];
	//glm::vec3 caustic[N * M];
	//glm::vec3 elements[(N - 1) * (M - 1) * 2];

	//float points_buffer[N * M * 3];
	//float normals_buffer[N * M * 3];
	//float caustic_buffer[N * M * 3];
	//int elements_buffer[(N - 1) * (M - 1) * 2 * 3];

	float** u_new = (float**)malloc(sizeof(float*) * width);
	float** control_point_heights = (float**)malloc(sizeof(float*) * width);

	glm::vec3* points = (glm::vec3*)malloc(sizeof(glm::vec3) * N * M);
	glm::vec3* normals = (glm::vec3*)malloc(sizeof(glm::vec3) * N * M);
	glm::vec3* caustic = (glm::vec3*)malloc(sizeof(glm::vec3) * N * M);
	glm::vec2* tex_coords = (glm::vec2*)malloc(sizeof(glm::vec2) * N * M);
	glm::vec3* elements = (glm::vec3*)malloc(sizeof(glm::vec3) * (N - 1) * (M - 1) * 2);

	float* points_buffer = (float*)malloc(sizeof(float) * N * M * 3);
	float* normals_buffer = (float*)malloc(sizeof(float) * N * M * 3);
	float* caustic_buffer = (float*)malloc(sizeof(float) * N * M * 3);
	float* tex_coords_buffer = (float*)malloc(sizeof(float) * N * M * 2);
	int* elements_buffer = (int*)malloc(sizeof(int) * (N - 1) * (M - 1) * 2 * 3);


	WaterSurface() {
		glGenVertexArrays(1, &vao);
		glGenVertexArrays(1, &caustic_vao);
		glGenBuffers(1, &elements_vbo);
		glGenBuffers(1, &normals_vbo);
		glGenBuffers(1, &points_vbo);
		glGenBuffers(1, &caustic_vbo);
		glGenBuffers(1, &tex_coords_vbo);

		// Initialize water heights
		for (int i = 0; i < this->width; i++) {
			for (int j = 0; j < this->length; j++) {
				this->u[i][j] = 0;
				this->v[i][j] = 0;
			}
			u_new[i] = (float*)malloc(sizeof(float) * length);
			control_point_heights[i] = (float*)malloc(sizeof(float) * length);
		}
	}

	~WaterSurface() {
		for (int i = 0; i < width; i++) {
			free(u_new[i]);
			free(control_point_heights[i]);
		}
		free(u_new);
		free(control_point_heights);
		free(points);
		free(normals);
		free(caustic);
		free(tex_coords);
		free(elements);
		free(points_buffer);
		free(normals_buffer);
		free(caustic_buffer);
		free(elements_buffer);
		free(tex_coords_buffer);
		glDeleteVertexArrays(1, &this->vao);
		glDeleteVertexArrays(1, &this->caustic_vao);
		glDeleteBuffers(1, &this->elements_vbo);
		glDeleteBuffers(1, &this->normals_vbo);
		glDeleteBuffers(1, &this->points_vbo);
		glDeleteBuffers(1, &this->caustic_vbo);
		glDeleteBuffers(1, &this->tex_coords_vbo);
	}

	void poke(glm::vec3 pos, float rippleHeight) {
		if (pos.x > -this->water_width / 2 &&
			pos.x < this->water_width / 2 &&
			pos.y > -this->water_length / 2 &&
			pos.y < this->water_length / 2) {
			int i = ((pos.x + this->water_width / 2) * this->width / this->water_width);
			int j = ((pos.y + this->water_length / 2) * this->length / this->water_length);

			if (i > 0 && j > 0 && i < this->width - 1 && j < this->length - 1) {
				this->u[i][j] = rippleHeight;
				this->u[i - 1][j - 1] = rippleHeight * 0.7;
				this->u[i - 1][j] = rippleHeight * 0.85;
				this->u[i - 1][j + 1] = rippleHeight * 0.7;
				this->u[i + 1][j - 1] = rippleHeight * 0.7;
				this->u[i + 1][j] = rippleHeight * 0.85;
				this->u[i + 1][j + 1] = rippleHeight * 0.7;
				this->u[i][j + 1] = rippleHeight * 0.85;
				this->u[i][j - 1] = rippleHeight * 0.85;
			}
		}
	}

	void update(float dt) {
		for (int i = 0; i < this->width; i++) {
			for (int j = 0; j < this->length; j++) {
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

				if (j == this->length - 1) {
					v4 = this->u[i][j];
				}
				else {
					v4 = this->u[i][j + 1];
				}

				float f = c * c * ((v1 + v2 + v3 + v4) - 4 * this->u[i][j]);
				this->v[i][j] += f * dt;
				this->v[i][j] *= this->damp;
				this->u_new[i][j] = u[i][j] + v[i][j] * dt;
			}
		}

		double sum_of_u = 0.0;
		for (int i = 0; i < this->width; i++) {
			for (int j = 0; j < this->length; j++) {
				sum_of_u += this->u_new[i][j];
			}
		}

		double avg_of_u = sum_of_u / (this->width * this->length);

		for (int i = 0; i < this->width; i++) {
			for (int j = 0; j < this->length; j++) {
				this->u[i][j] = this->u_new[i][j] - avg_of_u;
				this->control_point_heights[i][j] = this->u[i][j];
			}
		}

		static int x_delta[9] = { 0, -1, -1, 0, 1, 1, 1, 0, -1 };
		static int y_delta[9] = { 0, 0, -1, -1, -1, 0, 1, 1, 1 };

		for (int i = 3; i < this->width; i += 3) {
			for (int j = 3; j < this->length; j += 3) {
				glm::vec3 points[9];

				for (int k = 0; k < 9; k++) {
					int x_index = i + x_delta[k];
					int y_index = j + y_delta[k];

					points[k].x = -this->water_width / 2 + this->water_width * (1 - (x_index / (float)this->width));
					points[k].y = -this->water_length / 2 + this->water_length * (1 - (y_index / (float)this->length));
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
		this->pointsCalculate();
		// Calculate normals
		this->normalsCalculate();
		// Calculate the elements for each triangle
		this->elementsCalculate();
		// Calculate the caustic
		this->causticCalculate();
		// Calculate the texture of surface
		this->texCoordsCalculate();
		// Build scene
		this->buildScene();
	}
};

void InitShaders(GLuint* program, const char* vsSrc, const char* fsSrc) {
	std::vector<GLuint> shaderList;

	//load and compile shaders 	
	shaderList.push_back(CreateShader(GL_VERTEX_SHADER, LoadShader(vsSrc)));
	shaderList.push_back(CreateShader(GL_FRAGMENT_SHADER, LoadShader(fsSrc)));

	//create the shader program and attach the shaders to it
	*program = CreateProgram(shaderList);

	//delete shaders (they are on the GPU now)
	using namespace std;
	for_each(shaderList.begin(), shaderList.end(), glDeleteShader);
}


void BuildFloorScene(GLuint& VBO, GLuint& VAO, GLuint& Texture, char* file_name) {
	float vertices[] = {
		// positions          // normals		// texture coords
		 0.5f,  0.5f, 0.0f,   0.f, 0.f, 1.f,	1.0f, 1.0f, // top right
		 0.5f, -0.5f, 0.0f,   0.f, 0.f, 1.f,	1.0f, 0.0f, // bottom right
		-0.5f, -0.5f, 0.0f,   0.f, 0.f, 1.f,	0.0f, 0.0f, // bottom left
		-0.5f,  0.5f, 0.0f,   0.f, 0.f, 1.f,	0.0f, 1.0f  // top left 
	};
	unsigned int indices[] = {
		0, 1, 3, // first triangle
		1, 2, 3  // second triangle
	};
	unsigned int EBO;
	glGenVertexArrays(1, &VAO);
	glGenBuffers(1, &VBO);
	glGenBuffers(1, &EBO);

	glBindVertexArray(VAO);

	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

	// position attribute
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);
	// color attribute
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(3 * sizeof(float)));
	glEnableVertexAttribArray(1);
	// texture coord attribute
	glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(6 * sizeof(float)));
	glEnableVertexAttribArray(2);

	glBindVertexArray(0);

	// Bind Texture
	stbi_set_flip_vertically_on_load(true);
	int x, y, n;
	int force_channels = 3;
	unsigned char* image_data = stbi_load(file_name, &x, &y, &n, force_channels);

	if (!image_data) {
		printf("Could not load texture file %s\n", file_name);
		return;
	}

	glGenTextures(1, &Texture);
	glBindTexture(GL_TEXTURE_2D, Texture);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, x, y, 0, GL_RGB, GL_UNSIGNED_BYTE, image_data);
	glGenerateMipmap(GL_TEXTURE_2D);
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glBindTexture(GL_TEXTURE_2D, 0);

	stbi_image_free(image_data);
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
	GLFWwindow* window = glfwCreateWindow(1000, 1000, "Simple", NULL, NULL);
	//is all OK?
	if (window == NULL)
	{
		std::cout << "Cannot open GLFW window" << std::endl;
		glfwTerminate();
		return -1;
	}
	last_x = 500;
	last_y = 500;
	//Paste the window to the current context
	glfwMakeContextCurrent(window);

	//Load GLAD to configure OpenGL
	gladLoadGL();

	// waterSurface cannot be moved out of main function
	WaterSurface water_surface;
	GLuint shaderProg;
	InitShaders(&shaderProg, (char*)"shaders/water.vert", (char*)"shaders/water.frag");
	glUseProgram(shaderProg);
	GLint modelviewParameter = glGetUniformLocation(shaderProg, "modelview");
	GLint modelParameter = glGetUniformLocation(shaderProg, "model");
	GLint normalMatParameter = glGetUniformLocation(shaderProg, "normalMat");
	GLint lightPosParameter = glGetUniformLocation(shaderProg, "lightPos");
	GLint viewPosParameter = glGetUniformLocation(shaderProg, "viewPos");

	GLuint shaderProgCaustic;
	InitShaders(&shaderProgCaustic, (char*)"shaders/caustic.vert", (char*)"shaders/caustic.frag");
	glUseProgram(shaderProgCaustic);
	GLint causticModelviewParameter = glGetUniformLocation(shaderProgCaustic, "modelview");
	GLint causticModelParameter = glGetUniformLocation(shaderProgCaustic, "model");
	//GLint causticLightPosParameter = glGetUniformLocation(shaderProgCaustic, "lightPos");
	//GLint causticDepthParameter = glGetUniformLocation(shaderProgCaustic, "causticDepth");

	BuildFloorScene(floor_vbo, floor_vao, floor_texture, (char*)"floor4.jpg");
	GLuint shaderProgFloor;
	InitShaders(&shaderProgFloor, (char*)"shaders/floor.vert", (char*)"shaders/floor.frag");
	
	BuildLightScene(light_vbo, light_vao);
	GLuint shaderProgLight;
	InitShaders(&shaderProgLight, (char*)"shaders/light.vert", (char*)"shaders/light.frag");

	//Background color
	glClearColor(0.2f, 0.2f, 0.2f, 1.0f);
	//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glPointSize(point_size);
	glLineWidth(line_width);


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
	glEnable(GL_BLEND);

	float ripple_height = 0.7f;
	int raindrops_freq = 10;
	bool is_first_frame = true;

	float floorAmbientStrength = 0.65f;
	float floorDiffuseStrength = 0.3f;

	float causticsStrength = 0.99f;
	float causticsAlpha = 0.4f;

	float waterAmbientStrength = 0.2f;
	float waterDiffuseStrength = 0.6f;
	float waterTextureStrength = 1.f;
	float waterColorStrength = 0.6f;
	float waterAlpha = 0.65f;

	// Main while loop
	while (!glfwWindowShouldClose(window))
	{
		float current_frame = static_cast<float>(glfwGetTime());
		delta_time = current_frame - last_frame;
		last_frame = current_frame;

		//Clean the window
		processInput(window);
		glClear(GL_COLOR_BUFFER_BIT);
		//glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();
		ImGui::Begin("Water Surface");
		ImGui::SliderFloat("Diffusion Rate", &water_surface.c, 5.0f, 20.0f);
		ImGui::SliderFloat("Damp Rate", &water_surface.damp, 0.9f, 0.9999f);
		ImGui::SliderFloat("Ripple Height", &ripple_height, 0.1f, 1.2f);

		ImGui::SliderInt("Rain Intensity", &raindrops_freq, 1, 80);

		ImGui::ColorEdit3("Water Color", water_color);
		ImGui::SliderFloat("Light Position", &light_pos[2], 3.f, 10.0f);

		ImGui::SliderFloat("caustic strength", &causticsStrength, 0.1f, 1.f);
		ImGui::SliderFloat("caustic alpha", &causticsAlpha, 0.1f, 1.f);
		ImGui::SliderFloat("floor ambient strength", &floorAmbientStrength, 0.1f, 1.f);
		ImGui::SliderFloat("floor diffuse strength", &floorDiffuseStrength, 0.1f, 1.f);
		ImGui::SliderFloat("waterAmbientStrength", &waterAmbientStrength, 0.1f, 1.f);
		ImGui::SliderFloat("waterDiffuseStrength", &waterDiffuseStrength, 0.1f, 1.f);
		ImGui::SliderFloat("waterTextureStrength", &waterTextureStrength, 0.1f, 1.f);
		ImGui::SliderFloat("waterColorStrength", &waterColorStrength, 0.1f, 1.f);
		ImGui::SliderFloat("waterAlpha", &waterAlpha, 0.1f, 1.f);


		ImGui::End();
		//set the projection matrix
		glm::mat4 proj = glm::perspective(fov, 1.0f, 0.1f, 100.0f);
		//set the viewing matrix (looking from [0,0,5] to [0,0,0])
		glm::mat4 view = glm::lookAt(camera_pos, camera_pos + camera_front, camera_up);

		// Process poking on the water 
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
		glm::vec3 mouse_intersection = camera_pos + (-camera_pos.z / mouse_world.z) * mouse_world;

		ImGuiIO& io = ImGui::GetIO();
		io.AddMouseButtonEvent(GLFW_MOUSE_BUTTON_LEFT, glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT));
		if (!io.WantCaptureMouse)  //make sure you do not call this callback when over a menu
		{
			int mouse_state = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
			if (mouse_state == GLFW_PRESS && !is_mouse_down) {
				is_mouse_down = true;
			}
			if (mouse_state == GLFW_PRESS || mouse_state == GLFW_REPEAT) {
				water_surface.poke(mouse_intersection, ripple_height);
			}
			else if (mouse_state == GLFW_RELEASE && is_mouse_down) {
				is_mouse_down = false;
			}
		}

		if (!is_first_frame) {
			float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
			float rain_prob = min(delta_time, 0.016f) * raindrops_freq; // limit the delta time to be 0.05s at minimum (20fps)
			if (r < rain_prob) {
				double randx = (double)rand() / (RAND_MAX + 1) * water_surface.water_width - water_surface.water_width / 2;
				double randy = (double)rand() / (RAND_MAX + 1) * water_surface.water_length - water_surface.water_length / 2;
				water_surface.poke(glm::vec3(randx, randy, 0), ripple_height);
			}
		}
		is_first_frame = false;

		water_surface.update(delta_time);

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, floor_texture);

		// Draw walls and floor
		glUseProgram(shaderProgFloor);
		glBindVertexArray(floor_vao);

		glUniform3f(glGetUniformLocation(shaderProgFloor, "color"), water_color[0], water_color[1], water_color[2]);
		glUniform3f(glGetUniformLocation(shaderProgFloor, "lightColor"), light_color[0], light_color[1], light_color[2]);
		glUniform3f(glGetUniformLocation(shaderProgFloor, "lightPos"), light_pos[0], light_pos[1], light_pos[2]);
		glUniform1i(glGetUniformLocation(shaderProgFloor, "floorTexture"), 0);
		glUniform1f(glGetUniformLocation(shaderProgFloor, "ambientStrength"), floorAmbientStrength);
		glUniform1f(glGetUniformLocation(shaderProgFloor, "diffuseStrength"), floorDiffuseStrength);

		// bottom, right, left, back, front

		glm::vec3 scale[] = {
			glm::vec3(water_surface.water_width, water_surface.water_length, 1),
			glm::vec3(water_surface.water_width, water_surface.caustic_depth * water_surface.depth_overflow, 1),
			glm::vec3(water_surface.water_width, water_surface.caustic_depth * water_surface.depth_overflow, 1),
			glm::vec3(water_surface.caustic_depth * water_surface.depth_overflow, water_surface.water_length, 1),
			glm::vec3(water_surface.caustic_depth * water_surface.depth_overflow, water_surface.water_length, 1)
		};

		glm::vec4 rotation[] = {
			glm::vec4(0.f, 0.f, 0.f, 1.f),
			glm::vec4(90.f, 1.f, 0.f, 0.f),
			glm::vec4(90.f, 1.f, 0.f, 0.f),
			glm::vec4(90.f, 0.f, 1.f, 0.f),
			glm::vec4(90.f, 0.f, 1.f, 0.f),
		};

		glm::vec3 translation[] = {
			glm::vec3(0.f, 0.f, -water_surface.caustic_depth),
			glm::vec3(0.f, water_surface.water_length / 2, -water_surface.caustic_depth * (1 - water_surface.depth_overflow / 2)),
			glm::vec3(0.f, -water_surface.water_length / 2, -water_surface.caustic_depth * (1 - water_surface.depth_overflow / 2)),
			glm::vec3(-water_surface.water_width / 2, 0.f, -water_surface.caustic_depth * (1 - water_surface.depth_overflow / 2)),
			glm::vec3(water_surface.water_width / 2, 0.f, -water_surface.caustic_depth * (1 - water_surface.depth_overflow / 2)),
		};

		for (int i = 0; i < 4; i++) {
			glm::mat4 model = glm::mat4();
			model = glm::translate(model, translation[i]);
			model = glm::rotate(model, rotation[i][0], glm::vec3(rotation[i][1], rotation[i][2], rotation[i][3]));
			model = glm::scale(model, scale[i]);

			glm::mat4 model_view = proj * view * model;
			glm::mat3 normal_mat = glm::mat3(glm::transpose(glm::inverse(model)));

			glUniformMatrix4fv(glGetUniformLocation(shaderProgFloor, "modelview"), 1, GL_FALSE, glm::value_ptr(model_view));
			glUniformMatrix4fv(glGetUniformLocation(shaderProgFloor, "model"), 1, GL_FALSE, glm::value_ptr(model));
			glUniformMatrix3fv(glGetUniformLocation(shaderProgFloor, "normalMat"), 1, GL_FALSE, glm::value_ptr(normal_mat));

			glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
		}

		
		// Draw water and caustic
		glm::mat4 model = glm::mat4(1.0f);
		glm::mat4 model_view = proj * view * model;
		glm::mat3 normal_mat = glm::mat3(glm::transpose(glm::inverse(model)));

		// Caustic
		glUseProgram(shaderProgCaustic);
		glBindVertexArray(water_surface.caustic_vao);


		glUniform3f(glGetUniformLocation(shaderProgCaustic, "color"), water_color[0], water_color[1], water_color[2]);
		glUniform1f(glGetUniformLocation(shaderProgCaustic, "causticsStrength"), causticsStrength);
		glUniform1f(glGetUniformLocation(shaderProgCaustic, "causticsAlpha"), causticsAlpha);

		glUniformMatrix4fv(glGetUniformLocation(shaderProgCaustic, "modelview"), 1, GL_FALSE, glm::value_ptr(model_view));
		glUniformMatrix4fv(glGetUniformLocation(shaderProgCaustic, "model"), 1, GL_FALSE, glm::value_ptr(model));
		//glUniform3f(causticLightPosParameter, lightPos[0], lightPos[1], lightPos[2]);
		//glUniform1f(causticDepthParameter, waterSurface.causticDepth);

		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glDrawElements(GL_TRIANGLES, (water_surface.N - 1) * (water_surface.M - 1) * 2 * 3, GL_UNSIGNED_INT, 0);

		// Water
		glUseProgram(shaderProg);
		glBindVertexArray(water_surface.vao);

		//send the color to the fragment shader
		glUniform3f(glGetUniformLocation(shaderProg, "color"), water_color[0], water_color[1], water_color[2]);
		glUniform3f(glGetUniformLocation(shaderProg, "lightColor"), light_color[0], light_color[1], light_color[2]);
		glUniform1i(glGetUniformLocation(shaderProg, "floorTexture"), 0);
		glUniform1f(glGetUniformLocation(shaderProg, "ambientStrength"), waterAmbientStrength);
		glUniform1f(glGetUniformLocation(shaderProg, "diffuseStrength"), waterDiffuseStrength);
		glUniform1f(glGetUniformLocation(shaderProg, "textureStrength"), waterTextureStrength);
		glUniform1f(glGetUniformLocation(shaderProg, "colorStrength"), waterColorStrength);
		glUniform1f(glGetUniformLocation(shaderProg, "waterAlpha"), waterAlpha);


		glUniformMatrix4fv(modelviewParameter, 1, GL_FALSE, glm::value_ptr(model_view));
		glUniformMatrix4fv(modelParameter, 1, GL_FALSE, glm::value_ptr(model));
		glUniformMatrix3fv(normalMatParameter, 1, GL_FALSE, glm::value_ptr(normal_mat));
		glUniform3f(lightPosParameter, light_pos[0], light_pos[1], light_pos[2]);
		glUniform3f(viewPosParameter, camera_pos[0], camera_pos[1], camera_pos[2]);

		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glDrawElements(GL_TRIANGLES, (water_surface.N - 1) * (water_surface.M - 1) * 2 * 3, GL_UNSIGNED_INT, 0);

		// Draw light source
		
		glUseProgram(shaderProgLight);
		glBindVertexArray(light_vao);

		model = glm::mat4();
		model = glm::translate(model, glm::vec3(light_pos[0], light_pos[1], light_pos[2]));
		model = glm::scale(model, glm::vec3(0.2f));
		model_view = proj * view * model;
		glUniformMatrix4fv(glGetUniformLocation(shaderProgLight, "modelview"), 1, GL_FALSE, glm::value_ptr(model_view));

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
	glDeleteVertexArrays(1, &light_vao);
	glDeleteBuffers(1, &light_vbo);
	glDeleteVertexArrays(1, &floor_vao);
	glDeleteBuffers(1, &floor_vbo);
	glDeleteProgram(shaderProg);
	glDeleteProgram(shaderProgCaustic);
	glDeleteProgram(shaderProgLight);
	glfwDestroyWindow(window);
	glfwTerminate();
	return 0;
}

void mouse_callback(GLFWwindow* window, double xposIn, double yposIn)
{
	float xpos = static_cast<float>(xposIn);
	float ypos = static_cast<float>(yposIn);

	if (is_first_mouse)
	{
		last_x = xpos;
		last_y = ypos;
		is_first_mouse = false;
	}

	float xoffset = xpos - last_x;
	float yoffset = last_y - ypos; // reversed since y-coordinates go from bottom to top
	last_x = xpos;
	last_y = ypos;

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
	camera_front = glm::normalize(front);
}

void processInput(GLFWwindow* window)
{
	if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
		glfwSetWindowShouldClose(window, true);

	float camera_speed = player_speed * delta_time;
	if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
		//cout << "  w" << endl;
		camera_pos += camera_speed * camera_front / 4.0f;
	if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
		camera_pos -= camera_speed * camera_front / 4.0f;
	if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
		camera_pos -= glm::normalize(glm::cross(camera_front, camera_up)) * camera_speed;
	if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
		camera_pos += glm::normalize(glm::cross(camera_front, camera_up)) * camera_speed;
	if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS)
		camera_pos += glm::normalize(glm::normalize(camera_up)) * camera_speed;
	if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS)
		camera_pos -= glm::normalize(glm::normalize(camera_up)) * camera_speed;
	//cameraPos.y = 1.0f;
}