#include <iostream>
#include <SDL.h>
#include "Application2.hpp"
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
//#undef main
using namespace std;

int main(int argc, char* args[]) {

    Application2 app;

    app.loop();
    app.render();

    _CrtSetReportMode(_CRT_WARN, _CRTDBG_MODE_DEBUG);
    _CrtDumpMemoryLeaks();

    return 0;
}

//#include <GL/glew.h>
//#include <GLFW/glfw3.h>
//#include <iostream>
//#include <cmath>
//
//void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
//	glViewport(0, 0, width, height);
//}
//
//void DrawCircle(float x, float y, float radius, int num_segments) {
//	glBegin(GL_TRIANGLE_FAN);
//	for (int ii = 0; ii < num_segments; ii++) {
//		float theta = 2.0f * 3.1415 * float(ii) / float(num_segments);
//		float dx = radius * cos(theta);
//		float dy = radius * sin(theta);
//		glVertex2f(x + dx, y + dy);
//	}
//	glEnd();
//}
//
//int main() {
//	// Initialize GLFW
//	if (!glfwInit()) {
//		std::cerr << "Failed to initialize GLFW" << std::endl;
//		return -1;
//	}
//
//	// Create a windowed mode window and its OpenGL context
//	GLFWwindow* window = glfwCreateWindow(1280, 720, "Circle Example", NULL, NULL);
//	if (!window) {
//		std::cerr << "Failed to create GLFW window" << std::endl;
//		glfwTerminate();
//		return -1;
//	}
//
//	// Set framebuffer size callback
//	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
//
//	// Make the window's context current
//	glfwMakeContextCurrent(window);
//
//	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
//	glEnable(GL_BLEND);
//
//	// Loop until the user closes the window
//	while (!glfwWindowShouldClose(window)) {
//		// Render here
//		glClear(GL_COLOR_BUFFER_BIT);
//
//		// Get the size of the window
//		int width, height;
//		glfwGetFramebufferSize(window, &width, &height);
//
//		// Set up the viewport to maintain aspect ratio
//		glViewport(0, 0, width, height);
//
//		// Calculate the aspect ratio
//		float aspectRatio = static_cast<float>(width) / static_cast<float>(height);
//
//		// Apply aspect ratio to your projection matrix or adjust accordingly
//		// For example:
//		glMatrixMode(GL_PROJECTION);
//		glLoadIdentity();
//		glOrtho(-aspectRatio, aspectRatio, -1.0f, 1.0f, -1.0f, 1.0f);
//
//		float x = 640.0f, y = 360.0f, radius = 100.0f;
//		float ndcX = (2.0f * x) / width - 1.0f;
//		float ndcX2 = (2.0f * (x+100.0f)) / width - 1.0f;
//		float ndcY = 1.0f - (2.0f * y) / height;
//		float ndcRadius = (2.0f * radius) / width;
//
//		// Draw a circle
//		
//		glColor4f(0.0, 0.5, 0.5, 0.1);
//		DrawCircle(ndcX, ndcY, ndcRadius, 50);
//		DrawCircle(ndcX2, ndcY, ndcRadius, 50);
//
//		// Swap front and back buffers
//		glfwSwapBuffers(window);
//
//		// Poll for and process events
//		glfwPollEvents();
//	}
//
//	glfwTerminate();
//	return 0;
//}
