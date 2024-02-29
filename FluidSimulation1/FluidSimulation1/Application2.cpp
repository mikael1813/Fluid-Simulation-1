#include "Application2.hpp"
#include "Timer.hpp"



void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
	glViewport(0, 0, width, height);
}

Application2::Application2()
{

	m_is_running = true;


	// Initialize GLFW
	if (!glfwInit()) {
		std::cerr << "Failed to initialize GLFW" << std::endl;
		return;
	}

	// Create a windowed mode window and its OpenGL context
	m_window = glfwCreateWindow(1280, 720, "Circle Example", NULL, NULL);
	if (!m_window) {
		std::cerr << "Failed to create GLFW window" << std::endl;
		glfwTerminate();
		return;
	}

	// Set framebuffer size callback
	glfwSetFramebufferSizeCallback(m_window, framebuffer_size_callback);

	// Make the window's context current
	glfwMakeContextCurrent(m_window);

	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_BLEND);

	m_environment = new Environment();
}

Application2::~Application2()
{
	/*SDL_FreeSurface(m_window_surface);
	SDL_DestroyWindow(m_window);*/
	glfwTerminate();
	delete m_environment;
}

void Application2::events() {

	// Poll for and process events
	glfwPollEvents();
}

void Application2::loop()
{
	float time_passed = 0.0f;
	int frames = 0;
	std::chrono::steady_clock::time_point lastTime = std::chrono::steady_clock::now();
	// Loop until the user closes the window
	while (!glfwWindowShouldClose(m_window)) {

		double deltaTime = Timer::GetInstance()->GetTime();


		//std::cout << deltaTime << std::endl;

		//std::chrono::steady_clock::time_point time1 = std::chrono::steady_clock::now();

		this->update(deltaTime);

		//std::chrono::steady_clock::time_point time2 = std::chrono::steady_clock::now();
		//double tick = std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count();

		//time1 = std::chrono::steady_clock::now();

		// Render here
		this->render();

		//time2 = std::chrono::steady_clock::now();
		//tick = std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count();

		this->events();

		std::chrono::steady_clock::time_point time = std::chrono::steady_clock::now();
		double tick = std::chrono::duration_cast<std::chrono::microseconds>(time - lastTime).count() / 1000000.0f;
		lastTime = time;

		time_passed += tick;
		frames++;

		if (time_passed >= 1.0f) {
			std::cout << "FPS: " << frames << " " << std::endl;
			time_passed = 0.0f;
			frames = 0;
		}
	}
}

void Application2::mousePress(SDL_MouseButtonEvent& b) {
	if (b.button == SDL_BUTTON_LEFT) {

	}
	if (b.button)
		if (b.button == SDL_BUTTON_RIGHT) {
			int x = b.x;
			int y = b.y;
		}
}

void Application2::render()
{
	// Render here
	glClear(GL_COLOR_BUFFER_BIT);

	// Get the size of the window
	glfwGetFramebufferSize(m_window, &m_width, &m_height);

	// Set up the viewport to maintain aspect ratio
	glViewport(0, 0, m_width, m_height);

	// Calculate the aspect ratio
	float aspectRatio = static_cast<float>(m_width) / static_cast<float>(m_height);

	// Apply aspect ratio to your projection matrix or adjust accordingly
	// For example:
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-aspectRatio, aspectRatio, -1.0f, 1.0f, -1.0f, 1.0f);

	/*float x = 640.0f, y = 360.0f, radius = 100.0f;
	float ndcX = (2.0f * x) / m_width - 1.0f;
	float ndcX2 = (2.0f * (x + 100.0f)) / m_width - 1.0f;
	float ndcY = 1.0f - (2.0f * y) / m_height;
	float ndcRadius = (2.0f * radius) / m_width;*/

	// Draw a circle

	m_environment->render(m_width, m_height);

	// Swap front and back buffers
	glfwSwapBuffers(m_window);
}

void Application2::update(float dt)
{
	m_environment->update(dt);
}
