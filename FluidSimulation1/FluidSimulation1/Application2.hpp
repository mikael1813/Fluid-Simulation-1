#pragma once

#include <SDL.h>

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <cmath>

#include "Environment.hpp"

class Application2
{
public:
    Application2();
    ~Application2();

    void events();
    void loop();
    void render();

    void update(float dt);

    void mousePress(SDL_MouseButtonEvent& b);
private:

    GLFWwindow* m_window;

    Environment* m_environment;

    int m_width, m_height;
};