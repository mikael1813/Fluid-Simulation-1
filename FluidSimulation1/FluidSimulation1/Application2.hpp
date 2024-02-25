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

    void DrawCircle(float x, float y, float radius, int num_segments);

    void update(float dt);

    void mousePress(SDL_MouseButtonEvent& b);
private:
    /*SDL_Window* m_window;
    SDL_Surface* m_window_surface;
    SDL_Renderer* m_renderer;*/

    GLFWwindow* m_window;

    Environment* m_environment;

    int m_width, m_height;

    //SDL_Event m_window_event;

    bool m_is_running;
};