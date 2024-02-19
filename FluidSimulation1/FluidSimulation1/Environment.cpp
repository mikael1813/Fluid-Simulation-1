#include "Environment.hpp"

#include <SDL.h>

constexpr auto particleCount = 10;
constexpr auto particleRadius = 2;
constexpr auto particleDistance = 5;


Environment::Environment() {
	m_particles = std::vector<Particle>{};

	for (int i = 0; i < particleCount; i++) {
		for (int j = 0; j < particleCount; j++) {
			m_particles.push_back(Particle(200+i * particleDistance,200+ j * particleDistance));
		}
	}
}

Environment::~Environment()
{
}

void DrawCircle(SDL_Renderer* renderer, int centreX, int centreY, int radius)
{
    const int32_t diameter = (radius * 2);

    int32_t x = (radius - 1);
    int32_t y = 0;
    int32_t tx = 1;
    int32_t ty = 1;
    int32_t error = (tx - diameter);

    while (x >= y)
    {
        //  Each of the following renders an octant of the circle
        SDL_RenderDrawPoint(renderer, centreX + x, centreY - y);
        SDL_RenderDrawPoint(renderer, centreX + x, centreY + y);
        SDL_RenderDrawPoint(renderer, centreX - x, centreY - y);
        SDL_RenderDrawPoint(renderer, centreX - x, centreY + y);
        SDL_RenderDrawPoint(renderer, centreX + y, centreY - x);
        SDL_RenderDrawPoint(renderer, centreX + y, centreY + x);
        SDL_RenderDrawPoint(renderer, centreX - y, centreY - x);
        SDL_RenderDrawPoint(renderer, centreX - y, centreY + x);

        if (error <= 0)
        {
            ++y;
            error += ty;
            ty += 2;
        }

        if (error > 0)
        {
            --x;
            tx += 2;
            error += (tx - diameter);
        }
    }
}

void Environment::draw(SDL_Renderer* renderer)
{
    for (auto& particle : m_particles) {
        DrawCircle(renderer, particle.m_position_x, particle.m_position_y, particleRadius);
    }
}


