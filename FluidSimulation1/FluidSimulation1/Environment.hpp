#pragma once

#include "Particle.hpp"

#include "SDL.h"

#include <vector>


class Environment
{
public:
	Environment();
	~Environment();

	void render(SDL_Renderer* renderer);
	void update(float dt);

private:
	std::vector<Particle> m_particles;
	std::vector<Surface2D> m_obstacles;
};