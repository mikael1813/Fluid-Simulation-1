#pragma once

#include "Particle.hpp"

#include "SDL.h"

#include <vector>

struct MatrixComponenets {
	std::vector<Particle> particles;
	std::vector<Surface2D> obstacles;
};


class Environment
{
public:
	Environment();
	~Environment();

	void render(SDL_Renderer* renderer);
	void update(float dt);

private:
	std::vector<Particle> m_Particles;
	std::vector<Surface2D> m_Obstacles;

	std::vector<std::vector<MatrixComponenets>> m_InteractionsMatrix;
};