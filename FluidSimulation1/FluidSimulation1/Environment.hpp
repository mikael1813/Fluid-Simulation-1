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

	float calculateDensity(Vector2D point);
	float calculateProperty(Vector2D point);

	void updateInteractionMatrix();
	std::vector<Particle> getParticlesInCell(Vector2D particlePosition);

	void render(int width, int height);
	void update(float dt);

private:
	std::vector<Particle> m_Particles;
	std::vector<float> m_ParticleProperties;
	std::vector<float> m_ParticleDensities;

	Vector2D calculatePressureForce(int particleIndex);
	std::vector<Surface2D> m_Obstacles;

	std::vector<std::vector<MatrixComponenets>> m_InteractionsMatrix;
};