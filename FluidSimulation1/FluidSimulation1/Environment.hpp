#pragma once

#include "Particle.hpp"

#include "SDL.h"

#include <vector>
#include "thread"

struct MatrixComponenets {
	std::vector<Particle*> particles;
	std::vector<Surface2D*> obstacles;
};


class Environment
{
public:
	Environment();
	~Environment();


	void render(int width, int height);
	void update(float dt);

private:
	std::vector<Particle*> m_Particles;
	std::vector<float> m_ParticleProperties;
	std::vector<float> m_ParticleDensities;

	std::vector<std::thread> m_Threads;

	float calculateDensity(Vector2D point);
	float calculateProperty(Vector2D point);

	void renderParticles(int width, int height);

	void updateParticleDensities(int start, int end);
	void calculateFutureVelocitiesAndCheckCollisions(double dt, int start, int end);

	void parallelUpdateParticleDensities();
	void parallelUpdateParticles(double dt);

	Vector2D calculateViscosityForce(Particle* particle);

	Vector2D calculatePressureForce(Particle* particle);
	void addToInteractionMatrixCellSurroundingCells(int x, int y, std::vector<std::vector<MatrixComponenets>>& temporary);
	void updateInteractionMatrix();
	std::vector<Particle*> getParticlesInCell(Vector2D particlePosition);


	std::vector<Surface2D> m_Obstacles;

	std::vector<std::vector<MatrixComponenets>> m_InteractionsMatrix;
};