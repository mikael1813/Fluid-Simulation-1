#include "Pipe.hpp"
#include <random>
#include <list>
#include "Phisics.hpp"

void Pipe::update(float dt, std::vector<Particle*>& particles, std::vector<Particle*> surroundingParticles, float particleSize) {

	std::list<Particle*> particlesToAdd;

	// Seed the random number generator
	std::random_device rd;
	std::mt19937 gen(rd());

	for (int i = 0; i < constants::m_PI * m_InteractionRadius * m_InteractionRadius / (particleSize * 100); i++) {
		int posX = std::uniform_int_distribution<int>(m_Position.X - m_InteractionRadius, m_Position.X + m_InteractionRadius)(gen);
		int posY = std::uniform_int_distribution<int>(m_Position.Y - m_InteractionRadius, m_Position.Y + m_InteractionRadius)(gen);

		particlesToAdd.push_back(new Particle(posX, posY, m_ID++));
	}

	for (auto& particle : surroundingParticles) {
		Vector2D direction = particle->getPosition() - m_Position;

		if (direction.getMagnitude() == 0) {
			direction = Vector2D::getRandomDirection();
		}

		float distance = direction.getMagnitude();
		if (distance <= m_InteractionRadius) {
			particle->addForce(direction / distance * m_Pressure);
		}

		std::list<Particle*>::iterator otherParticle = particlesToAdd.begin();

		std::cout << "A" << std::endl;

		// here lies a error down below
		// to do: solve it

		while (otherParticle != particlesToAdd.end()) {
			float distance2 = (particle->getPosition() - (*otherParticle)->getPosition()).getMagnitude();

			if (distance2 <= m_InteractionRadius) {
				delete* otherParticle;
				particlesToAdd.erase(otherParticle);
			}

			otherParticle++;
		}

		std::cout << "b" << std::endl;;
	}

	for (auto& particle : particlesToAdd) {
		particles.push_back(particle);
	}
}