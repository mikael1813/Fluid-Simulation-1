#pragma once

#include "Particle.hpp"

#include <vector>

class Pipe {
public:
	void update(float dt, std::vector<Particle*>& particles, std::vector<Particle*> surroundingParticles, float particleSize);

	Pipe(Vector2D position) : m_Position(position) {
		m_InteractionRadius = 50.0f;
		m_Pressure = 0.0f;
	}

	Vector2D getPosition() {
		return m_Position;
	}

	float getInteractionRadius() {
		return m_InteractionRadius;
	}

private:
	Vector2D m_Position;
	float m_InteractionRadius;
	float m_Pressure;
	int m_ID = 10000;
};