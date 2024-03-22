#pragma once

#include "Particle.hpp"

#include <vector>

class Pipe {
public:
	virtual void update(float dt, std::vector<Particle*>& particles, std::vector<Particle*> surroundingParticles, float particleSize) = 0;

	Vector2D getPosition() {
		return m_Position;
	}

	float getInteractionRadius() {
		return m_InteractionRadius;
	}

protected:
	Vector2D m_Position;
	float m_InteractionRadius;
	float m_Pressure;
	int m_ID = 10000;
};

class GeneratorPipe : public Pipe {

public:

	GeneratorPipe(Vector2D position) {
		m_Position = position;
		m_InteractionRadius = 50.0f;
		m_Pressure = 0.0f;
	}
	void update(float dt, std::vector<Particle*>& particles, std::vector<Particle*> surroundingParticles, float particleSize);
};


class ConsumerPipe : public Pipe {

public:

	ConsumerPipe(Vector2D position) {
		m_Position = position;
		m_InteractionRadius = 50.0f;
		m_Pressure = 0.0f;
	}

	void update(float dt, std::vector<Particle*>& particles, std::vector<Particle*> surroundingParticles, float particleSize);
};