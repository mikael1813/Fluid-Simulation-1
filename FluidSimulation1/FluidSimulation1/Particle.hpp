#pragma once

#include "Phisics.hpp"

#include <iostream>

constexpr float GRAVITY = 9.0f;


class Particle {
public:
	Particle(float x, float y, int id) : m_Position(Vector2D(x, y)), m_ID(id) {}
	Vector2D m_Position;
	Vector2D m_LastSafePosition;
	Vector2D m_Velocity;
	float m_Density = 0.0f;
	int m_ID;

	void update(float dt) {
		if (dt == 0) {
			return;
		}

		m_LastSafePosition = m_Position;

		Vector2D gravity(0.0f, GRAVITY);

		//m_Velocity += gravity * dt;

		m_Position.X += m_Velocity.X * dt;
		m_Position.Y += m_Velocity.Y * dt;

		m_Velocity = m_Velocity * 0.99f;
	}

private:

	float visible_radius = 2.0f;
	float mass = 1.0f;
};