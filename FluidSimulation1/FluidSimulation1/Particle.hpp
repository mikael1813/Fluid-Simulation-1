#pragma once

#include "Phisics.hpp"

#include <iostream>

constexpr float gravity = 99.8f;


class Particle {
public:
	Particle(float x, float y) : m_Position(Point2D(x, y)) {}
	Point2D m_Position;
	Point2D m_LastSafePosition;
	Vector2D m_Velocity;

	void update(float dt) {
		m_LastSafePosition = m_Position;
		if (dt == 0) {
			return;
		}

		m_Velocity.Y += gravity * dt;

		m_Position.X += m_Velocity.X * dt;
		m_Position.Y += m_Velocity.Y * dt;
	}

private:
	
	float visible_radius = 2.0f;
	float mass = 1.0f;

};