#include "Environment.hpp"
#include "InteractionMatrixClass.hpp"

#include <algorithm>

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <cmath>
#include <omp.h>
#include <random>

//#include <SDL.h>
#include <chrono>

#include <mutex> 

constexpr auto particleCount = 2000;
constexpr auto particleRadius = 2;
constexpr auto particleRadiusOfRepel = 50;
constexpr auto particleDistance = 30;

constexpr auto particleRepulsionForce = 5.0f;

constexpr int SCREEN_WIDTH = 1280;
constexpr int SCREEN_HEIGHT = 720;

constexpr float viscosityStrength = 0.1f;

constexpr int THREAD_COUNT = 8;

float ExampleFunction(Vector2D point) {
	return cos(point.Y - 3 + sin(point.X));
}


Environment::Environment() {


	InteractionMatrixClass::getInstance()->initializeMatrix(SCREEN_WIDTH, SCREEN_HEIGHT, particleRadiusOfRepel);

	m_Particles = std::vector<Particle*>{};

	// Seed the random number generator
	std::random_device rd;
	std::mt19937 gen(rd());

	int count = 0;

	for (int i = 0; i < particleCount; i++) {

		bool ok;
		float posX;
		float posY;
		do {
			ok = true;

			posX = std::uniform_int_distribution<int>(100, SCREEN_WIDTH - 100)(gen);
			posY = std::uniform_int_distribution<int>(100, SCREEN_HEIGHT - 100)(gen);

			for (auto& particle : m_Particles) {
				if (particle->m_Position.X == posX && particle->m_Position.Y == posY) {
					ok = false;
					break;
				}
			}
		} while (!ok);

		m_Particles.push_back(new Particle(posX, posY, count++));
		m_ParticleProperties.push_back(ExampleFunction(Vector2D(posX, posY)));
		m_ParticleDensities.push_back(0.0f);

	}

	m_Obstacles.push_back(Surface2D(50, 10, 1200, 11));
	m_Obstacles.push_back(Surface2D(50, 10, 50, 700));
	m_Obstacles.push_back(Surface2D(50, 699, 1200, 700));
	m_Obstacles.push_back(Surface2D(1200, 10, 1200, 700));

	/*m_Obstacles.push_back(Surface2D(500, 400, 600, 300));
	m_Obstacles.push_back(Surface2D(600, 300, 700, 400));
	m_Obstacles.push_back(Surface2D(700, 400, 500, 400));*/


}

Environment::~Environment()
{
	for (auto& particle : m_Particles) {
		delete particle;
	}
}

//void DrawCircle(SDL_Renderer* renderer, int centreX, int centreY, int radius)
//{
//	const int32_t diameter = (radius * 2);
//
//	int32_t x = (radius - 1);
//	int32_t y = 0;
//	int32_t tx = 1;
//	int32_t ty = 1;
//	int32_t error = (tx - diameter);
//
//	while (x >= y)
//	{
//		//  Each of the following renders an octant of the circle
//		SDL_RenderDrawPoint(renderer, centreX + x, centreY - y);
//		SDL_RenderDrawPoint(renderer, centreX + x, centreY + y);
//		SDL_RenderDrawPoint(renderer, centreX - x, centreY - y);
//		SDL_RenderDrawPoint(renderer, centreX - x, centreY + y);
//		SDL_RenderDrawPoint(renderer, centreX + y, centreY - x);
//		SDL_RenderDrawPoint(renderer, centreX + y, centreY + x);
//		SDL_RenderDrawPoint(renderer, centreX - y, centreY - x);
//		SDL_RenderDrawPoint(renderer, centreX - y, centreY + x);
//
//		if (error <= 0)
//		{
//			++y;
//			error += ty;
//			ty += 2;
//		}
//
//		if (error > 0)
//		{
//			--x;
//			tx += 2;
//			error += (tx - diameter);
//		}
//	}
//}

void DrawCircle(int width, int height, float x, float y, float radius, int num_segments) {

	float raport = (float)width / (float)height;

	x = (2.0f * x) / width - 1.0f;
	x = x * raport;

	y = 1.0f - (2.0f * y) / height;

	radius = (2.0f * radius) / width;

	glBegin(GL_TRIANGLE_FAN);
	for (int ii = 0; ii < num_segments; ii++) {
		float theta = 2.0f * 3.1415 * float(ii) / float(num_segments);
		float dx = radius * cos(theta);
		float dy = radius * sin(theta);
		glVertex2f(x + dx, y + dy);
	}
	glEnd();
}

void DrawLine(int width, int height, Vector2D a, Vector2D b) {

	float raport = (float)width / (float)height;

	a.X = (2.0f * a.X) / width - 1.0f;
	a.X = a.X * raport;

	a.Y = 1.0f - (2.0f * a.Y) / height;

	b.X = (2.0f * b.X) / width - 1.0f;
	b.X = b.X * raport;

	b.Y = 1.0f - (2.0f * b.Y) / height;

	glBegin(GL_LINES);
	glVertex2f(a.X, a.Y);
	glVertex2f(b.X, b.Y);
	glEnd();
}

//void drawSurface2D(SDL_Renderer* renderer, Surface2D surface) {
//	SDL_RenderDrawLine(renderer, surface.Point1.X, surface.Point1.Y, surface.Point2.X, surface.Point2.Y);
//}

// Function to map velocity to color (blue for velocity 0, red for max velocity)
void velocityToColor(float velocity, float& red, float& green, float& blue) {
	float maxVelocity = 50.0f; // Maximum velocity in the simulation

	// Normalize velocity between 0 and 1
	float normalizedVelocity = std::min(std::max(velocity / maxVelocity, 0.0f), 1.0f);

	// Linear interpolation between blue and red based on normalized velocity
	red = normalizedVelocity;      // Red component increases with velocity
	green = 0.5f - normalizedVelocity / 2;                  // No green component
	blue = 1.0f - normalizedVelocity; // Blue component decreases with velocity
}

void Environment::renderParticles(int width, int height) {
	for (auto& particle : m_Particles) {

		//float density = particle->m_Density;

		Vector2D vc = particle->m_Velocity;

		//float color = density / maxDensity;

		float blue, green, red;
		velocityToColor(particle->m_Velocity.magnitude(), red, green, blue);

		glColor4f(red, green, blue, 1.0f);
		DrawCircle(width, height, particle->m_Position.X, particle->m_Position.Y, particleRadius * 2, 20);

		/*glColor4f(1.0, 1.0, 1.0, 0.4f);
		DrawLine(width, height, particle->m_Position, particle->m_Position + vc);*/
	}
}

void Environment::render(int width, int height)
{
	/*float maxDensity = 0.0f;
	for (auto& particle : m_Particles) {

		if (particle->m_Density > maxDensity) {
			maxDensity = particle->m_Density;
		}
	}*/

	this->renderParticles(width, height);

	for (auto& obstacle : m_Obstacles) {
		glColor3f(1.0, 1.0, 1.0);
		DrawLine(width, height, obstacle.Point1, obstacle.Point2);
	}
}

//float smoothing_kernel(float radius, float distance) {
//
//	//float volume = 3.1415f * pow(radius, 2);
//
//	float value = std::max(0.0f, (radius - distance) / radius);
//
//	return value * value * value;
//}

float smoothingKernel(float radius, float distance) {
	if (distance >= radius) {
		return 0.0f;
	}

	float x = (radius - distance) / radius;
	return x * x;
}

float smoothingKernelDerivative(float radius, float distance) {
	if (distance >= radius) {
		return 0.0f;
	}
	float x = (radius - distance) / radius;
	return 2 * x;
}

float viscositySmoothingKernel(float radius, float distance) {
	if (distance >= radius) {
		return 0.0f;
	}
	float x = (radius * radius - distance * distance) / (radius * radius);
	return x * x * x;
}

float Environment::calculateDensity(Vector2D point) {

	constexpr auto scalar = 1000;

	float density = 0.0f;
	const float mass = 1.0f;

	for (auto& particle : InteractionMatrixClass::getInstance()->getParticlesInCell(point, particleRadiusOfRepel)) {
		float distance = sqrt(Math::squared_distance(point, particle->m_Position));
		float influence = smoothingKernel(particleRadiusOfRepel, distance);
		density += mass * influence;
	}

	float volume = 3.1415f * pow(particleRadiusOfRepel, 2);

	return density / volume * scalar;

}

float Environment::calculateProperty(Vector2D point) {
	float property = 0.0f;
	const float mass = 1.0f;

	for (int i = 0; i < m_Particles.size(); i++) {
		float distance = sqrt(Math::squared_distance(point, m_Particles.at(i)->m_Position));

		float influence = smoothingKernel(particleRadiusOfRepel, distance);

		float density = m_ParticleDensities.at(i);

		property += m_ParticleProperties.at(i) * influence * mass / density;
	}

	return property;
}

//Vector2D Environment::calculatePropertyGradient(Vector2D point) {
//	Vector2D propertyGradient = Vector2D();
//	const float mass = 1.0f;
//
//	for (int i = 0; i < m_Particles.size(); i++) {
//		float distance = sqrt(squared_distance(point, m_Particles.at(i).m_Position));
//		if (distance == 0) {
//			continue;
//		}
//		Vector2D dir = (m_Particles.at(i).m_Position - point) / -distance;
//		float slope = smoothingKernelDerivative(particleRadiusOfRepel, distance);
//
//		float density = m_ParticleDensities.at(i);
//
//		propertyGradient += dir * slope * mass / density;
//	}
//
//	return propertyGradient;
//}

Vector2D Environment::calculateViscosityForce(Particle* particle) {

	Vector2D viscosityForce = Vector2D();
	Vector2D position = particle->m_Position;

	for (auto& otherParticle : InteractionMatrixClass::getInstance()->getParticlesInCell(particle->m_Position, particleRadiusOfRepel)) {
		float distance = sqrt(Math::squared_distance(position, otherParticle->m_Position));
		float influence = viscositySmoothingKernel(particleRadiusOfRepel, distance);

		viscosityForce += (otherParticle->m_Velocity - particle->m_Velocity) * influence;
	}

	return viscosityForce * viscosityStrength;
}

float convertDensityToPressure(float density) {
	const float targetDensity = 0.5f;
	//const float pressureConstant = 10.0f;
	const float pressureConstant = 30.0f;

	float densityError = density - targetDensity;
	float pressure = pressureConstant * densityError;
	return pressure;
}

float calculateSharedPressure(float density1, float density2) {
	float pressure1 = convertDensityToPressure(density1);
	float pressure2 = convertDensityToPressure(density2);
	return (pressure1 + pressure2) / 2;
}

Vector2D Environment::calculatePressureForce(Particle* particle) {
	Vector2D pressureForce = Vector2D();
	const float mass = 1.0f;

	for (auto& otherParticle : InteractionMatrixClass::getInstance()->getParticlesInCell(particle->m_Position, particleRadiusOfRepel)) {
		//for (auto& otherParticle :m_Particles) {

		if (particle->m_ID == otherParticle->m_ID) {
			continue;
		}

		float distance = sqrt(Math::squared_distance(particle->m_Position, otherParticle->m_Position));
		if (distance < particleRadius) {
			int tt = 0;
		}
		Vector2D dir = distance < particleRadius ? Vector2D::getRandomDirection() : (otherParticle->m_Position - particle->m_Position) / distance;
		float slope = smoothingKernelDerivative(particleRadiusOfRepel, distance);


		float density = otherParticle->m_Density;

		float sharedPressure = calculateSharedPressure(density, otherParticle->m_Density);

		pressureForce += -sharedPressure * dir * slope * mass / density;
	}

	return pressureForce;
}

void Environment::updateParticleDensities(int start, int end) {
	for (int i = start; i < end; i++) {

		Particle* particle = m_Particles.at(i);

		/*Vector2D futurePosition = particle->m_Position + particle->m_Velocity * dt;
		particle->m_Density = calculateDensity(futurePosition);*/

		particle->m_Density = calculateDensity(particle->m_Position);

	}
}

void Environment::updateParticles(double dt, int start, int end) {
	for (int i = start; i < end; i++) {

		Particle* particle = m_Particles.at(i);

		Vector2D pressureForce = calculatePressureForce(particle);
		Vector2D pressureAcceleration = pressureForce / particle->m_Density;

		particle->m_Velocity += pressureAcceleration * dt;

		Vector2D viscosityForce = calculateViscosityForce(particle);

		particle->m_Velocity += viscosityForce * dt;

		particle->update(dt);


		for (auto& obstacle : m_Obstacles) {
			if (Math::check_line_segment_circle_intersection(obstacle.Point1, obstacle.Point2, particle->m_Position, particleRadius)) {
				/*Vector2D normalVector = Math::calculateNormalVector(Math::calculateSlope(obstacle.Point1, obstacle.Point2));
				Vector2D reflectionVector = Math::calculateReflectionVector(particle->m_Velocity, normalVector);*/

				// magnitude of reflection vector
				//float magnitude = sqrt(reflectionVector.X * reflectionVector.X + reflectionVector.Y * reflectionVector.Y);

				// normalize the reflection vector
				/*reflectionVector.X /= magnitude;
				reflectionVector.Y /= magnitude;*/

				//particle->m_Velocity = reflectionVector * 0.1f;
				particle->m_Velocity = Vector2D();

				/*particle.m_Velocity.Y = -particle.m_Velocity.Y;*/

				particle->m_Position = particle->m_LastSafePosition;

				break;
			}
		}

		//for (auto& otherParticle : m_Particles) {
		for (auto& otherParticle : InteractionMatrixClass::getInstance()->getParticlesInCell(particle->m_Position, particleRadiusOfRepel)) {
			if (particle->m_ID == otherParticle->m_ID) {
				continue;
			}
			//if (squared_distance(m_Particles.at(i).m_Position, m_Particles.at(j).m_Position) <= (particleRadius * particleRadius * 4)) {

				//Vector2D normalVector = Vector2D(m_Particles.at(j).m_Position.X - m_Particles.at(i).m_Position.X, m_Particles.at(j).m_Position.Y - m_Particles.at(i).m_Position.Y);

				////magnitude of normal vector
				//float magnitude = sqrt(normalVector.X * normalVector.X + normalVector.Y * normalVector.Y);

				//// normalize the normal vector
				//normalVector.X /= magnitude;
				//normalVector.Y /= magnitude;

				//Vector2D reflectionVector = Math::calculateReflectionVector(m_Particles.at(i).m_Velocity, normalVector);

				////m_Particles.at(i).m_Velocity = reflectionVector * 0.1f;


				//m_Particles.at(i).m_Position = m_Particles.at(i).m_LastSafePosition;
			//}
			if (Math::squared_distance(particle->m_Position, otherParticle->m_Position) <= (particleRadius * particleRadius) * 4) {
				Vector2D normalVector = Vector2D(otherParticle->m_Position.X - particle->m_Position.X, otherParticle->m_Position.Y - particle->m_Position.Y);

				//magnitude of normal vector
				float magnitude = -1 * sqrt(normalVector.X * normalVector.X + normalVector.Y * normalVector.Y);

				// normalize the normal vector
				normalVector.X /= magnitude;
				normalVector.Y /= magnitude;

				//float power = smoothing_kernel(particleRadiusOfRepel - particleRadius, sqrt(squared_distance(m_Particles.at(i).m_Position, m_Particles.at(j).m_Position)) - particleRadius * 2);
				float power = 1;
				power *= particleRepulsionForce;

				particle->m_Velocity += normalVector * power;

				otherParticle->m_Velocity -= normalVector * power;

				int x = 0;

				//particle.m_Position = particle.m_LastSafePosition;
			}
		}
	}
}

void Environment::parallelUpdateParticleDensities() {
	for (int i = 0; i < THREAD_COUNT; i++) {
		m_Threads.push_back(std::thread(&Environment::updateParticleDensities, this, i * m_Particles.size() / THREAD_COUNT, (i + 1) * m_Particles.size() / THREAD_COUNT));
	}

	for (int i = 0; i < THREAD_COUNT; i++) {
		m_Threads.at(i).join();
	}

	m_Threads.clear();
}
void Environment::parallelUpdateParticles(double dt) {
	for (int i = 0; i < THREAD_COUNT; i++) {
		m_Threads.push_back(std::thread(&Environment::updateParticles, this, dt, i * m_Particles.size() / THREAD_COUNT, (i + 1) * m_Particles.size() / THREAD_COUNT));
	}

	for (int i = 0; i < THREAD_COUNT; i++) {
		m_Threads.at(i).join();
	}

	m_Threads.clear();
}

void Environment::update(float dt) {

	std::chrono::steady_clock::time_point time1 = std::chrono::steady_clock::now();
	InteractionMatrixClass::getInstance()->updateInteractionMatrix(m_Particles, particleRadiusOfRepel);

	std::chrono::steady_clock::time_point time2 = std::chrono::steady_clock::now();
	double tick = std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count();

	time1 = std::chrono::steady_clock::now();

	this->parallelUpdateParticleDensities();

	time2 = std::chrono::steady_clock::now();
	tick = std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count();

	/*time1 = std::chrono::steady_clock::now();

	this->updateParticleDensities(0, m_Particles.size());

	time2 = std::chrono::steady_clock::now();
	tick = std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count();*/

	time1 = std::chrono::steady_clock::now();

	this->parallelUpdateParticles(dt);

	time2 = std::chrono::steady_clock::now();
	tick = std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count();

	/*time1 = std::chrono::steady_clock::now();

	this->updateParticles(dt, 0, m_Particles.size());

	time2 = std::chrono::steady_clock::now();
	tick = std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count();*/

	int aba = 0;
}
