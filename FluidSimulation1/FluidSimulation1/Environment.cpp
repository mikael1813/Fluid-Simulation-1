#include "Environment.hpp"

#include <algorithm>

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <cmath>
#include <omp.h>
#include <random>

#include <SDL.h>

constexpr auto particleCount = 10;
constexpr auto particleRadius = 5;
constexpr auto particleRadiusOfRepel = 100;
constexpr auto particleDistance = 30;

constexpr auto maximumSpeed = 400.0f;

constexpr auto particleRepulsionForce = 5.0f;

//constexpr auto SCREEN_WIDTH = 1280;
constexpr auto SCREEN_WIDTH = 2300;

constexpr auto MIN_WIDTH = -500;
constexpr auto MAX_WIDTH = 1800;

constexpr auto SCREEN_HEIGHT = 720;

float ExampleFunction(Vector2D point) {
	return cos(point.Y - 3 + sin(point.X));
}


Environment::Environment() {
	for (int i = 0; i < SCREEN_HEIGHT / particleRadiusOfRepel; i++) {
		std::vector<MatrixComponenets> row;
		for (int j = 0; j < SCREEN_WIDTH / particleRadiusOfRepel; j++) {
			row.push_back(MatrixComponenets());
		}
		m_InteractionsMatrix.push_back(row);
	}

	m_Particles = std::vector<Particle>{};

	// Seed the random number generator
	std::random_device rd;
	std::mt19937 gen(rd());

	for (int i = 0; i < particleCount * 2; i++) {
		for (int j = 0; j < particleCount; j++) {
			/*float posX = 200 + i * particleDistance;
			float posY = 200 + j * particleDistance;*/
			float posX = std::uniform_int_distribution<int>(50, 1200)(gen);
			float posY = std::uniform_int_distribution<int>(10, 600)(gen);
			m_Particles.push_back(Particle(posX, posY));
			m_ParticleProperties.push_back(ExampleFunction(Vector2D(posX, posY)));
			m_ParticleDensities.push_back(0.0f);
		}
	}

	m_Obstacles.push_back(Surface2D(50, 10, 1200, 11));
	m_Obstacles.push_back(Surface2D(50, 10, 50, 700));
	m_Obstacles.push_back(Surface2D(50, 699, 1200, 700));
	m_Obstacles.push_back(Surface2D(1200, 10, 1200, 700));


}

Environment::~Environment()
{
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

// Function to calculate the color based on speed
void GetSpeedColor(float speed, float& red, float& green, float& blue) {
	// Normalize the speed between 0 and 1
	float normalizedSpeed = speed / maximumSpeed; // Adjust maximumSpeed as needed

	// Interpolate between cyan and bright red based on the normalized speed
	red = normalizedSpeed;  // Faster objects appear more red
	green = 1.0f - normalizedSpeed;        // Slower objects appear more green
	blue = 1.0f;                    // Set blue component to maximum (cyan)

	// Clamp color values to [0, 1]
	red = std::max(0.0f, std::min(1.0f, red));
	green = std::max(0.0f, std::min(1.0f, green));
	blue = std::max(0.0f, std::min(1.0f, blue));
}

void Environment::render(int width, int height)
{
	float maxDensity = 0.0f;
	for (auto& particleDensity : m_ParticleDensities) {

		if (particleDensity > maxDensity) {
			maxDensity = particleDensity;
		}
	}
	for (auto& particle : m_Particles) {

		float density = calculateDensity(particle.m_Position);

		Vector2D vc = particle.m_Velocity;

		float color = density / maxDensity;

		glColor4f(color, color, color, 1.0f);
		DrawCircle(width, height, particle.m_Position.X, particle.m_Position.Y, particleRadius * 2, 20);

		glColor4f(1.0, 1.0, 1.0, 0.4f);
		DrawLine(width, height, particle.m_Position, particle.m_Position + vc);
		//DrawCircle(renderer, particle.m_Position.X, particle.m_Position.Y, particleRadius);

		//float red, green, blue;
		//// magnitude of velocity
		//float magnitude = sqrt(particle.m_Velocity.X * particle.m_Velocity.X + particle.m_Velocity.Y * particle.m_Velocity.Y);

		//GetSpeedColor(magnitude, red, green, blue);

		//glColor4f(red, green, blue, 0.1);
		//DrawCircle(width, height, particle.m_Position.X, particle.m_Position.Y, particleRadiusOfRepel, 20);

		//DrawCircle(renderer, particle.m_Position.X, particle.m_Position.Y, particleRadiusOfRepel);
	}

	for (auto& obstacle : m_Obstacles) {
		glColor3f(1.0, 1.0, 1.0);
		DrawLine(width, height, obstacle.Point1, obstacle.Point2);
	}

	/*SDL_RenderDrawLine(renderer, 10, 10, 10, 700);

	SDL_RenderDrawLine(renderer, 10, 700, 1200, 700);

	SDL_RenderDrawLine(renderer, 1200, 10, 1200, 700);*/
}

// Function to compute the squared distance between two points
float squared_distance(const Vector2D& p1, const Vector2D& p2) {
	float dx = p1.X - p2.X;
	float dy = p1.Y - p2.Y;
	return dx * dx + dy * dy;
}

// Check if the line segment AB intersects the circle with center C and radius R
bool check_line_segment_circle_intersection(const Vector2D& A, const Vector2D& B, const Vector2D& C, float radius) {
	// Compute squared distances
	float dist_AB_squared = squared_distance(A, B);
	float dist_AC_squared = squared_distance(A, C);
	float dist_BC_squared = squared_distance(B, C);

	// Check if any of the endpoints (A or B) is inside the circle
	if (dist_AC_squared <= radius * radius || dist_BC_squared <= radius * radius) {
		return true;
	}

	// Check if the line segment intersects the circle
	float dot_product = (C.X - A.X) * (B.X - A.X) + (C.Y - A.Y) * (B.Y - A.Y);
	if (dot_product < 0 || dot_product > dist_AB_squared) {
		return false; // Closest point to C is outside the line segment AB
	}

	// Compute the closest point P on the line segment AB to C
	Vector2D P;
	P.X = A.X + (B.X - A.X) * dot_product / dist_AB_squared;
	P.Y = A.Y + (B.Y - A.Y) * dot_product / dist_AB_squared;

	// Check if the distance from P to C is less than or equal to the radius
	float dist_CP_squared = squared_distance(C, P);
	return dist_CP_squared <= radius * radius;
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

float Environment::calculateDensity(Vector2D point) {

	constexpr auto scalar = 1000;

	float density = 0.0f;
	const float mass = 1.0f;

	for (auto& particle : m_Particles) {
		float distance = sqrt(squared_distance(point, particle.m_Position));
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
		float distance = sqrt(squared_distance(point, m_Particles.at(i).m_Position));

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

float convertDensityToPressure(float density) {
	const float targetDensity = 0.05f;
	const float pressureConstant = 1.5f;

	float densityError = density - targetDensity;
	float pressure = pressureConstant * densityError;
	return pressure;
}

float calculateSharedPressure(float density1, float density2) {
	float pressure1 = convertDensityToPressure(density1);
	float pressure2 = convertDensityToPressure(density2);
	return (pressure1 + pressure2) / 2;
}

Vector2D getRandomDir() {
	// Seed the random number generator
	std::random_device rd;
	std::mt19937 gen(rd());

	// Create a uniform distribution for the angle (in radians)
	std::uniform_real_distribution<float> dist(0, 2 * M_PI); // Range: [0, 2 * pi]

	// Generate a random angle
	float angle = dist(gen);

	// Calculate the x and y components of the direction vector
	float x = std::cos(angle);
	float y = std::sin(angle);

	// Create and output the direction vector
	Vector2D direction = { x, y };

	return direction;
}

Vector2D Environment::calculatePressureForce(int particleIndex) {
	Vector2D pressureForce = Vector2D();
	const float mass = 1.0f;

	for (int i = 0; i < m_Particles.size(); i++) {

		if (i == particleIndex) {
			continue;
		}

		float distance = sqrt(squared_distance(m_Particles.at(particleIndex).m_Position, m_Particles.at(i).m_Position));
		if (distance < particleRadius) {
			int tt = 0;
		}
		Vector2D dir = distance < particleRadius ? getRandomDir() : (m_Particles.at(i).m_Position - m_Particles.at(particleIndex).m_Position) / distance;
		float slope = smoothingKernelDerivative(particleRadiusOfRepel, distance);

		float density = m_ParticleDensities.at(i);

		float sharedPressure = calculateSharedPressure(density, m_ParticleDensities.at(i));

		pressureForce += -sharedPressure * dir * slope * mass / density;
	}

	return pressureForce;
}

void Environment::update(float dt) {
	// Parallelize the loop using OpenMP
#pragma omp parallel for reduction(+:total)
	for (int i = 0; i < m_Particles.size(); i++) {

		//Vector2D futurePosition = m_Particles.at(i).m_Position + m_Particles.at(i).m_Velocity * dt;
		//m_ParticleDensities.at(i) = calculateDensity(futurePosition);

		m_ParticleDensities.at(i) = calculateDensity(m_Particles.at(i).m_Position);

	}
	for (int i = 0; i < m_Particles.size(); i++) {

		Vector2D pressureForce = calculatePressureForce(i);
		Vector2D pressureAcceleration = pressureForce / m_ParticleDensities.at(i);

		m_Particles.at(i).m_Velocity += pressureAcceleration * dt;

		m_Particles.at(i).update(dt);


		for (auto& obstacle : m_Obstacles) {
			if (check_line_segment_circle_intersection(obstacle.Point1, obstacle.Point2, m_Particles.at(i).m_Position, particleRadius)) {
				Vector2D normalVector = Math::calculateNormalVector(Math::calculateSlope(obstacle.Point1, obstacle.Point2));
				Vector2D reflectionVector = Math::calculateReflectionVector(m_Particles.at(i).m_Velocity, normalVector);

				// magnitude of reflection vector
				float magnitude = sqrt(reflectionVector.X * reflectionVector.X + reflectionVector.Y * reflectionVector.Y);

				// normalize the reflection vector
				/*reflectionVector.X /= magnitude;
				reflectionVector.Y /= magnitude;*/

				m_Particles.at(i).m_Velocity = reflectionVector * 0.1f;

				/*particle.m_Velocity.Y = -particle.m_Velocity.Y;*/

				m_Particles.at(i).m_Position = m_Particles.at(i).m_LastSafePosition;
			}
		}

		for (int j = 0; j < m_Particles.size(); j++) {
			if (i == j) {
				continue;
			}
			if (squared_distance(m_Particles.at(i).m_Position, m_Particles.at(j).m_Position) <= (particleRadius * particleRadius * 4)) {

				//Vector2D normalVector = Vector2D(m_Particles.at(j).m_Position.X - m_Particles.at(i).m_Position.X, m_Particles.at(j).m_Position.Y - m_Particles.at(i).m_Position.Y);

				////magnitude of normal vector
				//float magnitude = sqrt(normalVector.X * normalVector.X + normalVector.Y * normalVector.Y);

				//// normalize the normal vector
				//normalVector.X /= magnitude;
				//normalVector.Y /= magnitude;

				//Vector2D reflectionVector = Math::calculateReflectionVector(m_Particles.at(i).m_Velocity, normalVector);

				////m_Particles.at(i).m_Velocity = reflectionVector * 0.1f;


				//m_Particles.at(i).m_Position = m_Particles.at(i).m_LastSafePosition;
			}
			if (squared_distance(m_Particles.at(i).m_Position, m_Particles.at(j).m_Position) <= (particleRadius * particleRadius) * 4) {
				Vector2D normalVector = Vector2D(m_Particles.at(j).m_Position.X - m_Particles.at(i).m_Position.X, m_Particles.at(j).m_Position.Y - m_Particles.at(i).m_Position.Y);

				//magnitude of normal vector
				float magnitude = -1 * sqrt(normalVector.X * normalVector.X + normalVector.Y * normalVector.Y);

				// normalize the normal vector
				normalVector.X /= magnitude;
				normalVector.Y /= magnitude;

				//float power = smoothing_kernel(particleRadiusOfRepel - particleRadius, sqrt(squared_distance(m_Particles.at(i).m_Position, m_Particles.at(j).m_Position)) - particleRadius * 2);
				float power = 1;
				power *= particleRepulsionForce;

				m_Particles.at(i).m_Velocity += normalVector * power;

				m_Particles.at(j).m_Velocity -= normalVector * power;

				int x = 0;

				//particle.m_Position = particle.m_LastSafePosition;
			}
		}
	}
}
