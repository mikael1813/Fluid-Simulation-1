#include "Environment.hpp"

#include <SDL.h>

constexpr auto particleCount = 10;
constexpr auto particleRadius = 10;
constexpr auto particleDistance = 30;

constexpr auto SCREEN_WIDTH = 1280;
constexpr auto SCREEN_HEIGHT = 720;


Environment::Environment() {
	for (int i = 0; i < SCREEN_HEIGHT / particleRadius; i++) {
		std::vector<MatrixComponenets> row;
		for (int j = 0; j < SCREEN_WIDTH / particleRadius; j++) {
			row.push_back(MatrixComponenets());
		}
		m_InteractionsMatrix.push_back(row);
	}

	m_Particles = std::vector<Particle>{};

	for (int i = 0; i < particleCount; i++) {
		for (int j = 0; j < particleCount; j++) {
			m_Particles.push_back(Particle(200 + i * particleDistance, 200 + j * particleDistance));
		}
	}

	m_Obstacles.push_back(Surface2D(10, 10, 10, 700));
	m_Obstacles.push_back(Surface2D(10, 700, 1200, 700));
	m_Obstacles.push_back(Surface2D(1200, 10, 1200, 700));


}

Environment::~Environment()
{
}

void DrawCircle(SDL_Renderer* renderer, int centreX, int centreY, int radius)
{
	const int32_t diameter = (radius * 2);

	int32_t x = (radius - 1);
	int32_t y = 0;
	int32_t tx = 1;
	int32_t ty = 1;
	int32_t error = (tx - diameter);

	while (x >= y)
	{
		//  Each of the following renders an octant of the circle
		SDL_RenderDrawPoint(renderer, centreX + x, centreY - y);
		SDL_RenderDrawPoint(renderer, centreX + x, centreY + y);
		SDL_RenderDrawPoint(renderer, centreX - x, centreY - y);
		SDL_RenderDrawPoint(renderer, centreX - x, centreY + y);
		SDL_RenderDrawPoint(renderer, centreX + y, centreY - x);
		SDL_RenderDrawPoint(renderer, centreX + y, centreY + x);
		SDL_RenderDrawPoint(renderer, centreX - y, centreY - x);
		SDL_RenderDrawPoint(renderer, centreX - y, centreY + x);

		if (error <= 0)
		{
			++y;
			error += ty;
			ty += 2;
		}

		if (error > 0)
		{
			--x;
			tx += 2;
			error += (tx - diameter);
		}
	}
}

void drawSurface2D(SDL_Renderer* renderer, Surface2D surface) {
	SDL_RenderDrawLine(renderer, surface.Point1.X, surface.Point1.Y, surface.Point2.X, surface.Point2.Y);
}

void Environment::render(SDL_Renderer* renderer)
{
	for (auto& particle : m_Particles) {
		DrawCircle(renderer, particle.m_Position.X, particle.m_Position.Y, particleRadius);
	}

	for (auto& obstacle : m_Obstacles) {
		drawSurface2D(renderer, obstacle);
	}

	/*SDL_RenderDrawLine(renderer, 10, 10, 10, 700);

	SDL_RenderDrawLine(renderer, 10, 700, 1200, 700);

	SDL_RenderDrawLine(renderer, 1200, 10, 1200, 700);*/
}

// Function to compute the squared distance between two points
float squared_distance(const Point2D& p1, const Point2D& p2) {
	float dx = p1.X - p2.X;
	float dy = p1.Y - p2.Y;
	return dx * dx + dy * dy;
}

// Check if the line segment AB intersects the circle with center C and radius R
bool check_line_segment_circle_intersection(const Point2D& A, const Point2D& B, const Point2D& C, float radius) {
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
	Point2D P;
	P.X = A.X + (B.X - A.X) * dot_product / dist_AB_squared;
	P.Y = A.Y + (B.Y - A.Y) * dot_product / dist_AB_squared;

	// Check if the distance from P to C is less than or equal to the radius
	float dist_CP_squared = squared_distance(C, P);
	return dist_CP_squared <= radius * radius;
}

void Environment::update(float dt) {
	for (auto& particle : m_Particles) {
		particle.update(dt);

		for (auto& obstacle : m_Obstacles) {
			if (check_line_segment_circle_intersection(obstacle.Point1, obstacle.Point2, particle.m_Position, particleRadius)) {
				particle.m_Velocity.Y = -particle.m_Velocity.Y;
				particle.m_Position = particle.m_LastSafePosition;
			}
		}
	}
}