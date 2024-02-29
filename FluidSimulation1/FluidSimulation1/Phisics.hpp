#pragma once

#include <cmath> // Include <cmath> for handling infinity

struct Vector2D {

	Vector2D() {
		X = 0;
		Y = 0;
	}

	Vector2D(float x, float y) {
		X = x;
		Y = y;
	}

	// Function to calculate the magnitude of a 2D vector (x, y)
	float magnitude() {
		return std::sqrt(X * X + Y * Y);
	}

	// add 2 vectors
	Vector2D operator+(const Vector2D& other) const {
		return Vector2D(X + other.X, Y + other.Y);
	}

	// multiply by a scalar
	Vector2D operator*(float scalar) const {
		return Vector2D(X * scalar, Y * scalar);
	}

	// add 1 vector to current vector
	Vector2D& operator+=(const Vector2D& other) {
		X += other.X;
		Y += other.Y;
		return *this;
	}

	// substract 2 vectors
	Vector2D operator-(const Vector2D& other) const {
		return Vector2D(X - other.X, Y - other.Y);
	}

	// subtract a vector from the current vector
	Vector2D& operator-=(const Vector2D& other) {
		X -= other.X;
		Y -= other.Y;
		return *this;
	}

	// divide by a scalar
	Vector2D operator/(float scalar) const {
		return Vector2D(X / scalar, Y / scalar);
	}

	// float times vector
	friend Vector2D operator*(float scalar, const Vector2D& vector) {
		return vector * scalar;
	}

	// minus vector
	Vector2D operator-() const {
		return Vector2D(-X, -Y);
	}

	float X;
	float Y;
};

struct Surface2D {

	Surface2D(Vector2D point1, Vector2D point2) {
		Point1 = point1;
		Point2 = point2;
	}

	Surface2D(float x1, float y1, float x2, float y2) {
		Point1 = Vector2D(x1, y1);
		Point2 = Vector2D(x2, y2);
	}

	Vector2D Point1;
	Vector2D Point2;
};



class Math {

public:

	// Function to calculate the slope of a line given two points
	static double calculateSlope(Vector2D a, Vector2D b) {
		// Ensure x2 is not equal to x1 to avoid division by zero
		if (a.X == b.X) {
			/*std::cerr << "Error: Division by zero (x2 - x1 = 0)" << std::endl;*/
			return INFINITY;
		}

		// Calculate the slope using the formula (y2 - y1) / (x2 - x1)
		return (b.Y - a.Y) / (b.X - a.X);
	}

	// Function to calculate the normal vector given the slope of the surface line
	static Vector2D calculateNormalVector(double surfaceLineSlope) {
		// Calculate the slope of the perpendicular line
		double perpendicularLineSlope;

		if (surfaceLineSlope == 0.0) {
			// Handle the case when the surface line is horizontal (slope is 0)
			perpendicularLineSlope = INFINITY; // Treat the perpendicular line slope as infinity
		}
		else {
			// Calculate the slope of the perpendicular line (negative reciprocal)
			perpendicularLineSlope = -1.0 / surfaceLineSlope;
		}

		// The normal vector is represented by the coefficients (1, m), where m is the perpendicular line slope
		double normalX = 1.0;
		double normalY = perpendicularLineSlope;

		// Calculate the magnitude of the normal vector
		double magnitude = sqrt(normalX * normalX + normalY * normalY);

		// Normalize the components to obtain the direction of the normal vector
		normalX /= magnitude;
		normalY /= magnitude;

		/*if (abs(normalY) > 1.0) {
			normalX /= normalY;
			if (normalY == INFINITY) {
				normalY = 1.0;
			}
			else {
				normalY /= normalY;
			}
		}*/

		return Vector2D(normalX, normalY);
	}

	// Function to calculate the reflection vector given the incident vector and the normal vector
	static Vector2D calculateReflectionVector(const Vector2D& incidentVector, const Vector2D& normalVector) {
		// Calculate the dot product of the incident vector and the normal vector
		double dotProduct = incidentVector.X * normalVector.X + incidentVector.Y * normalVector.Y;

		// Calculate the reflection vector
		Vector2D reflectionVector = Vector2D(incidentVector.X - 2 * dotProduct * normalVector.X,
			incidentVector.Y - 2 * dotProduct * normalVector.Y);

		return reflectionVector;
	}
};