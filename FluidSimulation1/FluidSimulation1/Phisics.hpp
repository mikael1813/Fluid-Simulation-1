#pragma once

struct Point2D {
	Point2D() {
		X = 0;
		Y = 0;
	}
	Point2D(float x, float y) {
		X = x;
		Y = y;
	}
	float X;
	float Y;
};

struct Surface2D {

	Surface2D(Point2D point1, Point2D point2) {
		Point1 = point1;
		Point2 = point2;
	}

	Surface2D(float x1, float y1, float x2, float y2) {
		Point1 = Point2D(x1, y1);
		Point2 = Point2D(x2, y2);
	}

	Point2D Point1;
	Point2D Point2;
};

struct Vector2D {

	Vector2D() {
		X = 0;
		Y = 0;
	}

	Vector2D(float x, float y) {
		X = x;
		Y = y;
	}

	float X;
	float Y;
};

// Function to calculate the slope of a line given two points
double calculateSlope(Point2D a, Point2D b) {
	// Ensure x2 is not equal to x1 to avoid division by zero
	if (a.X == b.X) {
		std::cerr << "Error: Division by zero (x2 - x1 = 0)" << std::endl;
		return 0.0;
	}

	// Calculate the slope using the formula (y2 - y1) / (x2 - x1)
	return (b.Y - a.Y) / (b.X - a.X);
}