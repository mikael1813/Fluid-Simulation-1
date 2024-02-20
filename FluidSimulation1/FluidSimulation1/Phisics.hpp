#pragma once

struct Point {
	Point() {
		X = 0;
		Y = 0;
	}
	Point(float x, float y) {
		X = x;
		Y = y;
	}
	float X;
	float Y;
};

struct Surface2D {

	Surface2D(Point point1, Point point2) {
		Point1 = point1;
		Point2 = point2;
	}

	Surface2D(float x1, float y1, float x2, float y2) {
		Point1 = Point(x1, y1);
		Point2 = Point(x2, y2);
	}

	Point Point1;
	Point Point2;
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