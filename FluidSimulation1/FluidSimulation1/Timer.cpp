#include "Timer.hpp"

#include <iostream>

Timer* Timer::s_instance = nullptr;

float Timer::GetTime() {
	std::chrono::steady_clock::time_point current_time = std::chrono::steady_clock::now();
	double tick = std::chrono::duration_cast<std::chrono::microseconds>(current_time - m_lastTime).count();

	tick = tick / 1000000.0f;

	if (tick > 0.1f) {
		tick = 0.1f;
	}

	m_lastTime = current_time;

	return tick*4;
}
