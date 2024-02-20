#pragma once

#include <chrono>

class Timer {
public:

	static Timer* GetInstance() {
		return s_instance = (s_instance != nullptr) ? s_instance : new Timer();
	}

	float GetTime();
private:
	static Timer* s_instance;

	std::chrono::steady_clock::time_point m_lastTime = std::chrono::steady_clock::now();

	Timer() {}
};
