#include "Particle.hpp"

#include "SDL.h"

#include <vector>

class Environment
{
public:
		Environment();
		~Environment();

		void draw(SDL_Renderer* renderer);

private:
	std::vector<Particle> m_particles;
};