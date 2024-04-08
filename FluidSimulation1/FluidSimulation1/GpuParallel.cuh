#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "Particle.hpp"

#include <vector>
#include "InteractionMatrixClass.hpp"


void GpuParallelUpdateParticleDensities(std::vector<Particle*> particles, InteractionMatrixClass* interactionMatrix, int particleRadiusOfRepel);