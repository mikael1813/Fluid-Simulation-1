#include "GpuParallel.cuh"

float smoothingKernel(float radius, float distance) {
	if (distance >= radius) {
		return 0.0f;
	}

	float x = (radius - distance) / radius;
	return x * x;
}

__global__ void updateParticleDensitiesKernel(Particle* particles, InteractionMatrixClass* interactionMatrix, int particleRadiusOfRepel) {

	int index = threadIdx.x + blockIdx.x * blockDim.x;
	Particle* particle = &particles[index];

	Vector2D point = particle->m_PredictedPosition;

	std::vector<Particle*> particlesInCell = interactionMatrix->getParticlesInCell(point, particleRadiusOfRepel);

	constexpr auto scalar = 1000;

	float density = 0.0f;
	const float mass = 1.0f;

	for (int i = 0; i < particlesInCell.size(); i++) {
		float distance = sqrt(Math::squared_distance(point, particle->m_PredictedPosition));
		float influence = smoothingKernel(particleRadiusOfRepel, distance);
		density += mass * influence;
	}

	float volume = 3.1415f * pow(particleRadiusOfRepel, 2);

	density = density / volume * scalar;

	particle->m_Density = density;
}



void GpuParallelUpdateParticleDensities(std::vector<Particle*> particles, InteractionMatrixClass* interactionMatrix, int particleRadiusOfRepel) {


	std::vector<Particle*>* cudaParticles;
	InteractionMatrixClass* cudaInteractionMatrix;
	int* cudaParticleRadiusOfRepel;

	cudaMalloc(&cudaParticles, particles.size() * sizeof(Particle*));
	cudaMalloc(&cudaInteractionMatrix, sizeof(InteractionMatrixClass));
	cudaMalloc(&cudaParticleRadiusOfRepel, sizeof(int));

	cudaMemcpy(cudaParticles, particles.data(), particles.size() * sizeof(Particle*), cudaMemcpyHostToDevice);
	cudaMemcpy(cudaInteractionMatrix, interactionMatrix, sizeof(InteractionMatrixClass), cudaMemcpyHostToDevice);
	cudaMemcpy(cudaParticleRadiusOfRepel, &particleRadiusOfRepel, sizeof(int), cudaMemcpyHostToDevice);

	updateParticleDensitiesKernel << <1, particles.size() >> > (cudaParticles, cudaInteractionMatrix, *cudaParticleRadiusOfRepel);

	int x = 0;
}