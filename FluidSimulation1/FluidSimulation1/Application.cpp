#include "Application.hpp"
#include "Timer.hpp"

#include <iostream>

Application::Application()
{

	m_is_running = true;


	m_window = SDL_CreateWindow("SDL2 Window",
		SDL_WINDOWPOS_CENTERED,
		SDL_WINDOWPOS_CENTERED,
		1280, 720,
		SDL_WINDOW_RESIZABLE);

	if (!m_window)
	{
		std::cout << "Failed to create window\n";
		std::cout << "SDL2 Error: " << SDL_GetError() << "\n";
		return;
	}

	/*m_window_surface = SDL_GetWindowSurface(m_window);

	if (!m_window_surface)
	{
		std::cout << "Failed to get window's surface\n";
		std::cout << "SDL2 Error: " << SDL_GetError() << "\n";
		return;
	}*/

	m_renderer = SDL_CreateRenderer(
		m_window,
		-1,
		SDL_RENDERER_ACCELERATED |
		SDL_RENDERER_PRESENTVSYNC);
	if (m_renderer == nullptr) {
		SDL_Log("Failed to create Renderer: %s", SDL_GetError());
		return;
	}

	m_environment = new Environment();
}

Application::~Application()
{
	SDL_FreeSurface(m_window_surface);
	SDL_DestroyWindow(m_window);

}

void Application::events() {

	while (SDL_PollEvent(&m_window_event) > 0)
	{
		switch (m_window_event.type)
		{
		case SDL_QUIT:
			m_is_running = false;
			break;
		case SDL_MOUSEBUTTONDOWN:
			mousePress(m_window_event.button);
			break;
		case SDL_MOUSEBUTTONUP:
			if (m_window_event.button.button == SDL_BUTTON_LEFT) {
				//this->chessWindow->releaseLeftClick();
			}
			break;
		}
		switch (m_window_event.window.event) {
		case SDL_WINDOWEVENT_SIZE_CHANGED:
			m_window_surface = SDL_GetWindowSurface(m_window);
			break;
		}

	}
}

void Application::loop()
{
	while (m_is_running) {
		events();

		double deltaTime = Timer::GetInstance()->GetTime();

		//std::cout << deltaTime << std::endl;

		update(deltaTime);

		render();
	}
}

void Application::mousePress(SDL_MouseButtonEvent& b) {
	if (b.button == SDL_BUTTON_LEFT) {

	}
	if (b.button)
		if (b.button == SDL_BUTTON_RIGHT) {
			int x = b.x;
			int y = b.y;
		}
}

void Application::render()
{
	/*SDL_UpdateWindowSurface(m_window);
	SDL_FillRect(m_window_surface, NULL, SDL_MapRGB(m_window_surface->format, 0, 0, 0));*/
	SDL_RenderClear(m_renderer);

	SDL_SetRenderDrawColor(m_renderer, 255, 255, 255, 255);
	//SDL_RenderClear(m_renderer);

	//m_environment->render(m_renderer);

	SDL_SetRenderDrawColor(m_renderer, 0, 0, 0, 255);

	SDL_RenderPresent(m_renderer);
	//this->chessWindow->draw(m_window_surface);
}

void Application::update(float dt)
{
	m_environment->update(dt);
}
