
class Particle {
public:
	Particle(float x, float y) : m_position_x(x), m_position_y(y) {}
	float m_position_x, m_position_y;
private:
	float m_velocity_x, m_velocity_y;
	float visible_radius = 2.0f;


};