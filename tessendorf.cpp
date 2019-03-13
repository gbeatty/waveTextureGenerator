#define _USE_MATH_DEFINES
#include <cmath>
#include "fftw-3.3.5/fftw3.h"
#include "tessendorf.h"

std::mutex tessendorf::fftCreationMutex;

tessendorf::tessendorf(float amplitude, float speed, Cartesian3 direction, float choppiness, float time, float phaseDuration, int resX, int resY, float scaleX, float scaleY, float waveSizeLimit, unsigned long rngSeed)
{
    // Parameters.
    A = amplitude;
    V = speed;
    w_hat = direction.Normalize();
    t = time;
    lambda = choppiness;
    M = resX;
    N = resY;
    Lx = scaleX;
    Ly = scaleY;
    l = waveSizeLimit;
    T = phaseDuration;
    omega_0 = 2. * M_PI / T;
    seed = rngSeed;
    engine = std::mt19937();

    // Precalculate known constants.
    P_h__L = pow(V, 2) / GRAVITY;
    P_h__l_2 = pow(l, 2);
}

tessendorf::~tessendorf()
{

}

float tessendorf::omega(Cartesian3 k)
{
    return floor(sqrt(GRAVITY * k.Magnitude()) / omega_0) * omega_0;
}

float tessendorf::P_h(Cartesian3 k)
{
    float k_length = k.Magnitude();

    if (k_length < DBL_EPSILON) {
        return 0.; // Avoid divison by zero error.
    }

    Cartesian3 k_hat = k.Normalize();

    float nomin = exp(-1. / pow(k_length * P_h__L, 2));
    float denom = pow(k_length, 4);
    float scale = exp(-pow(k_length, 2) * P_h__l_2);

    return A * nomin / denom * pow(k_hat.Dot(w_hat), 2.0) * scale;
}

complex tessendorf::h_tilde_0(Cartesian3 k)
{
    return complex(dist(engine), dist(engine)) * (float)sqrt(P_h(k) / 2.);
}

complex tessendorf::h_tilde(Cartesian3 k)
{
    complex h_tilde_0_k = h_tilde_0(k);
    complex h_tilde_0_k_star = h_tilde_0(-k);

    float omega_k_t = omega(k) * t;

    float cos_omega_k_t = cos(omega_k_t);
    float sin_omega_k_t = sin(omega_k_t);

    complex c0(cos_omega_k_t, sin_omega_k_t);
    complex c1(cos_omega_k_t, -sin_omega_k_t);

    return h_tilde_0_k * c0 + h_tilde_0_k_star * c1;
}

std::vector<VertexData> tessendorf::simulate()
{
    engine.seed(seed);

    const size_t fftSize = M * N;

    std::vector<VertexData> vertices(fftSize);

    std::vector<complex> h_tildes_in(fftSize);
    std::vector<complex> displacement_x_in(fftSize);
    std::vector<complex> displacement_y_in(fftSize);
    std::vector<complex> gradient_h_in(fftSize);
    std::vector<complex> gradient_x_in(fftSize);
    std::vector<complex> gradient_y_in(fftSize);

    std::vector<complex> h_tildes_out(fftSize);
    std::vector<complex> displacement_x_out(fftSize);
    std::vector<complex> displacement_y_out(fftSize);
    std::vector<complex> gradient_h_out(fftSize);
    std::vector<complex> gradient_x_out(fftSize);
    std::vector<complex> gradient_y_out(fftSize);

    for (int m = 0; m < M; m++) {
        for (int n = 0; n < N; n++) {
            int index = m * N + n;

            int m_ = m - M / 2;  // m coord offsetted.
            int n_ = n - N / 2; // n coord offsetted.

            Cartesian3 k(2. * M_PI * n_ / Lx, 2. * M_PI * m_ / Ly, 0.);


            complex h_tilde_k = h_tilde(k);
            h_tildes_in[index] = h_tilde_k;

            Cartesian3 k_hat = k.Normalize();
            displacement_x_in[index] = complex(0., -k_hat.x / k_hat.Magnitude()) * h_tilde_k; // Displacement by equation (29).
            displacement_y_in[index] = complex(0., -k_hat.y / k_hat.Magnitude()) * h_tilde_k;



            // Gradient by equation (20).
            gradient_h_in[index] = complex(0., k.z) * h_tilde_k;
            gradient_x_in[index] = complex(0., k.x) * h_tilde_k;
            gradient_y_in[index] = complex(0., k.y) * h_tilde_k;
        }
    }

    // fftwf_plan_dft_1d is not thread safe. Guard it with a lock
    std::unique_lock<std::mutex> lock(fftCreationMutex);
    fftwf_plan p_h = fftwf_plan_dft_1d(fftSize, reinterpret_cast<fftwf_complex*>(&h_tildes_in[0]), reinterpret_cast<fftwf_complex*>(&h_tildes_out[0]), FFTW_FORWARD, FFTW_ESTIMATE);
    fftwf_plan p_dx = fftwf_plan_dft_1d(fftSize, reinterpret_cast<fftwf_complex*>(&displacement_x_in[0]), reinterpret_cast<fftwf_complex*>(&displacement_x_out[0]), FFTW_FORWARD, FFTW_ESTIMATE);
    fftwf_plan p_dy = fftwf_plan_dft_1d(fftSize, reinterpret_cast<fftwf_complex*>(&displacement_y_in[0]), reinterpret_cast<fftwf_complex*>(&displacement_y_out[0]), FFTW_FORWARD, FFTW_ESTIMATE);
    fftwf_plan p_gradient_h = fftwf_plan_dft_1d(fftSize, reinterpret_cast<fftwf_complex*>(&gradient_h_in[0]), reinterpret_cast<fftwf_complex*>(&gradient_h_out[0]), FFTW_FORWARD, FFTW_ESTIMATE);
    fftwf_plan p_gradient_x = fftwf_plan_dft_1d(fftSize, reinterpret_cast<fftwf_complex*>(&gradient_x_in[0]), reinterpret_cast<fftwf_complex*>(&gradient_x_out[0]), FFTW_FORWARD, FFTW_ESTIMATE);
    fftwf_plan p_gradient_y = fftwf_plan_dft_1d(fftSize, reinterpret_cast<fftwf_complex*>(&gradient_y_in[0]), reinterpret_cast<fftwf_complex*>(&gradient_y_out[0]), FFTW_FORWARD, FFTW_ESTIMATE);
    lock.unlock();

    // fftwf_execute is thread safe
    fftwf_execute(p_h);
    fftwf_execute(p_dx);
    fftwf_execute(p_dy);
    fftwf_execute(p_gradient_h);
    fftwf_execute(p_gradient_x);
    fftwf_execute(p_gradient_y);

    // fftwf_destroy_plan not thread safe
    lock.lock();
    fftwf_destroy_plan(p_h);
    fftwf_destroy_plan(p_dx);
    fftwf_destroy_plan(p_dy);
    fftwf_destroy_plan(p_gradient_h);
    fftwf_destroy_plan(p_gradient_x);
    fftwf_destroy_plan(p_gradient_y);
    lock.unlock();

    float signs[2] = { -1., 1. };

    for (int m = 0; m < M; m++) {
        for (int n = 0; n < N; n++) {
            int index = m * N + n;
            int sign = signs[(m + n) & 1]; // Sign-flip all of the odd coefficients.

            float m_ = m - M / 2.0f;  // m coord offsetted.
            float n_ = n - N / 2.0f;  // n coord offsetted.

            const Cartesian3 displacement(
                n_ * Lx / N + real(displacement_x_out[index]) * lambda * sign,
                m_ * Ly / M + real(displacement_y_out[index]) * lambda * sign,
                real(h_tildes_out[index]) * sign);

            const Cartesian3 gradient(
                real(gradient_x_out[index]) * sign,
                real(gradient_y_out[index]) * sign,
                real(gradient_h_out[index]) * sign);

            Cartesian3 normal(-gradient.x, -gradient.y, 1);
			normal = normal.Normalize();

            vertices[index] = VertexData(displacement, normal);

        }
    }

    return vertices;
}