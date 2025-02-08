#include <fft.h>
#include <cmath>
#include <vector>
#include <fftw3.h>

class DctCalculator::Impl {
public:
    Impl(size_t width, std::vector<double>* input, std::vector<double>* output)
        : width_(width), input_(input), output_(output) {
        if (width_ == 0 || !input_ || !output_ || input_->size() != width_ * width_ ||
            output_->size() != width_ * width_) {
            throw std::invalid_argument("Invalid args");
        }

        plan_ = fftw_plan_r2r_2d(static_cast<int>(width_), static_cast<int>(width_), input_->data(),
                                 output_->data(), FFTW_REDFT01, FFTW_REDFT01, FFTW_MEASURE);

        if (!plan_) {
            throw std::runtime_error("FFTW plan creation failed");
        }
    }

    ~Impl() {
        if (plan_) {
            fftw_destroy_plan(plan_);
        }
    }

    void Inverse() {
        for (int u = 0; u < static_cast<int>(width_); ++u) {
            for (int v = 0; v < static_cast<int>(width_); ++v) {
                (*input_)[u * width_ + v] /= 16.0;

                if (u == 0) {
                    (*input_)[u * width_ + v] *= sqrt(2.0);
                }
                if (v == 0) {
                    (*input_)[u * width_ + v] *= sqrt(2.0);
                }
            }
        }

        fftw_execute(plan_);
    }

private:
    size_t width_;
    std::vector<double>* input_;
    std::vector<double>* output_;
    fftw_plan plan_ = nullptr;
};

DctCalculator::DctCalculator(size_t width, std::vector<double>* input, std::vector<double>* output)
    : impl_(new Impl(width, input, output)) {
}

void DctCalculator::Inverse() {
    impl_->Inverse();
}

DctCalculator::~DctCalculator() = default;
