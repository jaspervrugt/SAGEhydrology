// ========================================================================
// SAGEhydrology | Model source
// File: Native/models/CaRTPlus/cartplus_routing.hpp
//
// CaRT / CaRT+ numerical model and supporting routines.
//
// SAGE adaptations: Jasper A. Vrugt, University of California, Irvine.
// Upstream model: Xia (2026). Preserve the adjacent license notices.
// ========================================================================

// CaRT+ routing updates adapted from Xia's MIT-licensed reference.
// Upstream revision 83682732728646669a190835d95ab3e60a8b5cb5.
// See LICENSE in this directory.
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <type_traits>
#include <vector>

#include "cartplus_tangent.hpp"

namespace cartplus {

template <std::size_t Length, typename Scalar = double>
class UnitHydrograph {
public:
    /**
     * UnitHydrograph: Initialize the numerical scalar, routing object, or
     * solver workspace.
     *
     * INPUT ARGUMENTS:
     * time_base : Hydrograph time base (days).
     * dt : Integration or observation interval duration (days).
     * skip_zero_inputs : Whether a zero input can use the sparse routing
     * update.
     *
     * OUTPUT ARGUMENTS:
     * Numerical result or flux structure in the argument time base;
     * non-const state references are advanced in place.
     */
    explicit UnitHydrograph(Scalar time_base, double dt = 1.0,
                            bool skip_zero_inputs = false)
        : weights_(dt == 1.0 ? Length : static_cast<std::size_t>(
              std::ceil(primal(time_base) * (Length == 6 ? 1.0 : 2.0) / dt))
                  ),
          pending_(weights_.size()), skip_zero_inputs_(skip_zero_inputs) {
        static_assert(Length == 6 || Length == 13);
        const auto cumulative = [time_base](Scalar lag) {
            const Scalar relative = lag / time_base;
            if constexpr (Length == 6) {
                return pow(minimum(relative, 1.0), 2.5);
            }
            return relative <= 1.0 ? 0.5 * pow(relative, 2.5)
                : 1.0 - 0.5 * pow(maximum(2.0-relative, 0.0), 2.5);
        };
        for (std::size_t i = 0; i < weights_.size(); ++i) {
            weights_[i] = cumulative((i + 1) * dt) - cumulative(i * dt);
            weight_sum_ += weights_[i];
        }
        if constexpr (std::is_same_v<Scalar, Tangent>) {
            for (std::size_t j = 0; j < parameter_count; ++j) {
                if (time_base.slope[j] != 0.0) weight_seeds_.push_back(j);
            }
        }
    }

    /**
     * advance: Advance routing memory and return its outlet depth.
     *
     * INPUT ARGUMENTS:
     * inflow : Input series to be routed, in interval depths (mm).
     *
     * OUTPUT ARGUMENTS:
     * Numerical result or flux structure in the argument time base;
     * non-const state references are advanced in place.
     */
    Scalar advance(Scalar inflow) {
        if constexpr (std::is_same_v<Scalar, Tangent>) {
            if (skip_zero_inputs_) return advance_subdaily(inflow);
        }
        const Scalar released = pending_[0] + weights_[0] * inflow;
        for (std::size_t i = 1; i < weights_.size(); ++i) {
            if constexpr (std::is_same_v<Scalar, Tangent>) {
                // Routing weights depend only on the time-base seeds.
                const auto& weight = weights_[i];
                const auto& next = pending_[i];
                auto& previous = pending_[i - 1];
                previous.value = next.value + weight.value * inflow.value;
                for (std::size_t j = 0; j < parameter_count; ++j) {
                    previous.slope[j] = next.slope[j]
                        + weight.value * inflow.slope[j];
                }
                for (auto j : weight_seeds_) {
                    previous.slope[j] += weight.slope[j] * inflow.value;
                }
            } else {
                pending_[i - 1] = pending_[i] + weights_[i] * inflow;
            }
        }
        return released;
    }

    /**
     * storage_change: Compute retained routing water from input minus
     * output.
     *
     * INPUT ARGUMENTS:
     * inflow : Input series to be routed, in interval depths (mm).
     * outflow : Water leaving the routing operator.
     *
     * OUTPUT ARGUMENTS:
     * Numerical result or flux structure in the argument time base;
     * non-const state references are advanced in place.
     */
    Scalar storage_change(Scalar inflow, Scalar outflow) const {
        return weight_sum_ * inflow - outflow;
    }

private:
    /**
     * advance_subdaily: Advance a circular routing buffer for one subdaily
     * input.
     *
     * INPUT ARGUMENTS:
     * inflow : Input series to be routed, in interval depths (mm).
     *
     * OUTPUT ARGUMENTS:
     * Numerical result or flux structure in the argument time base;
     * non-const state references are advanced in place.
     */
    Scalar advance_subdaily(const Scalar& inflow) {
        const Scalar released = pending_[head_] + weights_[0] * inflow;
        pending_[head_] = Scalar{};
        if (++head_ == pending_.size()) head_ = 0;

        // No contributions to future intervals from an identically zero
        // input. A zero primal with nonzero sensitivities still
        // contributes.
        if (inflow.value == 0.0 &&
            std::all_of(inflow.slope.begin(), inflow.slope.end(),
                        [](double slope) { return slope == 0.0; })) {
            return released;
        }

        // Accumulate in the same time order into fixed slots. Rotating the
        // head avoids shifting a full hydrograph of tangent vectors per
        // step.
        std::size_t slot = head_;
        for (std::size_t i = 1; i < weights_.size(); ++i) {
            auto& next = pending_[slot];
            const auto& weight = weights_[i];
            next.value += weight.value * inflow.value;
            for (std::size_t j = 0; j < parameter_count; ++j) {
                next.slope[j] += weight.value * inflow.slope[j];
            }
            for (auto j : weight_seeds_) {
                next.slope[j] += weight.slope[j] * inflow.value;
            }
            if (++slot == pending_.size()) slot = 0;
        }
        return released;
    }

    std::vector<Scalar> weights_;
    std::vector<Scalar> pending_;
    std::vector<std::size_t> weight_seeds_;
    Scalar weight_sum_{};
    bool skip_zero_inputs_;
    std::size_t head_{};
};

template <typename Scalar>
struct QuickFluxes {
    Scalar discharge;
    Scalar exchange;
};

/**
 * advance_quick: Advance the quick store for the selected observation.
 *
 * INPUT ARGUMENTS:
 * storage : Water store (mm); reference arguments are updated in place.
 * routed : Hydrograph-routed input to the quick store.
 * bypass : Bypass input after hydrograph routing.
 * release_rate : Quick-store release coefficient.
 * exchange_rate : Inter-store exchange coefficient.
 * capacity : Store or infiltration capacity in this operator.
 * connectivity : Dimensionless connected-area fraction.
 * storage_exponent : Storage exponent used by this numerical operator.
 * direct_output : Optional output pointer receiving direct-branch
 * discharge.
 *
 * OUTPUT ARGUMENTS:
 * Numerical result or flux structure in the argument time base; non-const
 * state references are advanced in place.
 */
template <typename Scalar>
inline QuickFluxes<Scalar> advance_quick(
    Scalar& storage, Scalar routed, Scalar bypass, Scalar release_rate,
    Scalar exchange_rate, Scalar capacity, Scalar connectivity,
    Scalar storage_exponent, Scalar* direct_output = nullptr
) {
    const Scalar relative = maximum(storage, 0.0) / 30.0;
    const Scalar exchange = exchange_rate * relative * relative * relative
        * sqrt(relative);
    const Scalar available = maximum(storage + routed + 0.9 * exchange, 0.0)
        ;
    const Scalar direct = maximum(bypass + 0.1 * exchange, 0.0);
    const Scalar bounded_log_multiplier = 6.0 * tanh(
        (connectivity + storage_exponent
            * log(maximum(available, 1e-3) / capacity)) / 6.0
    );
    const Scalar coefficient = release_rate * exp(bounded_log_multiplier);
    const Scalar squared = available * available;
    const Scalar next = available * exp(
        -0.25 * log1p(4.0 * coefficient * squared * squared)
    );
    if (direct_output) *direct_output = direct;
    const QuickFluxes<Scalar> fluxes{
        available - next + direct,
        (available - storage - routed) + (direct - bypass)
    };
    storage = next;
    return fluxes;
}

}  // namespace cartplus
