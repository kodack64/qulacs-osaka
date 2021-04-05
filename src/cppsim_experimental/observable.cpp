#include "observable.hpp"

#include "state.hpp"

CPPCTYPE MultiQubitPauliOperator::get_expectation_value(
    const QuantumStateBase* state) const {
    if (state->get_device_type() == DeviceType::Cpu) {
        if (state->is_state_vector()) {
            return expectation_value_multi_qubit_Pauli_operator_partial_list(
                this->_target_index.data(), this->_pauli_id.data(),
                (UINT)this->_target_index.size(), state->data_c(), state->dim);
        } else {
            return dm_expectation_value_multi_qubit_Pauli_operator_partial_list(
                this->_target_index.data(), this->_pauli_id.data(),
                (UINT)this->_target_index.size(), state->data_c(), state->dim);
        }
    } else if (state->get_device_type() == DeviceType::Gpu) {
#ifdef _USE_GPU
        if (state->is_state_vector()) {
            return expectation_value_multi_qubit_Pauli_operator_partial_list_host(
                this->get_index_list().data(), this->get_pauli_id_list().data(),
                (UINT)this->get_index_list().size(), state->data(), state->dim,
                state->get_cuda_stream(), state->device_number);
        } else {
            throw std::runtime_error(
                "Get expectation value for DensityMatrix on GPU is not "
                "supported");
        }
#else
        throw std::invalid_argument("GPU is not supported in this build");
#endif
    } else {
        throw std::invalid_argument("Unsupported device type");
    }
}

CPPCTYPE MultiQubitPauliOperator::get_transition_amplitude(
    const QuantumStateBase* state_bra,
    const QuantumStateBase* state_ket) const {
    if (state_bra->get_device_type() != state_ket->get_device_type())
        throw std::invalid_argument("Device type is different");
    if (state_bra->is_state_vector() != state_ket->is_state_vector())
        throw std::invalid_argument("is_state_vector is not matched");
    if (state_bra->dim != state_ket->dim)
        throw std::invalid_argument("state_bra->dim != state_ket->dim");

    if (state_bra->get_device_type() == DeviceType::Cpu) {
        if (state_bra->is_state_vector()) {
            return transition_amplitude_multi_qubit_Pauli_operator_partial_list(
                this->_target_index.data(), this->_pauli_id.data(),
                (UINT)this->_target_index.size(), state_bra->data_c(),
                state_ket->data_c(), state_bra->dim);
        } else {
            throw std::invalid_argument(
                "TransitionAmplitude for density matrix is not implemtend");
        }
    } else if (state_bra->get_device_type() == DeviceType::Gpu) {
#ifdef _USE_GPU
        if (state_bra->is_state_vector()) {
            return transition_amplitude_multi_qubit_Pauli_operator_partial_list_host(
                this->get_index_list().data(), this->get_pauli_id_list().data(),
                (UINT)this->get_index_list().size(), state_bra->data(),
                state_ket->data(), state_bra->dim, state_bra->get_cuda_stream(),
                state_bra->device_number);
        } else {
            throw std::runtime_error(
                "Get expectation value for DensityMatrix on GPU is not "
                "supported");
        }
#else
        throw std::invalid_argument("GPU is not supported in this build");
#endif
    } else {
        throw std::invalid_argument("Unsupported device");
    }
}

CPPCTYPE Observable::get_expectation_value(
    const QuantumStateBase* state) const {
    CPPCTYPE sum = 0;
    for (UINT index = 0; index < this->_pauli_terms.size(); ++index) {
        sum += this->_coef_list.at(index) *
               this->_pauli_terms.at(index).get_expectation_value(state);
    }
    return sum;
}

CPPCTYPE Observable::get_transition_amplitude(const QuantumStateBase* state_bra,
    const QuantumStateBase* state_ket) const {
    CPPCTYPE sum = 0;
    for (UINT index = 0; index < this->_pauli_terms.size(); ++index) {
        sum += this->_coef_list.at(index) *
               this->_pauli_terms.at(index).get_transition_amplitude(
                   state_bra, state_ket);
    }
    return sum;
}
