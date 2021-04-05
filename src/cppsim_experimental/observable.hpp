
/**
 * @file pauli_operator.hpp
 * @brief Definition and basic functions for MultiPauliTerm
 */

#pragma once

#include <cassert>
#include <csim/stat_ops.hpp>
#include <csim/stat_ops_dm.hpp>
#include <iostream>
#include <regex>
#include <vector>

#ifdef _USE_GPU
#include <gpusim/stat_ops.h>
#endif

#include "type.hpp"

class QuantumStateBase;

enum {
    PAULI_ID_I = 0,
    PAULI_ID_X = 1,
    PAULI_ID_Y = 2,
    PAULI_ID_Z = 3,
};

class DllExport MultiQubitPauliOperator {
private:
    std::vector<UINT> _target_index;
    std::vector<UINT> _pauli_id;

public:
    MultiQubitPauliOperator(){};
    MultiQubitPauliOperator(const std::vector<UINT>& target_qubit_index_list,
        const std::vector<UINT>& pauli_id_list)
        : _target_index(target_qubit_index_list), _pauli_id(pauli_id_list){};

    virtual const std::vector<UINT>& get_pauli_id_list() const {
        return _pauli_id;
    }

    virtual const std::vector<UINT>& get_index_list() const {
        return _target_index;
    }

    explicit MultiQubitPauliOperator(std::string pauli_string) {
        std::string pattern = "([IXYZ])\\s*([0-9]+)\\s*";
        std::regex re(pattern);
        std::cmatch result;
        while (std::regex_search(pauli_string.c_str(), result, re)) {
            std::string pauli = result[1].str();
            UINT index = (UINT)std::stoul(result[2].str());
            _target_index.push_back(index);
            if (pauli == "I")
                _pauli_id.push_back(PAULI_ID_I);
            else if (pauli == "X")
                _pauli_id.push_back(PAULI_ID_X);
            else if (pauli == "Y")
                _pauli_id.push_back(PAULI_ID_Y);
            else if (pauli == "Z")
                _pauli_id.push_back(PAULI_ID_Z);
            else
                assert(false && "Error in regex");
            pauli_string = result.suffix();
        }
        assert(_target_index.size() == _pauli_id.size());
    }

    virtual ~MultiQubitPauliOperator(){};
    virtual void add_single_Pauli(UINT qubit_index, UINT pauli_type) {
        if (pauli_type >= 4)
            throw std::invalid_argument("pauli type must be any of 0,1,2,3");
        _target_index.push_back(qubit_index);
        _pauli_id.push_back(pauli_type);
    }

    virtual CPPCTYPE get_expectation_value(const QuantumStateBase* state) const;

    virtual CPPCTYPE get_transition_amplitude(const QuantumStateBase* state_bra,
        const QuantumStateBase* state_ket) const;
};

using PauliOperator = MultiQubitPauliOperator;

class DllExport Observable {
private:
    std::vector<MultiQubitPauliOperator> _pauli_terms;
    std::vector<CPPCTYPE> _coef_list;
    UINT _qubit_count;

public:
    Observable(){};

    /**
     * @param[in] qubit_count qubit数
     * @return Observableのインスタンス
     */
    explicit Observable(const UINT qubit_count) : _qubit_count(qubit_count){};

    virtual UINT get_term_count() const { return (UINT)_pauli_terms.size(); }

    /**
     * \~japanese-en
     * @return Observable のqubit数
     */
    virtual UINT get_qubit_count() const { return this->_qubit_count; }

    /**
     * \~japanese-en
     * @return Observable の行列表現の次元
     */
    virtual ITYPE get_state_dim() const { return (1ULL) << this->_qubit_count; }

    virtual std::pair<CPPCTYPE, MultiQubitPauliOperator> get_term(
        UINT index) const {
        return std::make_pair(_coef_list.at(index), _pauli_terms.at(index));
    }
    virtual void add_term(CPPCTYPE coef, MultiQubitPauliOperator op) {
        _coef_list.push_back(coef);
        _pauli_terms.push_back(op);
    }
    virtual void add_term(CPPCTYPE coef, std::string s) {
        _coef_list.push_back(coef);
        _pauli_terms.push_back(MultiQubitPauliOperator(s));
    }
    virtual void remove_term(UINT index) {
        _coef_list.erase(_coef_list.begin() + index);
        _pauli_terms.erase(_pauli_terms.begin() + index);
    }

    /**
     * \~japanese-en
     * ランダムなパウリ演算子を observable に追加する
     * @param [in] operator_count observable に追加するパウリ演算子数
     */
    void add_random_operator(const UINT operator_count);

    /**
     * \~japanese-en
     * GeneralQuantumOperatorのある量子状態に対応するエネルギー(期待値)を計算して返す
     *
     * @param[in] state 期待値をとるときの量子状態
     * @return 入力で与えた量子状態に対応するGeneralQuantumOperatorの期待値
     */
    virtual CPPCTYPE get_expectation_value(const QuantumStateBase* state) const;

    /**
     * \~japanese-en
     * GeneralQuantumOperatorによってある状態が別の状態に移る遷移振幅を計算して返す
     *
     * @param[in] state_bra 遷移先の量子状態
     * @param[in] state_ket 遷移前の量子状態
     * @return 入力で与えた量子状態に対応するGeneralQuantumOperatorの遷移振幅
     */
    virtual CPPCTYPE get_transition_amplitude(const QuantumStateBase* state_bra,
        const QuantumStateBase* state_ket) const;

    //     /**
    //      * \~japanese-en
    //      * Observable の基底状態の固有値を arnordi method により求める
    //      * (A - \mu I)
    //      の絶対値最大固有値を求めることで基底状態の固有値を求める．
    //      * @param[in] state 固有値を求めるための量子状態
    //      * @param[in] n_iter 計算の繰り返し回数
    //      * @param [in] mu 固有値をシフトするための係数
    //      * @return Observable の基底状態の固有値
    //      */
    //     virtual CPPCTYPE solve_ground_state_eigenvalue_by_arnoldi_method(
    //         QuantumStateBase* state, const UINT iter_count,
    //         const CPPCTYPE mu = 0.0) const;

    //     /**
    //      * \~japanese-en
    //      * Observable の基底状態の固有値を lanczos method により求める
    //      * (A - \mu I)
    //      の絶対値最大固有値を求めることで基底状態の固有値を求める．
    //      * @param[in] state 固有値を求めるための量子状態
    //      * @param[in] n_iter 計算の繰り返し回数
    //      * @param [in] mu 固有値をシフトするための係数
    //      * @return Observable の基底状態の固有値
    //      */
    //     CPPCTYPE solve_ground_state_eigenvalue_by_lanczos_method(
    //         QuantumStateBase* state, const UINT iter_count,
    //         const CPPCTYPE mu = 0.0) const;

    //     /**
    //      * \~japanese-en
    //      * Observable の基底状態の固有値を power method により求める
    //      * (A - \mu I)
    //      の絶対値最大固有値を求めることで基底状態の固有値を求める．
    //      * @param[in] state 固有値を求めるための量子状態
    //      * @param[in] n_iter 計算の繰り返し回数
    //      * @param [in] mu 固有値をシフトするための係数
    //      *  @return Observable の基底状態の固有値
    //      */
    //     CPPCTYPE
    //     solve_ground_state_eigenvalue_by_power_method(QuantumStateBase*
    //     state,
    //         const UINT iter_count, const CPPCTYPE mu = 0.0) const;

    //     /**
    //      * \~japanese-en
    //      * state_to_be_multiplied に Observable を作用させる．
    //      * 結果は dst_state に格納される．dst_state
    //      * はすべての要素を0に初期化してから計算するため，
    //      任意の状態を渡してよい．
    //      * @param [in] state_to_be_multiplied 作用を受ける状態
    //      * @param [in] dst_state 結果を格納する状態
    //      */
    //     void apply_to_state(QuantumStateBase* work_state,
    //         const QuantumStateBase& state_to_be_multiplied,
    //         QuantumStateBase* dst_state) const;

    // private:
    //     /**
    //      * \~japanese-en
    //      * solve_ground_state_eigenvalue_by_power_method の mu
    //      * のデフォルト値を計算する．
    //      */
    //     CPPCTYPE calculate_default_mu() const;
};
