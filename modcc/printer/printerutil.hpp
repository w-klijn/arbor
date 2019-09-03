#pragma once

// Convenience routines/helpers for source printers.

#include <ostream>
#include <string>
#include <vector>

#include "blocks.hpp"
#include "error.hpp"
#include "expression.hpp"
#include "module.hpp"

std::vector<std::string> namespace_components(const std::string& qualified_namespace);

inline const char* arb_header_prefix() {
    static const char* prefix = "arbor/";
    return prefix;
}

// TODO: this function will be obsoleted once arbor private/public headers are
// properly split.

inline const char* arb_private_header_prefix() {
    static const char* prefix = "";
    return prefix;
}

struct namespace_declaration_open {
    const std::vector<std::string>& ids;
    namespace_declaration_open(const std::vector<std::string>& ids): ids(ids) {}

    friend std::ostream& operator<<(std::ostream& o, const namespace_declaration_open& n) {
        for (auto& id: n.ids) {
            o << "namespace " << id << " {\n";
        }
        return o;
    }
};

struct namespace_declaration_close {
    const std::vector<std::string>& ids;
    namespace_declaration_close(const std::vector<std::string>& ids): ids(ids) {}

    friend std::ostream& operator<<(std::ostream& o, const namespace_declaration_close& n) {
        for (auto i = n.ids.rbegin(); i!=n.ids.rend(); ++i) {
            o << "} // namespace " << *i << "\n";
        }
        return o;
    }
};

// Enum representation:

inline const char* module_kind_str(const Module& m) {
    return m.kind()==moduleKind::density?
        "::arb::mechanismKind::density":
        "::arb::mechanismKind::point";
}

// Check expression non-null and scoped, or else throw.

inline void assert_has_scope(Expression* expr, const std::string& context) {
    return
        !expr? throw compiler_exception("missing expression for "+context):
        !expr->scope()? throw compiler_exception("printer invoked before semantic pass for "+context):
        void();
}


// Scope query functions:

// All local variables in scope with `is_indexed()` true.
std::vector<LocalVariable*> indexed_locals(scope_ptr scope);

// All local variables in scope with `is_arg()` and `is_indexed()` false.
std::vector<LocalVariable*> pure_locals(scope_ptr scope);


// Module state query functions:

// Normal (not API, net_receive) procedures in module:

std::vector<ProcedureExpression*> normal_procedures(const Module&);

struct public_variable_ids_t {
    std::vector<Id> state_ids;
    std::vector<Id> global_parameter_ids;
    std::vector<Id> range_parameter_ids;
};

// Public module variables by role.

public_variable_ids_t public_variable_ids(const Module&);

struct module_variables_t {
    std::vector<VariableExpression*> scalars;
    std::vector<VariableExpression*> arrays;
};

// Scalar and array variables with local linkage.

module_variables_t local_module_variables(const Module&);

// "normal" procedures in a module.
// A normal procedure is one that has been declared with the
// PROCEDURE keyword in NMODL.

std::vector<ProcedureExpression*> module_normal_procedures(const Module& m);

// Extract key procedures from module.

APIMethod* find_api_method(const Module& m, const char* which);

NetReceiveExpression* find_net_receive(const Module& m);

struct indexed_variable_info {
    std::string data_var;
    std::string index_var;

    bool accumulate = true; // true => add with weight_ factor on assignment
    bool readonly = false;  // true => can never be assigned to by a mechanism

    // Scale is the conversion factor from the data variable
    // to the NMODL value.
    double scale = 1;
    bool scalar() const { return index_var.empty(); }
};

indexed_variable_info decode_indexed_variable(IndexedVariable* sym);
