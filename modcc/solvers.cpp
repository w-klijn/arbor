#include <map>
#include <set>
#include <string>
#include <vector>

#include "astmanip.hpp"
#include "expression.hpp"
#include "parser.hpp"
#include "solvers.hpp"
#include "symdiff.hpp"
#include "symge.hpp"
#include "visitor.hpp"

// Cnexp solver visitor implementation.

void CnexpSolverVisitor::visit(BlockExpression* e) {
    // Do a first pass to extract variables comprising ODE system
    // lhs; can't really trust 'STATE' block.

    for (auto& stmt: e->statements()) {
        if (stmt && stmt->is_assignment() && stmt->is_assignment()->lhs()->is_derivative()) {
            auto id = stmt->is_assignment()->lhs()->is_derivative();
            dvars_.push_back(id->name());
        }
    }

    BlockRewriterBase::visit(e);
}

void CnexpSolverVisitor::visit(AssignmentExpression *e) {
    auto loc = e->location();
    scope_ptr scope = e->scope();

    auto lhs = e->lhs();
    auto rhs = e->rhs();
    auto deriv = lhs->is_derivative();

    if (!deriv) {
        statements_.push_back(e->clone());
        return;
    }

    auto s = deriv->name();
    linear_test_result r = linear_test(rhs, dvars_);

    if (!r.monolinear(s)) {
        error({"System not diagonal linear for cnexp", loc});
        return;
    }

    Expression* coef = r.coef[s].get();
    if (!coef || is_zero(coef)) {
        // s' = b becomes s = s + b*dt; use b_ as a local variable for
        // the constant term b.

        auto local_b_term = make_unique_local_assign(scope, r.constant.get(), "b_");
        statements_.push_back(std::move(local_b_term.local_decl));
        statements_.push_back(std::move(local_b_term.assignment));
        auto b_ = local_b_term.id->is_identifier()->spelling();

        std::string s_update = pprintf("% = %+%*dt", s, s, b_);
        statements_.push_back(Parser{s_update}.parse_line_expression());
        return;
    }
    else if (r.is_homogeneous) {
        // s' = a*s becomes s = s*exp(a*dt); use a_ as a local variable
        // for the coefficient.
        auto local_a_term = make_unique_local_assign(scope, coef, "a_");
        statements_.push_back(std::move(local_a_term.local_decl));
        statements_.push_back(std::move(local_a_term.assignment));
        auto a_ = local_a_term.id->is_identifier()->spelling();

        std::string s_update = pprintf("% = %*exp_pade_11(%*dt)", s, s, a_);
        statements_.push_back(Parser{s_update}.parse_line_expression());
        return;
    }
    else {
        // s' = a*s + b becomes s = -b/a + (s+b/a)*exp(a*dt); use
        // a_ as a local variable for the coefficient and ba_ for the
        // quotient.
        //
        // Note though this will be numerically bad for very small
        // (or zero) a. Perhaps re-implement as:
        //     s = s + exprel(a*dt)*(s*a+b)*dt
        // where exprel(x) = (exp(x)-1)/x and can be well approximated
        // by e.g. a Taylor expansion for small x.
        //
        // Special case ('gating variable') when s' = (b-s)/a; rather
        // than implement more general algebraic simplification, jump
        // straight to simplified update: s = b + (s-b)*exp(-dt/a).

        // Check for 'gating' case:
        if (rhs->is_binary() && rhs->is_binary()->op()==tok::divide) {
            auto denom = rhs->is_binary()->rhs();
            if (involves_identifier(denom, s)) {
                goto not_gating;
            }
            auto numer = rhs->is_binary()->lhs();
            linear_test_result r = linear_test(numer, {s});
            if (expr_value(r.coef[s]) != -1) {
                goto not_gating;
            }

            auto local_a_term = make_unique_local_assign(scope, denom, "a_");
            auto a_ = local_a_term.id->is_identifier()->spelling();
            auto local_b_term = make_unique_local_assign(scope, r.constant, "b_");
            auto b_ = local_b_term.id->is_identifier()->spelling();

            statements_.push_back(std::move(local_a_term.local_decl));
            statements_.push_back(std::move(local_a_term.assignment));
            statements_.push_back(std::move(local_b_term.local_decl));
            statements_.push_back(std::move(local_b_term.assignment));

            std::string s_update = pprintf("% = %+(%-%)*exp_pade_11(-dt/%)", s, b_, s, b_, a_);
            statements_.push_back(Parser{s_update}.parse_line_expression());
            return;
        }

not_gating:
        auto local_a_term = make_unique_local_assign(scope, coef, "a_");
        auto a_ = local_a_term.id->is_identifier()->spelling();

        auto ba_expr = make_expression<DivBinaryExpression>(loc,
                           r.constant->clone(), local_a_term.id->clone());
        auto local_ba_term = make_unique_local_assign(scope, ba_expr, "ba_");
        auto ba_ = local_ba_term.id->is_identifier()->spelling();

        statements_.push_back(std::move(local_a_term.local_decl));
        statements_.push_back(std::move(local_a_term.assignment));
        statements_.push_back(std::move(local_ba_term.local_decl));
        statements_.push_back(std::move(local_ba_term.assignment));

        std::string s_update = pprintf("% = -%+(%+%)*exp_pade_11(%*dt)", s, ba_, s, ba_, a_);
        statements_.push_back(Parser{s_update}.parse_line_expression());
        return;
    }
}


// Sparse solver visitor implementation.

static expression_ptr as_expression(symge::symbol_term term) {
    Location loc;
    if (term.is_zero()) {
        return make_expression<IntegerExpression>(loc, 0);
    }
    else {
        return make_expression<MulBinaryExpression>(loc,
            make_expression<IdentifierExpression>(loc, name(term.left)),
            make_expression<IdentifierExpression>(loc, name(term.right)));
    }
}

static expression_ptr as_expression(symge::symbol_term_diff diff) {
    Location loc;
    if (diff.left.is_zero() && diff.right.is_zero()) {
        return make_expression<IntegerExpression>(loc, 0);
    }
    else if (diff.right.is_zero()) {
        return as_expression(diff.left);
    }
    else if (diff.left.is_zero()) {
        return make_expression<NegUnaryExpression>(loc,
            as_expression(diff.right));
    }
    else {
        return make_expression<SubBinaryExpression>(loc,
            as_expression(diff.left),
            as_expression(diff.right));
    }
}

void SparseSolverVisitor::visit(BlockExpression* e) {
    // Do a first pass to extract variables comprising ODE system
    // lhs; can't really trust 'STATE' block.

    for (auto& stmt: e->statements()) {
        if (stmt && stmt->is_assignment() && stmt->is_assignment()->lhs()->is_derivative()) {
            auto id = stmt->is_assignment()->lhs()->is_derivative();
            dvars_.push_back(id->name());
        }
    }
    scale_factor_.resize(dvars_.size());

    BlockRewriterBase::visit(e);
}

void SparseSolverVisitor::visit(CompartmentExpression *e) {
    auto loc = e->location();

    for (auto& s: e->is_compartment()->state_vars()) {
        auto it = std::find(dvars_.begin(), dvars_.end(), s->is_identifier()->spelling());
        if (it == dvars_.end()) {
            error({"COMPARTMENT variable is not used", loc});
            return;
        }
        auto idx = it - dvars_.begin();
        scale_factor_[idx] = e->scale_factor()->clone();
    }
}

void SparseSolverVisitor::visit(AssignmentExpression *e) {
    if (A_.empty()) {
        unsigned n = dvars_.size();
        A_ = symge::sym_matrix(n, n);
    }

    auto loc = e->location();
    scope_ptr scope = e->scope();

    auto lhs = e->lhs();
    auto rhs = e->rhs();
    auto deriv = lhs->is_derivative();

    if (!deriv) {
        statements_.push_back(e->clone());

        auto id = lhs->is_identifier();
        if (id) {
            auto expand = substitute(rhs, local_expr_);
            if (involves_identifier(expand, dvars_)) {
                local_expr_[id->spelling()] = std::move(expand);
            }
        }
        return;
    }

    if (conserve_ && !A_[deq_index_].empty()) {
        deq_index_++;
        return;
    }

    auto s = deriv->name();
    auto expanded_rhs = substitute(rhs, local_expr_);
    linear_test_result r = linear_test(expanded_rhs, dvars_);
    if (!r.is_homogeneous) {
        error({"System not homogeneous linear for sparse", loc});
        return;
    }

    // Populate sparse symbolic matrix for GE.
    if (s!=dvars_[deq_index_]) {
        error({"ICE: inconsistent ordering of derivative assignments", loc});
        return;
    }

    auto dt_expr = make_expression<IdentifierExpression>(loc, "dt");
    auto one_expr = make_expression<NumberExpression>(loc, 1.0);
    for (unsigned j = 0; j<dvars_.size(); ++j) {
        expression_ptr expr;

        // For zero coefficient and diagonal element, the matrix entry is 1.
        // For non-zero coefficient c and diagonal element, the entry is 1-c*dt.
        // Otherwise, for non-zero coefficient c, the entry is -c*dt.

        if (r.coef.count(dvars_[j])) {
            expr = make_expression<MulBinaryExpression>(loc,
                       r.coef[dvars_[j]]->clone(),
                       dt_expr->clone());

            if (scale_factor_[j]) {
                expr =  make_expression<DivBinaryExpression>(loc, std::move(expr), scale_factor_[j]->clone());
            }
        }

        if (j==deq_index_) {
            if (expr) {
                expr = make_expression<SubBinaryExpression>(loc,
                           one_expr->clone(),
                           std::move(expr));
            }
            else {
                expr = one_expr->clone();
            }
        }
        else if (expr) {
                expr = make_expression<NegUnaryExpression>(loc, std::move(expr));
        }

        if (!expr) continue;

        auto local_a_term = make_unique_local_assign(scope, expr.get(), "a_");
        auto a_ = local_a_term.id->is_identifier()->spelling();

        statements_.push_back(std::move(local_a_term.local_decl));
        statements_.push_back(std::move(local_a_term.assignment));

        A_[deq_index_].push_back({j, symtbl_.define(a_)});
    }
    ++deq_index_;
}

void SparseSolverVisitor::visit(ConserveExpression *e) {
    if (A_.empty()) {
        unsigned n = dvars_.size();
        A_ = symge::sym_matrix(n, n);
    }
    conserve_ = true;

    auto loc = e->location();
    scope_ptr scope = e->scope();

    int row_idx;

    // Find a row that contains one of the state variables in the conserve statement
    auto& l = e->lhs()->is_stoich()->terms().front();
    auto ident = l->is_stoich_term()->ident()->is_identifier();
    if (ident) {
        auto it = std::find(dvars_.begin(), dvars_.end(), ident->name());
        if (it!=dvars_.end()) {
            row_idx = it - dvars_.begin();
        } else {
            error({"CONSERVE statement unknown is not a state variable", loc});
            return;
        }
    }
    else {
         error({"ICE: coefficient in state variable is not an identifier", loc});
         return;
    }

    // Replace that row with the conserve statement
    A_[row_idx].clear();

    for (unsigned j = 0; j < dvars_.size(); ++j) {
        auto state = dvars_[j];

        auto& terms = e->lhs()->is_stoich()->terms();
        auto it = std::find_if(terms.begin(), terms.end(), [&state](expression_ptr& p)
            { return p->is_stoich_term()->ident()->is_identifier()->name() == state;});

        if (it != terms.end()) {
            auto expr = (*it)->is_stoich_term()->coeff()->clone();
            if (scale_factor_[j]) {
                expr =  make_expression<MulBinaryExpression>(loc, scale_factor_[j]->clone(), std::move(expr));
            }

            auto local_a_term = make_unique_local_assign(scope, expr.get(), "a_");
            auto a_ = local_a_term.id->is_identifier()->spelling();

            statements_.push_back(std::move(local_a_term.local_decl));
            statements_.push_back(std::move(local_a_term.assignment));

            A_[row_idx].push_back({j, symtbl_.define(a_)});
        }
    }


    expression_ptr expr = e->rhs()->clone();
    auto local_a_term = make_unique_local_assign(scope, expr.get(), "a_");
    auto a_ = local_a_term.id->is_identifier()->spelling();

    statements_.push_back(std::move(local_a_term.local_decl));
    statements_.push_back(std::move(local_a_term.assignment));

    conserve_rhs_.push_back(a_);
    conserve_idx_.push_back(row_idx);

}

void SparseSolverVisitor::finalize() {
    std::vector<symge::symbol> rhs;
    for (const auto& var: dvars_) {
        rhs.push_back(symtbl_.define(var));
    }
    if (conserve_) {
        for (unsigned i = 0; i < conserve_idx_.size(); ++i) {
            rhs[conserve_idx_[i]] = symtbl_.define(conserve_rhs_[i]);
        }
    }
    A_.augment(rhs);

    symge::gj_reduce(A_, symtbl_);

    // Create and assign intermediate variables.
    for (unsigned i = 0; i<symtbl_.size(); ++i) {
        symge::symbol s = symtbl_[i];

        if (primitive(s)) continue;

        auto expr = as_expression(definition(s));
        auto local_t_term = make_unique_local_assign(block_scope_, expr.get(), "t_");
        auto t_ = local_t_term.id->is_identifier()->spelling();
        symtbl_.name(s, t_);

        statements_.push_back(std::move(local_t_term.local_decl));
        statements_.push_back(std::move(local_t_term.assignment));
    }

    // State variable updates given by rhs/diagonal for reduced matrix.
    Location loc;
    for (unsigned i = 0; i<A_.nrow(); ++i) {
        const symge::sym_row& row = A_[i];
        unsigned rhs_col = A_.augcol();
        unsigned lhs_col;
        for (unsigned r = 0; r < A_.nrow(); r++) {
            if (row[r]) {
                lhs_col = r;
                break;
            }
        }

        auto expr =
            make_expression<AssignmentExpression>(loc,
                make_expression<IdentifierExpression>(loc, dvars_[lhs_col]),
                make_expression<DivBinaryExpression>(loc,
                    make_expression<IdentifierExpression>(loc, symge::name(A_[i][rhs_col])),
                    make_expression<IdentifierExpression>(loc, symge::name(A_[i][lhs_col]))));

        statements_.push_back(std::move(expr));
    }

    BlockRewriterBase::finalize();
}

void LinearSolverVisitor::visit(BlockExpression* e) {
    BlockRewriterBase::visit(e);
}

void LinearSolverVisitor::visit(AssignmentExpression *e) {
    statements_.push_back(e->clone());
    return;
}

void LinearSolverVisitor::visit(LinearExpression *e) {
    auto loc = e->location();
    scope_ptr scope = e->scope();

    if (A_.empty()) {
        unsigned n = dvars_.size();
        A_ = symge::sym_matrix(n, n);
    }

    linear_test_result r = linear_test(e->lhs(), dvars_);
    if (!r.is_homogeneous) {
        error({"System not homogeneous linear for sparse", loc});
        return;
    }

    for (unsigned j = 0; j<dvars_.size(); ++j) {
        expression_ptr expr;

        if (r.coef.count(dvars_[j])) {
            expr = r.coef[dvars_[j]]->clone();
        }

        if (!expr) continue;

        auto a_ = expr->is_identifier()->spelling();

        A_[deq_index_].push_back({j, symtbl_.define(a_)});
    }
    rhs_.push_back(symtbl_.define(e->rhs()->is_identifier()->spelling()));
    ++deq_index_;
}
void LinearSolverVisitor::finalize() {
    A_.augment(rhs_);

    symge::gj_reduce(A_, symtbl_);

    // Create and assign intermediate variables.
    for (unsigned i = 0; i<symtbl_.size(); ++i) {
        symge::symbol s = symtbl_[i];

        if (primitive(s)) continue;

        auto expr = as_expression(definition(s));
        auto local_t_term = make_unique_local_assign(block_scope_, expr.get(), "t_");
        auto t_ = local_t_term.id->is_identifier()->spelling();
        symtbl_.name(s, t_);

        statements_.push_back(std::move(local_t_term.local_decl));
        statements_.push_back(std::move(local_t_term.assignment));
    }

    // State variable updates given by rhs/diagonal for reduced matrix.
    Location loc;
    for (unsigned i = 0; i < A_.nrow(); ++i) {
        const symge::sym_row& row = A_[i];
        unsigned rhs = A_.augcol();
        unsigned lhs;
        for (unsigned r = 0; r < A_.nrow(); r++) {
            if (row[r]) {
                lhs = r;
                break;
            }
        }

        auto expr =
            make_expression<AssignmentExpression>(loc,
                    make_expression<IdentifierExpression>(loc, dvars_[lhs]),
                    make_expression<DivBinaryExpression>(loc,
                            make_expression<IdentifierExpression>(loc, symge::name(A_[i][rhs])),
                            make_expression<IdentifierExpression>(loc, symge::name(A_[i][lhs]))));

        statements_.push_back(std::move(expr));
    }

    BlockRewriterBase::finalize();
}

// Implementation for `remove_unused_locals`: uses two visitors,
// `UnusedVisitor` and `RemoveVariableVisitor` below.

class UnusedVisitor : public Visitor {
protected:
    std::multimap<std::string, std::string> deps;
    std::set<std::string> unused_ids;
    std::set<std::string> used_ids;
    Symbol* lhs_sym = nullptr;

public:
    using Visitor::visit;

    UnusedVisitor() {}

    virtual void visit(Expression* e) override {}

    virtual void visit(BlockExpression* e) override {
        for (auto& s: e->statements()) {
            s->accept(this);
        }
    }

    virtual void visit(AssignmentExpression* e) override {
        auto lhs = e->lhs()->is_identifier();
        if (!lhs) return;

        lhs_sym = lhs->symbol();
        e->rhs()->accept(this);
        lhs_sym = nullptr;
    }

    virtual void visit(UnaryExpression* e) override {
        e->expression()->accept(this);
    }

    virtual void visit(BinaryExpression* e) override {
        e->lhs()->accept(this);
        e->rhs()->accept(this);
    }

    virtual void visit(CallExpression* e) override {
        for (auto& a: e->args()) {
            a->accept(this);
        }
    }

    virtual void visit(IfExpression* e) override {
        e->condition()->accept(this);
        e->true_branch()->accept(this);
        e->false_branch()->accept(this);
    }

    virtual void visit(IdentifierExpression* e) override {
        if (lhs_sym && lhs_sym->is_local_variable()) {
            deps.insert({lhs_sym->name(), e->name()});
        }
        else {
            used_ids.insert(e->name());
        }
    }

    virtual void visit(LocalDeclaration* e) override {
        for (auto& v: e->variables()) {
            unused_ids.insert(v.first);
        }
    }

    std::set<std::string> unused_locals() {
        if (!computed_) {
            for (auto& id: used_ids) {
                remove_deps_from_unused(id);
            }
            computed_ = true;
        }
        return unused_ids;
    }

    void reset() {
        deps.clear();
        unused_ids.clear();
        used_ids.clear();
        computed_ = false;
    }

private:
    bool computed_ = false;

    void remove_deps_from_unused(const std::string& id) {
        auto range = deps.equal_range(id);
        for (auto i = range.first; i != range.second; ++i) {
            if (unused_ids.count(i->second)) {
                remove_deps_from_unused(i->second);
            }
        }
        unused_ids.erase(id);
    }
};

class RemoveVariableVisitor: public BlockRewriterBase {
    std::set<std::string> remove_;

public:
    using BlockRewriterBase::visit;

    RemoveVariableVisitor(std::set<std::string> ids):
        remove_(std::move(ids)) {}

    RemoveVariableVisitor(std::set<std::string> ids, scope_ptr enclosing):
        BlockRewriterBase(enclosing), remove_(std::move(ids)) {}

    virtual void visit(LocalDeclaration* e) override {
        auto replacement = e->clone();
        auto& vars = replacement->is_local_declaration()->variables();

        for (const auto& id: remove_) {
            vars.erase(id);
        }
        if (!vars.empty()) {
            statements_.push_back(std::move(replacement));
        }
    }

    virtual void visit(AssignmentExpression* e) override {
        std::string lhs_id = e->lhs()->is_identifier()->name();
        if (!remove_.count(lhs_id)) {
            statements_.push_back(e->clone());
        }
    }
};

expression_ptr remove_unused_locals(BlockExpression* block) {
    UnusedVisitor unused_visitor;
    block->accept(&unused_visitor);

    RemoveVariableVisitor remove_visitor(unused_visitor.unused_locals());
    block->accept(&remove_visitor);
    return remove_visitor.as_block(false);
}
