/*
 * Decode-conditioned cofactoring for Yosys mux data ports.
 *
 * Copyright (c) 2026 PMUX competition project contributors
 * SPDX-License-Identifier: ISC
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
 * OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
 * CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */

#include "kernel/yosys.h"
#include "kernel/sigtools.h"
#include <algorithm>
#include <functional>
#include <map>
#include <set>
#include <tuple>
#include <utility>
#include <vector>

USING_YOSYS_NAMESPACE
PRIVATE_NAMESPACE_BEGIN

namespace PmuxSmallRom {
using namespace Yosys;

struct Statistics {
	int roms = 0, data_bits = 0, luts = 0, muxes = 0, functional_outputs = 0;
};

struct Table {
	// Each entry identifies a bit in the original little-endian address bus.
	std::vector<int> inputs, values;
};
using Key = std::pair<std::vector<int>, std::vector<int>>;

inline Table shrink(Table table)
{
	for (int bit = int(table.inputs.size()) - 1; bit >= 0; --bit) {
		std::vector<int> low, high;
		for (int code = 0; code < int(table.values.size()); ++code)
			((code >> bit) & 1 ? high : low).push_back(table.values[code]);
		if (low == high) {
			table.inputs.erase(table.inputs.begin() + bit);
			table.values = std::move(low);
		}
	}
	return table;
}

inline std::vector<int> complement_class(const std::vector<int> &values)
{
	std::vector<int> inverted;
	for (int value : values)
		inverted.push_back(1 - value);
	return std::min(values, inverted);
}

inline Table cofactor(const Table &table, int position, int value)
{
	Table result;
	result.inputs = table.inputs;
	result.inputs.erase(result.inputs.begin() + position);
	for (int code = 0; code < int(table.values.size()); ++code)
		if (((code >> position) & 1) == value)
			result.values.push_back(table.values[code]);
	return shrink(std::move(result));
}

struct Choice {
	int split = 0, cofactor_classes = 0, mux_classes = 0;
	int cost() const { return cofactor_classes + mux_classes; }
};

inline Choice choose_split(const std::vector<Table> &outputs, int abits)
{
	Choice best;
	best.cofactor_classes = 1000000;
	for (int split = 0; split < abits; ++split) {
		std::set<Key> lut_classes;
		std::set<std::vector<int>> mux_classes;
		auto account = [&](const Table &table) {
			// Constants, identity, and inversion are free in this estimate.
			// Emission still uses an exact LUT for an inverted single input.
			if (table.inputs.size() > 1)
				lut_classes.emplace(table.inputs, complement_class(table.values));
		};
		for (const Table &original : outputs) {
			Table whole = shrink(original);
			if (whole.inputs.size() <= 4) {
				account(whole);
			} else {
				// At most five original address bits; a five-input residual
				// therefore retains the original bit order and split index.
				mux_classes.insert(complement_class(original.values));
				account(cofactor(original, split, 0));
				account(cofactor(original, split, 1));
			}
		}
		Choice current;
		current.split = split;
		current.cofactor_classes = int(lut_classes.size());
		current.mux_classes = int(mux_classes.size());
		// Ascending iteration and strict comparison implement low-bit ties.
		if (current.cost() < best.cost())
			best = current;
	}
	return best;
}

inline bool exact_integer(Cell *cell, IdString name, int expected)
{
	if (!cell->hasParam(name))
		return false;
	const Const &value = cell->getParam(name);
	return GetSize(value) && value.is_fully_def() && value.convertible_to_int() &&
		value.as_int() == expected;
}

inline bool known_integer(Cell *cell, IdString name)
{
	if (!cell->hasParam(name))
		return false;
	const Const &value = cell->getParam(name);
	return GetSize(value) && value.is_fully_def() && value.convertible_to_int();
}

inline bool exact_ports(Cell *cell, std::initializer_list<IdString> names)
{
	if (cell->connections().size() != names.size())
		return false;
	for (IdString name : names)
		if (!cell->hasPort(name))
			return false;
	return true;
}

struct Candidate {
	RTLIL::Memory *memory = nullptr;
	Cell *read = nullptr, *init = nullptr;
	SigSpec address, original_data;
	std::vector<Table> outputs;
	Choice choice;
};

inline bool inspect(Candidate &candidate, const std::vector<Cell *> &users, SigMap &sigmap)
{
	auto memory = candidate.memory;
	if (memory->get_bool_attribute(ID::keep) || users.size() != 2)
		return false;
	for (Cell *cell : users) {
		if (cell->get_bool_attribute(ID::keep))
			return false;
		if (cell->type == ID($memrd_v2) && !candidate.read)
			candidate.read = cell;
		else if (cell->type.in(ID($meminit), ID($meminit_v2)) && !candidate.init)
			candidate.init = cell;
		else
			return false; // Includes every write, additional read/init, or unknown user.
	}
	Cell *rd = candidate.read, *init = candidate.init;
	if (!rd || !init || !exact_ports(rd, {ID::CLK, ID::EN, ID::ARST, ID::SRST, ID::ADDR, ID::DATA}))
		return false;
	bool init_v2 = init->type == ID($meminit_v2);
	if (init_v2 ? !exact_ports(init, {ID::ADDR, ID::DATA, ID::EN}) : !exact_ports(init, {ID::ADDR, ID::DATA}))
		return false;
	int abits = GetSize(rd->getPort(ID::ADDR));
	if (abits < 1 || abits > 5 || memory->start_offset != 0 || memory->size != (1 << abits) ||
			memory->width < 1 || memory->width > 64 ||
			!exact_integer(rd, ID::ABITS, abits) || !exact_integer(rd, ID::WIDTH, memory->width) ||
			!exact_integer(rd, ID::CLK_ENABLE, 0) || !exact_integer(rd, ID::CE_OVER_SRST, 0) ||
			!(exact_integer(rd, ID::CLK_POLARITY, 0) || exact_integer(rd, ID::CLK_POLARITY, 1)) ||
			GetSize(rd->getPort(ID::CLK)) != 1 || GetSize(rd->getPort(ID::DATA)) != memory->width ||
			sigmap(rd->getPort(ID::EN)) != SigSpec(State::S1) ||
			sigmap(rd->getPort(ID::ARST)) != SigSpec(State::S0) ||
			sigmap(rd->getPort(ID::SRST)) != SigSpec(State::S0))
		return false;
	for (IdString parameter : {ID::ARST_VALUE, ID::SRST_VALUE, ID::INIT_VALUE})
		if (!rd->hasParam(parameter) || GetSize(rd->getParam(parameter)) != memory->width)
			return false;
	for (IdString parameter : {ID::TRANSPARENCY_MASK, ID::COLLISION_X_MASK})
		if (!rd->hasParam(parameter) || GetSize(rd->getParam(parameter)) != 0)
			return false;
	int init_abits = GetSize(init->getPort(ID::ADDR));
	if (init_abits < 1 || init_abits > 32 || !exact_integer(init, ID::ABITS, init_abits) ||
			!exact_integer(init, ID::WIDTH, memory->width) || !exact_integer(init, ID::WORDS, memory->size) ||
			!known_integer(init, ID::PRIORITY) || !sigmap(init->getPort(ID::ADDR)).is_fully_zero() ||
			(init_v2 && sigmap(init->getPort(ID::EN)) != SigSpec(State::S1, memory->width)))
		return false;
	SigSpec data = sigmap(init->getPort(ID::DATA));
	if (GetSize(data) != memory->width * memory->size || !data.is_fully_def())
		return false;
	candidate.address = sigmap(rd->getPort(ID::ADDR));
	candidate.original_data = rd->getPort(ID::DATA);
	for (SigBit bit : candidate.address)
		if (!bit.wire && bit.data != State::S0 && bit.data != State::S1)
			return false;
	for (int output = 0; output < memory->width; ++output) {
		Table table;
		for (int bit = 0; bit < abits; ++bit)
			table.inputs.push_back(bit);
		for (int code = 0; code < memory->size; ++code)
			table.values.push_back(data[code * memory->width + output].data == State::S1);
		candidate.outputs.push_back(std::move(table));
	}
	candidate.choice = choose_split(candidate.outputs, abits);
	return true;
}

class Emitter {
	Module *module;
	SigSpec address;
	int split;
	std::string source;
	Statistics &stats;
	std::map<Key, SigBit> cache;

public:
	Emitter(Module *module, const Candidate &candidate, Statistics &stats) :
		module(module), address(candidate.address), split(candidate.choice.split),
		source(candidate.read->get_src_attribute()), stats(stats) { }

	SigBit emit(Table table) {
		table = shrink(std::move(table));
		Key key(table.inputs, table.values);
		auto found = cache.find(key);
		if (found != cache.end())
			return found->second;
		SigBit result;
		if (table.inputs.empty()) {
			result = table.values[0] ? State::S1 : State::S0;
		} else if (table.inputs.size() == 1 && table.values[0] == 0 && table.values[1] == 1) {
			result = address[table.inputs[0]];
		} else if (table.inputs.size() <= 4) {
			SigSpec inputs;
			for (int bit : table.inputs)
				inputs.append(address[bit]);
			Const lut(State::S0, int(table.values.size()));
			for (int code = 0; code < int(table.values.size()); ++code)
				lut.set(code, table.values[code] ? State::S1 : State::S0);
			result = SigBit(module->addWire(NEW_ID));
			module->addLut(NEW_ID, inputs, result, lut, source);
			++stats.luts;
		} else {
			// Five independent inputs are the only non-leaf case.
			int position = int(std::find(table.inputs.begin(), table.inputs.end(), split) - table.inputs.begin());
			log_assert(position < int(table.inputs.size()));
			SigBit low = emit(cofactor(table, position, 0));
			SigBit high = emit(cofactor(table, position, 1));
			result = module->Mux(NEW_ID, low, high, address[split], source)[0];
			++stats.muxes;
		}
		cache.emplace(std::move(key), result);
		return result;
	}

	// Exact Boolean decomposition before ABC9: f(B,S,F)=h(g(B,S),S,F).
	// A five-input function has two exclusive bound inputs B, two shared
	// inputs S, and one free input F. For each S, at most two distinct
	// (f at F=0, f at F=1) cofactors are allowed. This proves one bit g
	// suffices; h has four inputs. All 32 original assignments are retained.
	SigBit emit_functional(Table original) {
		Table table = shrink(original);
		if (table.inputs.size() != 5)
			return emit(std::move(table));
		Key key(table.inputs, table.values);
		auto found = cache.find(key);
		if (found != cache.end())
			return found->second;
		struct Plan { Table bound; std::vector<int> outer, truth; int cost = 100; } best;
		auto leaf_cost = [&](Table t) {
			t = shrink(std::move(t));
			if (cache.count(Key(t.inputs, t.values)) || t.inputs.empty() ||
					(t.inputs.size() == 1 && t.values[0] == 0 && t.values[1] == 1))
				return 0;
			return 1;
		};
		int split_position = int(std::find(table.inputs.begin(), table.inputs.end(), split) - table.inputs.begin());
		log_assert(split_position < 5);
		int fallback_cost = leaf_cost(cofactor(table, split_position, 0)) +
			leaf_cost(cofactor(table, split_position, 1)) + 1;
		for (int free = 0; free < 5; ++free) {
			std::vector<int> bound;
			for (int bit = 0; bit < 5; ++bit)
				if (bit != free) bound.push_back(bit);
			for (int s0 = 0; s0 < 4; ++s0) for (int s1 = s0 + 1; s1 < 4; ++s1) {
				std::vector<std::vector<int>> signatures(4);
				bool valid = true;
				for (int code = 0; code < 16; ++code) {
					int address_code = 0;
					for (int bit = 0; bit < 4; ++bit)
						address_code |= ((code >> bit) & 1) << bound[bit];
					int signature = table.values[address_code] | (table.values[address_code | (1 << free)] << 1);
					int shared = ((code >> s0) & 1) | (((code >> s1) & 1) << 1);
					auto &values = signatures[shared];
					if (std::find(values.begin(), values.end(), signature) == values.end()) values.push_back(signature);
					if (values.size() > 2) valid = false;
				}
				if (!valid) continue;
				for (auto &values : signatures) std::sort(values.begin(), values.end());
				// Complementing g independently in each shared-input cofactor
				// leaves f unchanged, and can expose a previously emitted g.
				for (int flips = 0; flips < 16; ++flips) {
					Plan current;
					for (int bit : bound) current.bound.inputs.push_back(table.inputs[bit]);
					current.outer = {table.inputs[bound[s0]], table.inputs[bound[s1]], table.inputs[free]};
					current.truth.assign(16, 0);
					for (int code = 0; code < 16; ++code) {
						int address_code = 0;
						for (int bit = 0; bit < 4; ++bit) address_code |= ((code >> bit) & 1) << bound[bit];
						int shared = ((code >> s0) & 1) | (((code >> s1) & 1) << 1);
						int signature = table.values[address_code] | (table.values[address_code | (1 << free)] << 1);
						const auto &values = signatures[shared];
						int encoding = int(std::find(values.begin(), values.end(), signature) - values.begin()) ^ ((flips >> shared) & 1);
						current.bound.values.push_back(encoding);
						current.truth[encoding | (shared << 1)] = signature & 1;
						current.truth[encoding | (shared << 1) | 8] = (signature >> 1) & 1;
					}
					current.bound = shrink(std::move(current.bound));
					current.cost = leaf_cost(current.bound) + 1;
					if (current.cost < best.cost) best = std::move(current);
				}
			}
		}
		if (best.cost >= fallback_cost)
			return emit(std::move(table));
		SigSpec inputs(emit(best.bound));
		for (int bit : best.outer) inputs.append(address[bit]);
		Const truth(State::S0, 16);
		for (int code = 0; code < 16; ++code) truth.set(code, best.truth[code] ? State::S1 : State::S0);
		SigBit result(module->addWire(NEW_ID));
		module->addLut(NEW_ID, inputs, result, truth, source);
		++stats.luts;
		++stats.functional_outputs;
		cache.emplace(std::move(key), result);
		return result;
	}

};

inline Statistics decompose(Module *module, bool verbose = false, bool functional = false)
{
	Statistics stats;
	// The usual PMUX-only module does not allocate a SigMap or scan any cells.
	if (module->memories.empty() || !module->design->selected_whole_module(module) ||
			module->get_bool_attribute(ID::keep) || !module->processes.empty())
		return stats;
	dict<IdString, RTLIL::Memory *> possible;
	for (auto item : module->memories) {
		auto mem = item.second;
		if (mem->width >= 1 && mem->width <= 64 && mem->start_offset == 0 && mem->size >= 2 &&
				mem->size <= 32 && !(mem->size & (mem->size - 1)) && !mem->get_bool_attribute(ID::keep))
			possible[item.first] = mem;
	}
	if (possible.empty())
		return stats;
	dict<IdString, std::vector<Cell *>> users;
	for (Cell *cell : module->cells()) {
		if (!cell->hasParam(ID::MEMID)) {
			if (cell->type.in(ID($memrd), ID($memrd_v2), ID($memwr), ID($memwr_v2), ID($meminit), ID($meminit_v2)))
				return stats; // Malformed memory bookkeeping: preserve the module.
			continue;
		}
		const Const &id = cell->getParam(ID::MEMID);
		if (!(id.flags & RTLIL::CONST_FLAG_STRING) || !id.is_fully_def())
			return stats;
		IdString memid(id.decode_string());
		if (possible.count(memid))
			users[memid].push_back(cell);
	}
	SigMap sigmap(module);
	std::vector<Candidate> candidates;
	for (auto item : possible) {
		Candidate candidate;
		candidate.memory = item.second;
		if (users.count(item.first) && inspect(candidate, users.at(item.first), sigmap))
			candidates.push_back(std::move(candidate));
	}
	if (candidates.empty())
		return stats;
	// Preflight all candidates against the original graph before adding aliases.
	dict<SigBit, Cell *> drivers;
	pool<SigBit> ambiguous, inputs;
	for (Wire *wire : module->wires())
		if (wire->port_input)
			for (SigBit bit : sigmap(wire))
				if (bit.wire)
					inputs.insert(bit);
	for (Cell *cell : module->cells()) {
		bool known = cell->known();
		for (auto connection : cell->connections()) {
			if (known && !cell->output(connection.first))
				continue;
			for (SigBit bit : sigmap(connection.second)) {
				if (!bit.wire)
					continue;
				if (!known || drivers.count(bit))
					ambiguous.insert(bit);
				else
					drivers[bit] = cell;
			}
		}
	}
	std::vector<Candidate> accepted;
	for (Candidate &candidate : candidates) {
		pool<SigBit> outputs;
		bool safe = true;
		for (SigBit bit : sigmap(candidate.original_data)) {
			if (!bit.wire || inputs.count(bit) || ambiguous.count(bit) || outputs.count(bit) ||
					!drivers.count(bit) || drivers.at(bit) != candidate.read)
				safe = false;
			outputs.insert(bit);
		}
		for (SigBit bit : candidate.address)
			if (outputs.count(bit) || ambiguous.count(bit) || (inputs.count(bit) && drivers.count(bit)))
				safe = false;
		if (safe)
			accepted.push_back(std::move(candidate));
	}
	std::sort(accepted.begin(), accepted.end(), [](const Candidate &a, const Candidate &b) {
		return a.memory->name < b.memory->name;
	});
	for (const Candidate &candidate : accepted) {
		Emitter emitter(module, candidate, stats);
		SigSpec replacement;
		for (const Table &table : candidate.outputs)
			replacement.append(functional ? emitter.emit_functional(table) : emitter.emit(table));
		if (verbose)
			log("  Decomposed small constant ROM %s.%s: %d address bits, %d data bits, split %d, estimated classes %d + %d.\n",
				log_id(module), log_id(candidate.memory), GetSize(candidate.address), candidate.memory->width,
				candidate.choice.split, candidate.choice.cofactor_classes, candidate.choice.mux_classes);
		module->connect(candidate.original_data, replacement);
		++stats.roms;
		stats.data_bits += candidate.memory->width;
		module->remove(candidate.read);
		module->remove(candidate.init);
		module->memories.erase(candidate.memory->name);
		delete candidate.memory;
	}
	return stats;
}

} // namespace PmuxSmallRom


// Facts are local to ONE data-port activation. They are never written to the
// module's global SigMap: a decoded signal remains variable outside the branch.
struct PmuxFacts
{
	dict<SigBit, State> values;
	bool consistent = true;
};

struct PmuxOptWorker
{
	Module *module;
	SigMap sigmap;
	dict<SigBit, Cell *> drivers;
	pool<SigBit> ambiguous;
	bool verbose;
	int changed_cells = 0;
	int replaced_bits = 0;
	int decoded_conditions = 0;
	int shared_muxes = 0;
	int shared_data_choices = 0;
	int priority_chains = 0;
	int priority_branches = 0;
	int projected_muxes = 0;
	int projected_data_bits = 0;
	int default_pruned_muxes = 0;
	int default_pruned_branches = 0;
	int rebuilt_muxes = 0, rebuilt_data_bits = 0, rebuilt_nodes = 0;
	int rebuilt_decoders = 0, rebuild_estimated_gain = 0;
	int rebuild_control_skipped = 0, rebuild_nonleaf_skipped = 0;
	int factored_muxes = 0, factored_operators = 0, factor_estimated_gain = 0;
	int flattened_roots = 0, flattened_muxes = 0, flattened_leaves = 0;
	int flattened_binary_roots = 0, flattened_binary_nodes = 0;
	int flattened_conditions = 0, flatten_estimated_gain = 0;
	int shared_arith_muxes = 0, shared_arith_operators = 0, arith_estimated_gain = 0;
	int shared_arith_removed_cells = 0, shared_arith_distinct_expressions = 0;

	PmuxOptWorker(Module *module, bool verbose) :
		module(module), sigmap(module), verbose(verbose)
	{ }

	// Only these cell types can add facts beyond the selector itself. Keep
	// this conservative superset in sync with learn() when adding inference.
	static bool can_infer_inputs(Cell *cell)
	{
		return cell->type.in(ID($eq), ID($ne), ID($logic_not), ID($reduce_or),
			ID($reduce_bool), ID($reduce_and), ID($not), ID($_NOT_),
			ID($and), ID($logic_and), ID($or), ID($logic_or));
	}

	bool may_cofactor()
	{
		pool<SigBit> selectors;
		for (auto cell : module->selected_cells()) {
			if (!cell->type.in(ID($mux), ID($pmux)) || cell->get_bool_attribute(ID::keep))
				continue;
			auto &a = cell->getPort(ID::A);
			auto &b = cell->getPort(ID::B);
			auto &s = cell->getPort(ID::S);
			int width = GetSize(a);
			if (!width || GetSize(b) != width * GetSize(s))
				continue;
			pool<SigBit> local_selectors;
			for (auto bit : sigmap(s))
				if (bit.wire)
					local_selectors.insert(bit);
			if (local_selectors.empty())
				continue;
			bool variable_data = false;
			for (auto port : {ID::A, ID::B})
				for (auto bit : sigmap(cell->getPort(port))) {
					if (!bit.wire)
						continue;
					// learn() always records its root selector, even when that
					// signal has an unknown or multiple driver. Preserve that case.
					if (local_selectors.count(bit))
						return true;
					variable_data = true;
				}
			if (variable_data)
				for (auto bit : local_selectors)
					selectors.insert(bit);
		}
		if (selectors.empty())
			return false;
		for (auto cell : module->cells()) {
			if (!can_infer_inputs(cell) || !cell->hasPort(ID::Y))
				continue;
			auto &y = cell->getPort(ID::Y);
			// learn() only traverses a driver whose first Y bit is the
			// canonical selector. A possible driver suffices here: ambiguity
			// and cell widths are checked later by the ordinary analysis.
			if (GetSize(y) && selectors.count(sigmap(y[0])))
				return true;
		}
		// No inference can leave a selector, and no selector occurs in its
		// mux's data. Thus no data substitution is possible in this module.
		return false;
	}

	void index_drivers()
	{
		for (auto cell : module->cells()) {
			// Native queries reuse Yosys's initialized builtin metadata and
			// consult module port flags for user cells, avoiding a private
			// CellTypes table and a scan of every library module on each pass.
			bool known = cell->known();
			for (auto &conn : cell->connections()) {
				// Unknown cell ports may be outputs; do not infer facts through a
				// signal also attached to an unknown cell.
				if (!known) {
					for (auto bit : sigmap(conn.second))
						if (bit.wire)
							ambiguous.insert(bit);
					continue;
				}
				if (!cell->output(conn.first))
					continue;
				for (auto bit : sigmap(conn.second)) {
					if (!bit.wire)
						continue;
					if (drivers.count(bit))
						ambiguous.insert(bit);
					else
						drivers[bit] = cell;
				}
			}
		}
	}

	void learn(SigBit bit, State value, PmuxFacts &facts, int depth = 0)
	{
		bit = sigmap(bit);
		if (!facts.consistent)
			return;
		if (!bit.wire) {
			// A constant X/Z is not evidence for either Boolean value.
			if (bit.data != value)
				facts.consistent = false;
			return;
		}
		if (facts.values.count(bit)) {
			if (facts.values.at(bit) != value)
				facts.consistent = false;
			return;
		}
		facts.values[bit] = value;
		// Bounded implication traversal is linear for ordinary decoders and
		// also prevents runaway recursion on malformed combinational cycles.
		if (depth >= 8 || ambiguous.count(bit) || !drivers.count(bit))
			return;
		auto cell = drivers.at(bit);
		if (!cell->hasPort(ID::Y))
			return;
		auto y = sigmap(cell->getPort(ID::Y));
		if (GetSize(y) == 0 || y[0] != bit)
			return;

		if (cell->type.in(ID($eq), ID($ne))) {
			bool equal = (value == State::S1) == (cell->type == ID($eq));
			if (!equal)
				return;
			SigSpec a = sigmap(cell->getPort(ID::A));
			SigSpec b = sigmap(cell->getPort(ID::B));
			int width = std::max(GetSize(a), GetSize(b));
			// Verilog compares signed operands only if BOTH are signed.
			bool is_signed = cell->getParam(ID::A_SIGNED).as_bool() &&
				cell->getParam(ID::B_SIGNED).as_bool();
			a.extend_u0(width, is_signed);
			b.extend_u0(width, is_signed);
			if (a.is_fully_def())
				std::swap(a, b);
			if (!b.is_fully_def())
				return;
			decoded_conditions++;
			for (int i = 0; i < width; i++)
				learn(a[i], b[i].data, facts, depth + 1);
			return;
		}

		if (cell->type.in(ID($logic_not), ID($reduce_or), ID($reduce_bool), ID($reduce_and))) {
			SigSpec a = sigmap(cell->getPort(ID::A));
			bool input_one = cell->type == ID($logic_not) ? value == State::S0 : value == State::S1;
			bool all_known = cell->type == ID($reduce_and) ? input_one : !input_one;
			if (GetSize(a) == 1 || all_known) {
				decoded_conditions++;
				for (auto input : a)
					learn(input, input_one ? State::S1 : State::S0, facts, depth + 1);
			}
			return;
		}

		if (cell->type.in(ID($not), ID($_NOT_)) && GetSize(y) == 1 && GetSize(cell->getPort(ID::A)) == 1) {
			learn(cell->getPort(ID::A)[0], value == State::S1 ? State::S0 : State::S1, facts, depth + 1);
			return;
		}
		if (cell->type.in(ID($and), ID($logic_and), ID($or), ID($logic_or)) && GetSize(y) == 1) {
			bool is_and = cell->type.in(ID($and), ID($logic_and));
			if ((is_and && value == State::S1) || (!is_and && value == State::S0)) {
				SigSpec a = cell->getPort(ID::A), b = cell->getPort(ID::B);
				if (GetSize(a) == 1 && GetSize(b) == 1) {
					learn(a[0], value, facts, depth + 1);
					learn(b[0], value, facts, depth + 1);
				}
			}
		}
	}

	int replace(SigSpec &data, int offset, int width, const PmuxFacts &facts)
	{
		if (!facts.consistent)
			return 0;
		int count = 0;
		for (int i = offset; i < offset + width; i++) {
			if (!data[i].wire)
				continue;
			SigBit bit = sigmap(data[i]);
			if (bit.wire && facts.values.count(bit)) {
				data[i] = facts.values.at(bit);
				count++;
			}
		}
		return count;
	}

	static bool has_variable(const SigSpec &data, int offset, int width)
	{
		for (int i = offset; i < offset + width; i++)
			if (data[i].wire)
				return true;
		return false;
	}

	void run()
	{
		for (auto cell : module->selected_cells()) {
			if (!cell->type.in(ID($mux), ID($pmux)) || cell->get_bool_attribute(ID::keep))
				continue;
			SigSpec a = cell->getPort(ID::A), b = cell->getPort(ID::B), s = cell->getPort(ID::S);
			int width = GetSize(a);
			if (!width || GetSize(b) != width * GetSize(s))
				continue;
			int count = 0;
			// Each B slice is observed only while its own selector is true.
			// Multiple asserted PMUX selectors are already undefined in RTLIL.
			for (int i = 0; i < GetSize(s); i++) {
				// A constant branch has no data bit that facts could replace.
				if (!has_variable(b, i * width, width))
					continue;
				PmuxFacts facts;
				learn(s[i], State::S1, facts);
				count += replace(b, i * width, width, facts);
			}
			// A is observed only when every selector is false. We infer only
			// individually necessary conditions; no one-hot assumption is made.
			if (has_variable(a, 0, width)) {
				PmuxFacts default_facts;
				for (auto bit : s)
					learn(bit, State::S0, default_facts);
				count += replace(a, 0, width, default_facts);
			}
			if (!count)
				continue;
			cell->setPort(ID::A, a);
			cell->setPort(ID::B, b);
			changed_cells++;
			replaced_bits += count;
			if (verbose)
				log("  Cofactored %d data bits in %s.%s.\n", count, log_id(module), log_id(cell));
		}
	}

	void prune_default_equivalent_branches()
	{
		for (auto cell : module->selected_cells()) {
			if (cell->type != ID($pmux) || cell->get_bool_attribute(ID::keep))
				continue;
			SigSpec a = cell->getPort(ID::A);
			SigSpec b = cell->getPort(ID::B);
			SigSpec s = cell->getPort(ID::S);
			int width = GetSize(a), choices = GetSize(s);
			if (width < 1 || width > 64 || choices < 2 || choices > 8 ||
					GetSize(b) != width * choices || GetSize(cell->getPort(ID::Y)) != width)
				continue;

			SigSpec canonical_a = sigmap(a);
			SigSpec kept_b, kept_s;
			int removed = 0;
			for (int branch = 0; branch < choices; branch++) {
				SigSpec data = b.extract(branch * width, width);
				SigSpec constant_data = sigmap(data);
				bool redundant = false;
				// Only the candidate slice must be constant. Other branches
				// may carry arbitrary data, including symbolic register holds.
				if (constant_data.is_fully_def()) {
					PmuxFacts facts;
					learn(s[branch], State::S1, facts);
					if (facts.consistent) {
						SigSpec conditional_a = canonical_a;
						replace(conditional_a, 0, width, facts);
						redundant = conditional_a.is_fully_def() && conditional_a == constant_data;
					}
				}
				if (redundant) {
					removed++;
				} else {
					kept_b.append(data);
					kept_s.append(s[branch]);
				}
			}
			if (!removed)
				continue;

			// S_i implies A == B_i for every removed branch. Any retained
			// active branch still selects its original B_i; with only removed
			// branches active the default A supplies the same defined value.
			// Keep this cell object and Y, so every cached driver stays valid.
			// No wire alias is added, and zero remaining choices are represented
			// by a legal ordinary mux rather than an empty parallel mux.
			if (GetSize(kept_s) <= 1) {
				cell->type = ID($mux);
				cell->parameters.erase(ID::S_WIDTH);
				if (GetSize(kept_s) == 0) {
					kept_s = State::S0;
					kept_b = a;
				}
			} else {
				cell->setParam(ID::S_WIDTH, GetSize(kept_s));
			}
			cell->setPort(ID::B, kept_b);
			cell->setPort(ID::S, kept_s);
			default_pruned_muxes++;
			default_pruned_branches += removed;
			if (verbose)
				log("  Removed %d default-equivalent branches in %s.%s.\n", removed, log_id(module), log_id(cell));
		}
	}

	void project_constants()
	{
		for (auto cell : module->selected_cells()) {
			if (cell->type != ID($pmux) || cell->get_bool_attribute(ID::keep))
				continue;
			int width = GetSize(cell->getPort(ID::A));
			int choices = GetSize(cell->getPort(ID::S));
			// Bound fact matching and only consider fully specified data.
			// The original default A may be arbitrary symbolic data.
			if (width < 2 || width > 64 || choices < 2 || choices > 8 ||
					GetSize(cell->getPort(ID::B)) != width * choices ||
					GetSize(cell->getPort(ID::Y)) != width)
				continue;
			SigSpec b = sigmap(cell->getPort(ID::B));
			if (!b.is_fully_def())
				continue;
			SigSpec projected = b.extract(0, width);
			std::vector<int> varying;
			for (int bit = 0; bit < width; bit++)
				for (int branch = 1; branch < choices; branch++)
					if (b[branch * width + bit] != projected[bit]) {
						varying.push_back(bit);
						break;
					}
			// Equal data vectors are already handled by normal mux reduction.
			if (varying.empty())
				continue;
			SigSpec s = sigmap(cell->getPort(ID::S));
			std::vector<PmuxFacts> facts(choices);
			bool consistent = true;
			for (int branch = 0; branch < choices; branch++) {
				learn(s[branch], State::S1, facts[branch]);
				if (!facts[branch].consistent) {
					consistent = false;
					break;
				}
			}
			if (!consistent)
				continue;
			pool<SigBit> output_bits;
			for (auto bit : sigmap(cell->getPort(ID::Y)))
				output_bits.insert(bit);
			facts.front().values.sort();
			std::vector<SigBit> common;
			for (auto &entry : facts.front().values) {
				if (ambiguous.count(entry.first) || output_bits.count(entry.first))
					continue;
				bool present = true;
				for (int branch = 1; branch < choices; branch++)
					if (!facts[branch].values.count(entry.first)) {
						present = false;
						break;
					}
				if (present)
					common.push_back(entry.first);
			}
			bool complete = true;
			for (auto bit : varying) {
				bool found = false;
				for (auto candidate : common) {
					bool matches = true;
					for (int branch = 0; branch < choices; branch++)
						if (facts[branch].values.at(candidate) != b[branch * width + bit].data) {
							matches = false;
							break;
						}
					if (matches) {
						projected[bit] = candidate;
						found = true;
						break;
					}
				}
				if (!found) {
					complete = false;
					break;
				}
			}
			if (!complete)
				continue;
			// S_i implies B_i == projected. With no active branch the
			// unchanged A remains selected; no don't-care default is needed.
			IdString reduce_name = NEW_ID;
			SigBit active = module->ReduceOr(reduce_name, s)[0];
			drivers[active] = module->cell(reduce_name);
			cell->type = ID($mux);
			cell->parameters.erase(ID::S_WIDTH);
			cell->setPort(ID::B, projected);
			cell->setPort(ID::S, active);
			// Mutating the original cell keeps every Y driver pointer valid.
			// No module connection is introduced, so SigMap is still current.
			projected_muxes++;
			projected_data_bits += GetSize(varying);
			log("  Projected %d data bits in %s.%s while preserving its default.\n",
				GetSize(varying), log_id(module), log_id(cell));
		}
	}

	struct DecodedMux {
		Cell *cell;
		SigSpec control;
		std::vector<SigSpec> table;
		int width;
		int original_choices;
	};

	// Normalize equality decoders to a sorted set of constrained signal bits.
	// An absent decode value is filled with A below, so sparse cases retain
	// their default behavior. Duplicate/overlapping codes are rejected.
	bool decode(SigBit select, SigSpec &control, int &value, int min_bits = 2, int max_bits = 8)
	{
		select = sigmap(select);
		if (!select.wire || ambiguous.count(select) || !drivers.count(select))
			return false;
		auto driver = drivers.at(select);
		if (!driver->hasPort(ID::Y) || GetSize(driver->getPort(ID::Y)) == 0 || sigmap(driver->getPort(ID::Y))[0] != select)
			return false;
		SigSpec a, b;
		if (driver->type == ID($eq)) {
			a = sigmap(driver->getPort(ID::A));
			b = sigmap(driver->getPort(ID::B));
			int width = std::max(GetSize(a), GetSize(b));
			bool is_signed = driver->getParam(ID::A_SIGNED).as_bool() &&
				driver->getParam(ID::B_SIGNED).as_bool();
			a.extend_u0(width, is_signed);
			b.extend_u0(width, is_signed);
			if (a.is_fully_def())
				std::swap(a, b);
			if (!b.is_fully_def())
				return false;
		} else if (driver->type == ID($logic_not)) {
			a = sigmap(driver->getPort(ID::A));
			b = SigSpec(State::S0, GetSize(a));
		} else {
			return false;
		}
		dict<SigBit, State> constraints;
		for (int i = 0; i < GetSize(a); i++) {
			if (!a[i].wire) {
				if (a[i] != b[i])
					return false;
				continue;
			}
			if (constraints.count(a[i]) && constraints.at(a[i]) != b[i].data)
				return false;
			constraints[a[i]] = b[i].data;
		}
		if (GetSize(constraints) < min_bits || GetSize(constraints) > max_bits)
			return false;
		constraints.sort();
		control = SigSpec();
		value = 0;
		for (auto &entry : constraints) {
			if (entry.second == State::S1)
				value |= 1 << GetSize(control);
			control.append(entry.first);
		}
		return true;
	}

	bool read_decoded_mux(Cell *cell, DecodedMux &decoded)
	{
		if (cell->type != ID($pmux) || cell->get_bool_attribute(ID::keep))
			return false;
		if (GetSize(cell->getPort(ID::A)) < 2 || GetSize(cell->getPort(ID::S)) < 4)
			return false;
		SigSpec a = sigmap(cell->getPort(ID::A));
		SigSpec b = sigmap(cell->getPort(ID::B));
		SigSpec s = sigmap(cell->getPort(ID::S));
		if (GetSize(a) < 2 || GetSize(s) < 4 || GetSize(b) != GetSize(a) * GetSize(s))
			return false;
		decoded.cell = cell;
		decoded.width = GetSize(a);
		decoded.original_choices = GetSize(s);
		pool<int> codes;
		for (int i = 0; i < GetSize(s); i++) {
			SigSpec control;
			int code;
			if (!decode(s[i], control, code))
				return false;
			if (i == 0) {
				decoded.control = control;
				decoded.table.assign(1 << GetSize(control), a);
			} else if (control != decoded.control) {
				return false;
			}
			if (codes.count(code))
				return false;
			codes.insert(code);
			decoded.table[code] = b.extract(i * decoded.width, decoded.width);
		}
		return true;
	}

	static std::vector<SigSpec> cofactor(const DecodedMux &mux, int bit, int value)
	{
		std::vector<SigSpec> table;
		for (int i = 0; i < GetSize(mux.table); i++)
			if (((i >> bit) & 1) == value)
				table.push_back(mux.table[i]);
		return table;
	}

	// Upper-bound style LUT4 estimate for a parallel sum of n selected terms:
	// pair AND terms into one LUT, then reduce those results with LUT4 ORs.
	// This is only an admission heuristic. Final ABC9 mapping is measured.
	static int sop_cost(int n)
	{
		if (n <= 0)
			return 0;
		int pairs = (n + 1) / 2;
		return pairs + (pairs + 1) / 3;
	}

	static int matching_choices(const std::vector<SigSpec> &a, const std::vector<SigSpec> &b)
	{
		int matches = 0;
		for (int i = 0; i < GetSize(a); i++)
			if (a[i] == b[i])
				matches++;
		return matches;
	}

	static int pair_cost(int choices, int matches)
	{
		if (matches < 2)
			return 2 * sop_cost(choices);
		int exceptions = choices - matches;
		// A single exception is one ordinary mux. More exceptions add a
		// default-value mux to the selected-term network.
		int wrapper = exceptions == 0 ? 0 : exceptions == 1 ? 1 : sop_cost(exceptions) + 1;
		return sop_cost(matches) + 2 * wrapper;
	}

	SigSpec make_pmux(const SigSpec &default_data, const std::vector<SigSpec> &data,
			const std::vector<SigBit> &selectors)
	{
		if (data.empty())
			return default_data;
		SigSpec b, s;
		for (int i = 0; i < GetSize(data); i++) {
			b.append(data[i]);
			s.append(selectors[i]);
		}
		if (GetSize(data) == 1)
			return module->Mux(NEW_ID, default_data, b, s);
		return module->Pmux(NEW_ID, default_data, b, s);
	}

	std::pair<SigSpec, SigSpec> emit_pair(const std::vector<SigSpec> &a,
			const std::vector<SigSpec> &b, const std::vector<SigBit> &selectors, int width)
	{
		SigSpec zero(State::S0, width);
		int matches = matching_choices(a, b);
		if (matches < 2)
			return {make_pmux(zero, a, selectors), make_pmux(zero, b, selectors)};
		std::vector<SigSpec> common_data, a_data, b_data;
		std::vector<SigBit> common_s, unique_s;
		for (int i = 0; i < GetSize(a); i++) {
			if (a[i] == b[i]) {
				common_data.push_back(a[i]);
				common_s.push_back(selectors[i]);
			} else {
				a_data.push_back(a[i]);
				b_data.push_back(b[i]);
				unique_s.push_back(selectors[i]);
			}
		}
		SigSpec common = make_pmux(zero, common_data, common_s);
		shared_data_choices += matches;
		return {make_pmux(common, a_data, unique_s), make_pmux(common, b_data, unique_s)};
	}

	struct RebuildMux {
		Cell *cell;
		SigSpec control, output, default_data;
		std::vector<SigSpec> table;
		std::vector<Cell *> decoders;
		int width, choices;
	};

	bool valid_rebuild_decoder(Cell *cell)
	{
		using namespace PmuxSmallRom;
		if (cell->type == ID($eq)) {
			if (!exact_ports(cell, {ID::A, ID::B, ID::Y}) ||
					!exact_integer(cell, ID::B_WIDTH, GetSize(cell->getPort(ID::B))) ||
					GetSize(cell->getPort(ID::B)) < 1 || GetSize(cell->getPort(ID::B)) > 64 ||
					(!exact_integer(cell, ID::B_SIGNED, 0) && !exact_integer(cell, ID::B_SIGNED, 1)))
				return false;
		} else if (cell->type == ID($logic_not)) {
			if (!exact_ports(cell, {ID::A, ID::Y}))
				return false;
		} else {
			return false;
		}
		return exact_integer(cell, ID::A_WIDTH, GetSize(cell->getPort(ID::A))) &&
			GetSize(cell->getPort(ID::A)) >= 1 && GetSize(cell->getPort(ID::A)) <= 64 &&
			(exact_integer(cell, ID::A_SIGNED, 0) || exact_integer(cell, ID::A_SIGNED, 1)) &&
			exact_integer(cell, ID::Y_WIDTH, 1) && GetSize(cell->getPort(ID::Y)) == 1;
	}

	bool read_rebuild_mux(Cell *cell, RebuildMux &mux,
			const pool<SigBit> &module_inputs, const pool<SigBit> &kept_bits)
	{
		using namespace PmuxSmallRom;
		if (cell->type != ID($pmux) || cell->get_bool_attribute(ID::keep) ||
				!exact_ports(cell, {ID::A, ID::B, ID::S, ID::Y}))
			return false;
		int width = GetSize(cell->getPort(ID::A)), choices = GetSize(cell->getPort(ID::S));
		if (width < 1 || width > 64 || choices < 2 || choices > 16 ||
				!exact_integer(cell, ID::WIDTH, width) || !exact_integer(cell, ID::S_WIDTH, choices) ||
				GetSize(cell->getPort(ID::B)) != width * choices || GetSize(cell->getPort(ID::Y)) != width)
			return false;
		mux.cell = cell;
		mux.width = width;
		mux.choices = choices;
		mux.output = sigmap(cell->getPort(ID::Y));
		mux.default_data = sigmap(cell->getPort(ID::A));
		pool<SigBit> outputs;
		for (SigBit bit : mux.output) {
			if (!bit.wire || module_inputs.count(bit) || kept_bits.count(bit) ||
					ambiguous.count(bit) || !drivers.count(bit) || drivers.at(bit) != cell || outputs.count(bit))
				return false;
			outputs.insert(bit);
		}
		SigSpec s = sigmap(cell->getPort(ID::S)), b = sigmap(cell->getPort(ID::B));
		pool<int> codes;
		for (int i = 0; i < choices; ++i) {
			SigBit select = s[i];
			if (!select.wire || module_inputs.count(select) || ambiguous.count(select) || !drivers.count(select))
				return false;
			Cell *decoder = drivers.at(select);
			if (!valid_rebuild_decoder(decoder))
				return false;
			SigSpec control;
			int code;
			if (!decode(select, control, code, 1, 4) || codes.count(code))
				return false;
			if (i == 0) {
				mux.control = control;
				mux.table.assign(1 << GetSize(control), mux.default_data);
			} else if (control != mux.control) {
				return false;
			}
			codes.insert(code);
			mux.table[code] = b.extract(i * width, width);
			mux.decoders.push_back(decoder);
		}
		for (SigBit bit : mux.control)
			if (outputs.count(bit))
				return false;
		// Only completed, reachable table rows matter. An undefined A is
		// permitted when every code has an explicit, fully defined branch.
		for (const SigSpec &row : mux.table)
			for (SigBit bit : row)
				if (outputs.count(bit) || (!bit.wire && bit.data != State::S0 && bit.data != State::S1))
					return false;
		return true;
	}

	// A reduced ordered decision DAG with symbolic signal-bit leaves. The
	// only inferred facts are the structurally proved decoder equalities.
	// Registers and arithmetic producing those leaves are never modified.
	struct RebuildDag {
		struct Node { int control, low, high; };
		SigSpec control;
		std::vector<int> order;
		std::vector<SigBit> leaves;
		std::vector<Node> nodes;
		std::vector<std::vector<int>> roots;
		std::map<SigBit, int> leaf_ids;
		std::map<std::tuple<int, int, int>, int> node_ids;

		int leaf(SigBit bit)
		{
			auto it = leaf_ids.find(bit);
			if (it != leaf_ids.end())
				return it->second;
			leaves.push_back(bit);
			int id = -GetSize(leaves);
			leaf_ids.emplace(bit, id);
			return id;
		}

		int build(const std::vector<SigBit> &column, int depth = 0, int code = 0)
		{
			if (depth == GetSize(order))
				return leaf(column[code]);
			int bit = order[depth];
			int low = build(column, depth + 1, code);
			int high = build(column, depth + 1, code | (1 << bit));
			if (low == high)
				return low;
			if (low < 0 && high < 0 && leaves[-low - 1] == SigBit(State::S0) && leaves[-high - 1] == SigBit(State::S1))
				return leaf(control[bit]);
			auto key = std::make_tuple(bit, low, high);
			auto it = node_ids.find(key);
			if (it != node_ids.end())
				return it->second;
			int id = GetSize(nodes);
			nodes.push_back({bit, low, high});
			node_ids.emplace(key, id);
			return id;
		}
	};

	bool rebuild_control_only(SigBit bit, const pool<SigBit> &control,
			dict<SigBit, bool> &cache, pool<SigBit> &visiting, int depth = 0)
	{
		if (!bit.wire || control.count(bit))
			return true;
		if (cache.count(bit))
			return cache.at(bit);
		if (depth >= 8 || visiting.count(bit) || ambiguous.count(bit) || !drivers.count(bit))
			return false;
		auto cell = drivers.at(bit);
		// Only stateless, deterministic combinational cells can establish
		// a pure control function. Register outputs remain symbolic data.
		if (!cell->type.in(ID($not), ID($pos), ID($neg), ID($and), ID($or), ID($xor), ID($xnor),
				ID($reduce_and), ID($reduce_or), ID($reduce_xor), ID($reduce_xnor), ID($reduce_bool),
				ID($logic_not), ID($logic_and), ID($logic_or), ID($eq), ID($ne), ID($eqx), ID($nex),
				ID($lt), ID($le), ID($gt), ID($ge), ID($add), ID($sub), ID($mul), ID($shl), ID($shr),
				ID($sshl), ID($sshr), ID($shift), ID($shiftx), ID($mux), ID($pmux), ID($lut)))
			return false;
		visiting.insert(bit);
		bool result = true;
		for (const auto &conn : cell->connections())
			if (cell->input(conn.first))
				for (SigBit input : sigmap(conn.second))
					if (!rebuild_control_only(input, control, cache, visiting, depth + 1))
						result = false;
		visiting.erase(bit);
		cache[bit] = result;
		return result;
	}

	bool admit_rebuild_structure(const RebuildMux &mux)
	{
		pool<SigBit> control, visiting;
		for (SigBit bit : mux.control)
			control.insert(bit);
		dict<SigBit, bool> cache;
		for (int bit = 0; bit < mux.width; ++bit) {
			bool pure_control = true;
			for (const SigSpec &row : mux.table)
				if (!rebuild_control_only(row[bit], control, cache, visiting))
					pure_control = false;
			// A Boolean function of <=4 control bits already fits one LUT4.
			// Do not credit the syntactic number of equality decoders as its
			// mapped area, or disturb sharing with unrelated control outputs.
			if (pure_control) {
				rebuild_control_skipped++;
				return false;
			}
		}
		for (const SigSpec &row : mux.table)
			for (SigBit bit : row)
				if (bit.wire && drivers.count(bit) && drivers.at(bit)->type.in(ID($mux), ID($pmux))) {
					// Preserve reconvergent select/data sharing above another mux.
					// Reconstruct only the direct leaf layer of the original graph.
					rebuild_nonleaf_skipped++;
					return false;
				}
		return true;
	}

	void rebuild_decoded()
	{
		pool<Cell *> selected;
		pool<SigBit> module_inputs, kept_bits, externally_observed;
		for (auto cell : module->selected_cells())
			selected.insert(cell);
		for (auto wire : module->wires())
			for (SigBit bit : sigmap(wire)) {
				if (wire->port_input)
					module_inputs.insert(bit);
				if (wire->get_bool_attribute(ID::keep))
					kept_bits.insert(bit);
				if (wire->port_output || wire->get_bool_attribute(ID::keep))
					externally_observed.insert(bit);
			}
		std::vector<RebuildMux> candidates;
		int table_bits = 0;
		for (auto cell : module->selected_cells()) {
			RebuildMux mux;
			if (!read_rebuild_mux(cell, mux, module_inputs, kept_bits))
				continue;
			if (!admit_rebuild_structure(mux))
				continue;
			table_bits += mux.width * GetSize(mux.table);
			candidates.push_back(std::move(mux));
			// Bound the entire invocation before changing any cell. This is
			// independent of input names, case identities and signal values.
			if (GetSize(candidates) > 128 || table_bits > 16384)
				return;
		}
		if (candidates.empty())
			return;
		std::sort(candidates.begin(), candidates.end(), [](const RebuildMux &a, const RebuildMux &b) {
			return a.cell->name < b.cell->name;
		});
		std::map<SigSpec, std::vector<int>> groups;
		for (int i = 0; i < GetSize(candidates); ++i)
			groups[candidates[i].control].push_back(i);

		// Port-sensitive fanout prevents claiming a decoder saving when its
		// result is also data, a module output, a kept wire or an unknown port.
		dict<SigBit, std::vector<std::pair<Cell *, IdString>>> users;
		for (auto cell : module->cells())
			for (const auto &conn : cell->connections())
				if (!cell->known() || cell->input(conn.first))
					for (SigBit bit : sigmap(conn.second))
						if (bit.wire)
							users[bit].emplace_back(cell, conn.first);

		for (const auto &group : groups) {
			pool<Cell *> members, decoder_set;
			for (int index : group.second) {
				const auto &mux = candidates[index];
				members.insert(mux.cell);
				for (Cell *decoder : mux.decoders)
					decoder_set.insert(decoder);
			}
			std::vector<Cell *> removable;
			for (Cell *decoder : decoder_set) {
				SigBit select = sigmap(decoder->getPort(ID::Y))[0];
				if (!selected.count(decoder) || decoder->get_bool_attribute(ID::keep) || externally_observed.count(select))
					continue;
				bool all_group_selects = true;
				for (const auto &user : users[select])
					if (!members.count(user.first) || user.second != ID::S)
						all_group_selects = false;
				if (all_group_selects)
					removable.push_back(decoder);
			}
			std::sort(removable.begin(), removable.end(), [](Cell *a, Cell *b) { return a->name < b->name; });

			// Do not double-credit old identical output columns that opt_merge
			// can already share. Preserve the cheapest old representation.
			std::map<std::vector<SigBit>, int> old_columns;
			for (int index : group.second) {
				const auto &mux = candidates[index];
				for (int bit = 0; bit < mux.width; ++bit) {
					std::vector<SigBit> column;
					for (const SigSpec &row : mux.table)
						column.push_back(row[bit]);
					bool equal = std::all_of(column.begin(), column.end(), [&](SigBit value) { return value == column.front(); });
					int cost = equal ? 0 : std::min(mux.choices,
						sop_cost(mux.choices) + (mux.default_data[bit] == SigBit(State::S0) ? 0 : 1));
					auto it = old_columns.find(column);
					if (it == old_columns.end())
						old_columns.emplace(std::move(column), cost);
					else
						it->second = std::min(it->second, cost);
				}
			}
			int original_cost = GetSize(removable);
			for (const auto &column : old_columns)
				original_cost += column.second;
			if (original_cost == 0)
				continue;

			std::vector<int> order;
			for (int bit = 0; bit < GetSize(group.first); ++bit)
				order.push_back(bit);
			RebuildDag best;
			bool have_best = false;
			do {
				RebuildDag current;
				current.control = group.first;
				current.order = order;
				for (int index : group.second) {
					const auto &mux = candidates[index];
					std::vector<int> roots;
					for (int bit = 0; bit < mux.width; ++bit) {
						std::vector<SigBit> column;
						for (const SigSpec &row : mux.table)
							column.push_back(row[bit]);
						roots.push_back(current.build(column));
					}
					current.roots.push_back(std::move(roots));
				}
				if (!have_best || current.nodes.size() < best.nodes.size()) {
					best = std::move(current);
					have_best = true;
				}
			} while (std::next_permutation(order.begin(), order.end()));
			int gain = original_cost - GetSize(best.nodes);
			if (gain <= 0)
				continue;

			// All analysis uses the original graph. New nodes are emitted in
			// dependency order; complete groups share the same interned DAG.
			std::vector<SigBit> emitted;
			auto value = [&](int id) { return id < 0 ? best.leaves[-id - 1] : emitted[id]; };
			for (const auto &node : best.nodes)
				emitted.push_back(module->Mux(NEW_ID, value(node.low), value(node.high), best.control[node.control])[0]);
			for (int i = 0; i < GetSize(group.second); ++i) {
				auto &mux = candidates[group.second[i]];
				SigSpec output;
				for (int root : best.roots[i])
					output.append(value(root));
				module->connect(mux.cell->getPort(ID::Y), output);
				if (verbose)
					log("  Rebuilt decoded mux %s.%s (%d data bits, %d control bits).\n",
						log_id(module), log_id(mux.cell), mux.width, GetSize(mux.control));
				rebuilt_data_bits += mux.width;
				module->remove(mux.cell);
			}
			for (Cell *decoder : removable)
				module->remove(decoder);
			rebuilt_muxes += GetSize(group.second);
			rebuilt_nodes += GetSize(best.nodes);
			rebuilt_decoders += GetSize(removable);
			rebuild_estimated_gain += gain;
		}
	}

	// Ephemeral proof handoff within ONE pmux_opt invocation. Only the local
	// flatten construction may produce these records after proving its leaf
	// partition. No user attribute or serialized netlist can provide one.
	struct PartitionProof {
		Cell *cell;
		IdString name;
		SigSpec a, b, s, y;
		bool exhaustive_selects;
		PartitionProof(Cell *cell, bool exhaustive_selects = false) : cell(cell), name(cell->name),
			a(cell->getPort(ID::A)), b(cell->getPort(ID::B)), s(cell->getPort(ID::S)), y(cell->getPort(ID::Y)),
			exhaustive_selects(exhaustive_selects) { }
	};
	std::vector<PartitionProof> continuation_proofs;

	bool read_partition_proof(Cell *cell, const PartitionProof &proof, RebuildMux &mux,
			const pool<SigBit> &module_inputs, const pool<SigBit> &kept_bits)
	{
		using namespace PmuxSmallRom;
		if (cell != proof.cell || cell->name != proof.name || cell->type != ID($pmux) || cell->get_bool_attribute(ID::keep) ||
				!exact_ports(cell, {ID::A, ID::B, ID::S, ID::Y}) ||
				cell->getPort(ID::A) != proof.a || cell->getPort(ID::B) != proof.b ||
				cell->getPort(ID::S) != proof.s || cell->getPort(ID::Y) != proof.y)
			return false;
		int width = GetSize(proof.a), choices = GetSize(proof.s);
		if (width < 1 || width > 64 || choices < 1 || choices > 32 || GetSize(proof.y) != width ||
				GetSize(proof.b) != width * choices || !exact_integer(cell, ID::WIDTH, width) ||
				!exact_integer(cell, ID::S_WIDTH, choices))
			return false;
		mux.cell = cell;
		mux.width = width;
		mux.choices = choices;
		mux.output = sigmap(proof.y);
		pool<SigBit> output;
		for (SigBit bit : mux.output) {
			if (!bit.wire || module_inputs.count(bit) || kept_bits.count(bit) || ambiguous.count(bit) ||
					!drivers.count(bit) || drivers.at(bit) != cell || output.count(bit))
				return false;
			output.insert(bit);
		}
		SigSpec data = sigmap(proof.b);
		if (!proof.exhaustive_selects)
			data.append(sigmap(proof.a));
		data.append(sigmap(proof.s));
		for (SigBit bit : data)
			if (output.count(bit) || (!bit.wire && bit.data != State::S0 && bit.data != State::S1))
				return false;
		return true;
	}

	// Boolean selector functions share exact residual truth tables. Small
	// functions become LUTs; larger ones are split by a fixed bounded order.
	// No data value or arithmetic leaf is evaluated by this representation.
	struct FlattenConditions {
		struct Node {
			PmuxSmallRom::Table table;
			int control = -1, low = 0, high = 0;
		};
		std::vector<Node> nodes;
		std::map<PmuxSmallRom::Key, int> cache;
		std::vector<int> order;
		// -1/-2 are zero/one; -3-bit is an identity control leaf.
		int build(PmuxSmallRom::Table table)
		{
			table = PmuxSmallRom::shrink(std::move(table));
			if (table.inputs.empty())
				return table.values[0] ? -2 : -1;
			if (table.inputs.size() == 1 && table.values[0] == 0 && table.values[1] == 1)
				return -3 - table.inputs[0];
			PmuxSmallRom::Key key(table.inputs, table.values);
			auto found = cache.find(key);
			if (found != cache.end())
				return found->second;
			Node node;
			if (table.inputs.size() <= 4) {
				node.table = table;
			} else {
				int position = -1;
				for (int bit : order) {
					auto it = std::find(table.inputs.begin(), table.inputs.end(), bit);
					if (it != table.inputs.end()) {
						position = int(it - table.inputs.begin());
						break;
					}
				}
				log_assert(position >= 0);
				node.control = table.inputs[position];
				node.low = build(PmuxSmallRom::cofactor(table, position, 0));
				node.high = build(PmuxSmallRom::cofactor(table, position, 1));
				if (node.low == node.high)
					return node.low;
			}
			int id = GetSize(nodes);
			nodes.push_back(std::move(node));
			cache.emplace(std::move(key), id);
			return id;
		}
	};


	// Internal proc guards may contain an undefined arm that is unreachable
	// under the parent's decoded path. They may only be traversed inside a
	// PMUX root; final full-control enumeration rejects every reachable X/Z.
	bool read_flatten_node(Cell *cell, RebuildMux &mux,
			const pool<SigBit> &module_inputs, const pool<SigBit> &kept_bits)
	{
		using namespace PmuxSmallRom;
		bool guard = cell->type == ID($mux);
		if ((!guard && cell->type != ID($pmux)) || cell->get_bool_attribute(ID::keep) ||
				!exact_ports(cell, {ID::A, ID::B, ID::S, ID::Y}))
			return false;
		int width = GetSize(cell->getPort(ID::A)), choices = GetSize(cell->getPort(ID::S));
		if (width < 1 || width > 64 || choices < (guard ? 1 : 2) || choices > (guard ? 1 : 16) ||
				!exact_integer(cell, ID::WIDTH, width) || (!guard && !exact_integer(cell, ID::S_WIDTH, choices)) ||
				GetSize(cell->getPort(ID::B)) != width * choices || GetSize(cell->getPort(ID::Y)) != width)
			return false;
		SigSpec a = sigmap(cell->getPort(ID::A)), b = sigmap(cell->getPort(ID::B));
		if (guard) {
			auto all_x = [](const SigSpec &value) {
				for (SigBit bit : value)
					if (bit != SigBit(State::Sx))
						return false;
				return true;
			};
			if (all_x(a) == all_x(b))
				return false;
		}
		mux.cell = cell;
		mux.width = width;
		mux.choices = choices;
		mux.output = sigmap(cell->getPort(ID::Y));
		mux.default_data = a;
		pool<SigBit> outputs;
		for (SigBit bit : mux.output) {
			if (!bit.wire || module_inputs.count(bit) || kept_bits.count(bit) || ambiguous.count(bit) ||
					!drivers.count(bit) || drivers.at(bit) != cell || outputs.count(bit))
				return false;
			outputs.insert(bit);
		}
		SigSpec s = sigmap(cell->getPort(ID::S));
		pool<int> codes;
		for (int branch = 0; branch < choices; ++branch) {
			SigBit select = s[branch];
			if (!select.wire || module_inputs.count(select) || ambiguous.count(select) || !drivers.count(select))
				return false;
			Cell *decoder = drivers.at(select);
			SigSpec local_control;
			int code;
			if (!valid_rebuild_decoder(decoder) || !decode(select, local_control, code, 1, 4) || codes.count(code))
				return false;
			if (branch == 0) {
				mux.control = local_control;
				mux.table.assign(1 << GetSize(local_control), a);
			} else if (local_control != mux.control) {
				return false;
			}
			codes.insert(code);
			mux.table[code] = b.extract(branch * width, width);
			mux.decoders.push_back(decoder);
		}
		for (SigBit bit : mux.control)
			if (outputs.count(bit))
				return false;
		for (const SigSpec &row : mux.table)
			for (SigBit bit : row)
				if (outputs.count(bit))
					return false;
		return true;
	}

	void flatten_decoded(std::vector<PartitionProof> &proofs, bool allow_binary = false)
	{
		pool<SigBit> module_inputs, kept_bits, observed;
		for (auto wire : module->wires())
			for (SigBit bit : sigmap(wire)) {
				if (wire->port_input)
					module_inputs.insert(bit);
				if (wire->get_bool_attribute(ID::keep))
					kept_bits.insert(bit);
				if (wire->port_output || wire->get_bool_attribute(ID::keep))
					observed.insert(bit);
			}
		std::vector<RebuildMux> muxes;
		pool<Cell *> guards;
		for (auto cell : module->selected_cells()) {
			RebuildMux mux;
			if (read_flatten_node(cell, mux, module_inputs, kept_bits)) {
				if (cell->type == ID($mux))
					guards.insert(cell);
				muxes.push_back(std::move(mux));
				if (muxes.size() > 128)
					return;
			}
		}
		if (muxes.size() < 2)
			return;
		std::sort(muxes.begin(), muxes.end(), [](const RebuildMux &a, const RebuildMux &b) {
			return a.cell->name < b.cell->name;
		});
		std::map<SigSpec, int> output_owner;
		dict<SigBit, int> bit_owner;
		dict<Cell *, int> cell_index;
		for (int i = 0; i < GetSize(muxes); ++i) {
			output_owner.emplace(muxes[i].output, i);
			cell_index[muxes[i].cell] = i;
			for (SigBit bit : muxes[i].output)
				bit_owner[bit] = i;
		}
		std::vector<std::vector<int>> children(GetSize(muxes));
		pool<int> has_parent;
		for (int i = 0; i < GetSize(muxes); ++i)
			for (const SigSpec &row : muxes[i].table) {
				auto it = output_owner.find(row);
				int child = it == output_owner.end() ? -1 : it->second;
				children[i].push_back(child);
				if (child >= 0)
					has_parent.insert(child);
			}
		dict<SigBit, std::vector<std::pair<Cell *, IdString>>> users;
		for (auto cell : module->cells())
			for (const auto &conn : cell->connections())
				if (!cell->known() || cell->input(conn.first))
					for (SigBit bit : sigmap(conn.second))
						if (bit.wire)
							users[bit].emplace_back(cell, conn.first);
		using ArithKind = std::tuple<IdString, int, int>;
		using ArithExpression = std::pair<SigSpec, SigSpec>;
		std::map<SigSpec, std::pair<ArithKind, ArithExpression>> arithmetic_leaves;
		if (allow_binary)
			for (Cell *op : module->selected_cells()) {
				if (!op->type.in(ID($add), ID($sub)) || op->get_bool_attribute(ID::keep) ||
						!PmuxSmallRom::exact_ports(op, {ID::A, ID::B, ID::Y}))
					continue;
				SigSpec left = sigmap(op->getPort(ID::A)), right = sigmap(op->getPort(ID::B));
				SigSpec output = sigmap(op->getPort(ID::Y));
				int width = GetSize(output);
				if (width < 1 || width > 64 || GetSize(left) != width || GetSize(right) != width ||
						!PmuxSmallRom::exact_integer(op, ID::A_WIDTH, width) ||
						!PmuxSmallRom::exact_integer(op, ID::B_WIDTH, width) ||
						!PmuxSmallRom::exact_integer(op, ID::Y_WIDTH, width) ||
						(!PmuxSmallRom::exact_integer(op, ID::A_SIGNED, 0) && !PmuxSmallRom::exact_integer(op, ID::A_SIGNED, 1)) ||
						(!PmuxSmallRom::exact_integer(op, ID::B_SIGNED, 0) && !PmuxSmallRom::exact_integer(op, ID::B_SIGNED, 1)) ||
						op->getParam(ID::A_SIGNED).as_bool() != op->getParam(ID::B_SIGNED).as_bool())
					continue;
				bool valid = true;
				for (SigBit bit : output)
					if (!bit.wire || ambiguous.count(bit) || !drivers.count(bit) || drivers.at(bit) != op)
						valid = false;
				if (!valid)
					continue;
				if (op->type == ID($add) && right < left)
					std::swap(left, right);
				arithmetic_leaves.emplace(output, std::make_pair(
					ArithKind(op->type, width, op->getParam(ID::A_SIGNED).as_bool()),
					ArithExpression(left, right)));
			}
		int expanded_bits = 0;
		for (int root = 0; root < GetSize(muxes); ++root) {
			if (has_parent.count(root) || guards.count(muxes[root].cell))
				continue;
			pool<int> cone, active;
			std::vector<int> preorder;
			bool bounded = true;
			std::function<void(int, int)> visit = [&](int index, int depth) {
				if (!bounded || active.count(index) || depth > 8) {
					bounded = false;
					return;
				}
				if (cone.count(index))
					return;
				cone.insert(index);
				active.insert(index);
				preorder.push_back(index);
				if (cone.size() > 32) {
					bounded = false;
					return;
				}
				for (int child : children[index])
					if (child >= 0)
						visit(child, depth + 1);
				active.erase(index);
			};
			visit(root, 0);
			if (!bounded || cone.size() < 2)
				continue;
			std::set<SigBit> control_set;
			for (int index : preorder)
				for (SigBit bit : muxes[index].control)
					control_set.insert(bit);
			if (control_set.empty() || control_set.size() > 11)
				continue;
			SigSpec control;
			dict<SigBit, int> control_index;
			for (SigBit bit : control_set) {
				control_index[bit] = GetSize(control);
				control.append(bit);
			}
			int combinations = 1 << GetSize(control), width = muxes[root].width;
			expanded_bits += combinations * width;
			if (expanded_bits > 65536)
				break;
			std::map<SigSpec, int> leaf_index;
			std::vector<SigSpec> leaves;
			std::vector<int> table;
			bool reachable_undefined = false;
			for (int code = 0; code < combinations; ++code) {
				int current = root;
				for (;;) {
					int local = 0;
					for (int bit = 0; bit < GetSize(muxes[current].control); ++bit)
						if (code & (1 << control_index.at(muxes[current].control[bit])))
							local |= 1 << bit;
					int child = children[current][local];
					if (child >= 0) {
						current = child;
						continue;
					}
					const SigSpec &leaf = muxes[current].table[local];
					for (SigBit bit : leaf)
						if (!bit.wire && bit.data != State::S0 && bit.data != State::S1)
							reachable_undefined = true;
					if (reachable_undefined)
						break;
					auto inserted = leaf_index.emplace(leaf, GetSize(leaves));
					if (inserted.second)
						leaves.push_back(leaf);
					table.push_back(inserted.first->second);
					break;
				}
			}
			if (reachable_undefined || leaves.size() < 2 || leaves.size() > 32)
				continue;
			// An inner node observed elsewhere remains live, including its
			// dependencies. Also retain any output used as an opaque data leaf
			// or as a selector control; the transformation never invents credit.
			pool<int> retained;
			for (int index : preorder) {
				if (index == root)
					continue;
				for (SigBit bit : muxes[index].output) {
					if (observed.count(bit))
						retained.insert(index);
					for (const auto &user : users[bit])
						if (!cell_index.count(user.first) || !cone.count(cell_index.at(user.first)) ||
								(user.second != ID::A && user.second != ID::B))
							retained.insert(index);
				}
			}
			auto retain_bit = [&](SigBit bit) {
				if (bit_owner.count(bit) && cone.count(bit_owner.at(bit)))
					retained.insert(bit_owner.at(bit));
			};
			for (SigBit bit : control)
				retain_bit(bit);
			for (const SigSpec &leaf : leaves)
				for (SigBit bit : leaf)
					retain_bit(bit);
			bool changed = true;
			while (changed) {
				changed = false;
				for (int index : preorder)
					if (retained.count(index))
						for (int child : children[index])
							if (child >= 0 && !retained.count(child)) {
								retained.insert(child);
								changed = true;
							}
			}
			if (retained.count(root) || cone.size() - retained.size() < 2)
				continue;
			pool<Cell *> removable, decoders;
			std::map<std::pair<SigSpec, std::vector<SigBit>>, int> old_columns;
			for (int index : preorder) {
				if (retained.count(index))
					continue;
				const auto &mux = muxes[index];
				removable.insert(mux.cell);
				for (Cell *decoder : mux.decoders)
					decoders.insert(decoder);
				for (int bit = 0; bit < width; ++bit) {
					std::vector<SigBit> column;
					for (const SigSpec &row : mux.table)
						column.push_back(row[bit]);
					bool same = std::all_of(column.begin(), column.end(), [&](SigBit value) { return value == column[0]; });
					int cost = same ? 0 : std::min(mux.choices,
						sop_cost(mux.choices) + (mux.default_data[bit] == SigBit(State::S0) ? 0 : 1));
					auto key = std::make_pair(mux.control, column);
					auto inserted = old_columns.emplace(key, cost);
					if (!inserted.second)
						inserted.first->second = std::min(inserted.first->second, cost);
				}
			}
			int original_cost = 0;
			for (const auto &entry : old_columns)
				original_cost += entry.second;
			for (Cell *decoder : decoders) {
				SigBit output = sigmap(decoder->getPort(ID::Y))[0];
				bool exclusive = !observed.count(output) && !decoder->get_bool_attribute(ID::keep);
				for (const auto &user : users[output])
					if (!removable.count(user.first) || user.second != ID::S)
						exclusive = false;
				if (exclusive)
					++original_cost;
			}
			std::vector<std::vector<int>> orders;
			std::vector<int> natural, graph;
			for (int bit = 0; bit < GetSize(control); ++bit)
				natural.push_back(bit);
			for (int index : preorder)
				for (SigBit bit : muxes[index].control) {
					int id = control_index.at(bit);
					if (std::find(graph.begin(), graph.end(), id) == graph.end())
						graph.push_back(id);
				}
			orders.push_back(natural);
			orders.push_back(graph);
			std::reverse(natural.begin(), natural.end());
			std::reverse(graph.begin(), graph.end());
			orders.push_back(natural);
			orders.push_back(graph);
			std::vector<std::vector<int>> distinct_orders;
			for (const auto &order : orders)
				if (std::find(distinct_orders.begin(), distinct_orders.end(), order) == distinct_orders.end())
					distinct_orders.push_back(order);
			orders.swap(distinct_orders);
			FlattenConditions best;
			std::vector<int> best_roots;
			std::vector<bool> best_needed;
			int best_default = -1, best_cost = original_cost;
			for (const auto &order : orders) {
				// Construct every leaf condition once per order. Trying another
				// default only changes reachability, never recomputes truth tables.
				FlattenConditions candidate;
				candidate.order = order;
				std::vector<int> all_roots;
				for (int leaf = 0; leaf < GetSize(leaves); ++leaf) {
					PmuxSmallRom::Table condition;
					for (int bit = 0; bit < GetSize(control); ++bit)
						condition.inputs.push_back(bit);
					for (int chosen : table)
						condition.values.push_back(chosen == leaf);
					all_roots.push_back(candidate.build(std::move(condition)));
				}
				for (int default_leaf = 0; default_leaf < GetSize(leaves); ++default_leaf) {
					std::vector<int> roots;
					std::vector<bool> needed(candidate.nodes.size(), false);
					int needed_count = 0;
					std::function<void(int)> mark = [&](int id) {
						if (id < 0 || needed[id])
							return;
						needed[id] = true;
						++needed_count;
						const auto &node = candidate.nodes[id];
						if (node.control >= 0) {
							mark(node.low);
							mark(node.high);
						}
					};
					for (int leaf = 0; leaf < GetSize(leaves); ++leaf)
						if (leaf != default_leaf) {
							roots.push_back(all_roots[leaf]);
							mark(all_roots[leaf]);
						}
					int choices = GetSize(leaves) - 1;
					int cost = needed_count + width *
						std::min(choices, sop_cost(choices) + (leaves[default_leaf].is_fully_zero() ? 0 : 1));
					// A nonzero default also needs the shared NOR of the select
					// vector. A LUT4 OR tree of n inputs uses ceil((n-1)/3)
					// nodes. This estimate is used only for the opt-in comparison.
					if (allow_binary && !leaves[default_leaf].is_fully_zero())
						cost += (choices + 1) / 3;
					if (cost < best_cost) {
						// Emission only needs nodes; the search cache and order are
						// temporary and have no role in the selected implementation.
						best.nodes = candidate.nodes;
						best_roots = std::move(roots);
						best_needed = std::move(needed);
						best_default = default_leaf;
						best_cost = cost;
					}
				}
			}
			// Keep operator groups visible to the later resource-sharing
			// rounds. This guard depends only on full-width expressions, never
			// source identifiers, case identity or sampled signal values.
			std::map<ArithKind, std::set<ArithExpression>> arithmetic_groups;
			bool preserve_arithmetic = false;
			if (allow_binary)
				for (const SigSpec &leaf : leaves) {
					auto expression = arithmetic_leaves.find(leaf);
					if (expression != arithmetic_leaves.end() &&
							arithmetic_groups[expression->second.first].insert(expression->second.second).second &&
							arithmetic_groups[expression->second.first].size() >= 3)
						preserve_arithmetic = true;
				}
			RebuildDag best_binary;
			bool use_binary = false;
			if (allow_binary && !preserve_arithmetic)
				for (const auto &order : orders) {
					RebuildDag candidate;
					candidate.control = control;
					candidate.order = order;
					std::vector<int> roots;
					for (int bit = 0; bit < width; ++bit) {
						std::vector<SigBit> column;
						column.reserve(table.size());
						for (int chosen : table)
							column.push_back(leaves[chosen][bit]);
						roots.push_back(candidate.build(column));
					}
					candidate.roots.push_back(std::move(roots));
					if (GetSize(candidate.nodes) < best_cost) {
						best_cost = GetSize(candidate.nodes);
						best_binary = std::move(candidate);
						use_binary = true;
					}
				}
			if ((!use_binary && best_default < 0) || original_cost - best_cost <= std::max(2, width / 2))
				continue;
			if (use_binary) {
				std::vector<SigBit> emitted;
				auto value = [&](int id) { return id < 0 ? best_binary.leaves[-id - 1] : emitted[id]; };
				for (const auto &node : best_binary.nodes)
					emitted.push_back(module->Mux(NEW_ID, value(node.low), value(node.high), control[node.control])[0]);
				SigSpec output;
				for (int node : best_binary.roots[0])
					output.append(value(node));
				module->connect(muxes[root].cell->getPort(ID::Y), output);
				log("  Binary-flattened %d nested decoded muxes at %s.%s into %d leaves and %d mux nodes (estimated gain %d).\n",
					GetSize(removable), log_id(module), log_id(muxes[root].cell), GetSize(leaves), GetSize(emitted), original_cost - best_cost);
				for (Cell *discard : removable)
					module->remove(discard);
				++flattened_roots;
				++flattened_binary_roots;
				flattened_binary_nodes += GetSize(emitted);
				flattened_muxes += GetSize(removable);
				flattened_leaves += GetSize(leaves);
				flatten_estimated_gain += original_cost - best_cost;
				continue;
			}
			std::vector<SigBit> emitted;
			auto value = [&](int id) -> SigBit {
				return id == -1 ? SigBit(State::S0) : id == -2 ? SigBit(State::S1) :
					id < -2 ? control[-id - 3] : emitted[id];
			};
			int emitted_count = 0;
			for (int index = 0; index < GetSize(best.nodes); ++index) {
				if (!best_needed[index]) {
					emitted.push_back(SigBit(State::S0));
					continue;
				}
				++emitted_count;
				const auto &node = best.nodes[index];
				if (node.control >= 0) {
					emitted.push_back(module->Mux(NEW_ID, value(node.low), value(node.high), control[node.control])[0]);
				} else {
					SigSpec inputs;
					for (int bit : node.table.inputs)
						inputs.append(control[bit]);
					Const bits(State::S0, GetSize(node.table.values));
					for (int code = 0; code < GetSize(node.table.values); ++code)
						bits.set(code, node.table.values[code] ? State::S1 : State::S0);
					SigBit output(module->addWire(NEW_ID));
					module->addLut(NEW_ID, inputs, output, bits);
					emitted.push_back(output);
				}
			}
			SigSpec data, select;
			int condition = 0;
			for (int leaf = 0; leaf < GetSize(leaves); ++leaf)
				if (leaf != best_default) {
					data.append(leaves[leaf]);
					select.append(value(best_roots[condition++]));
				}
			// The leaf functions form an exhaustive disjoint partition. Omitting
			// exactly one function makes that leaf the exact default, including
			// every original sparse-code/default path. Raw PMUX selectors never
			// qualify; externally multi-hot input bits retain ordinary values.
			Cell *replacement = module->addPmux(NEW_ID, leaves[best_default], data, select, muxes[root].cell->getPort(ID::Y));
			proofs.emplace_back(replacement, false);
			log("  Flattened %d nested decoded muxes at %s.%s into %d leaves and %d selector nodes (estimated gain %d).\n",
				GetSize(removable), log_id(module), log_id(muxes[root].cell), GetSize(leaves), emitted_count, original_cost - best_cost);
			// Remove exactly the cone nodes proved disposable above. This
			// exposes true arithmetic fanout for the following fresh worker.
			// Every inner node has_parent and is skipped before dereference;
			// nodes shared with another root are in retained, never here.
			for (Cell *discard : removable)
				module->remove(discard);
			++flattened_roots;
			flattened_muxes += GetSize(removable);
			flattened_leaves += GetSize(leaves);
			flattened_conditions += emitted_count;
			flatten_estimated_gain += original_cost - best_cost;
		}
	}

	// Factor identical-width commutative operators on structurally exclusive
	// decoded branches. This is a local algebraic rewrite, not an invocation
	// of opt_share; the symbolic default remains exact in every sparse case.
	void factor_decoded_ops(const std::vector<PartitionProof> &proofs = {}, bool allow_common = true, bool allow_arith = false)
	{
		continuation_proofs.clear();
		pool<Cell *> selected;
		pool<SigBit> module_inputs, kept_bits, observed;
		for (Cell *cell : module->selected_cells())
			selected.insert(cell);
		for (auto wire : module->wires())
			for (SigBit bit : sigmap(wire)) {
				if (wire->port_input)
					module_inputs.insert(bit);
				if (wire->get_bool_attribute(ID::keep))
					kept_bits.insert(bit);
				if (wire->port_output || wire->get_bool_attribute(ID::keep))
					observed.insert(bit);
			}
		dict<SigBit, std::vector<std::pair<Cell *, IdString>>> users;
		for (Cell *cell : module->cells())
			for (const auto &conn : cell->connections())
				if (!cell->known() || cell->input(conn.first))
					for (SigBit bit : sigmap(conn.second))
						if (bit.wire)
							users[bit].emplace_back(cell, conn.first);
		struct Arm { int branch; Cell *op; SigSpec other, second; };
		using Key = std::tuple<IdString, int, SigSpec>;
		struct Plan {
			Cell *mux;
			std::vector<Arm> arms;
			std::vector<Cell *> operators;
			SigSpec common;
			bool complete, dual, exhaustive;
			int gain, width, expressions;
		};
		std::vector<Plan> plans;
		std::vector<Cell *> possible;
		for (Cell *cell : module->selected_cells())
			if (cell->type == ID($pmux) && !cell->get_bool_attribute(ID::keep))
				possible.push_back(cell);
		if (GetSize(possible) > 128)
			return;
		std::sort(possible.begin(), possible.end(), [](Cell *a, Cell *b) { return a->name < b->name; });
		dict<Cell *, const PartitionProof *> proof_by_cell;
		for (const PartitionProof &proof : proofs)
			proof_by_cell[proof.cell] = &proof;
		for (Cell *cell : possible) {
			RebuildMux mux;
			bool codes_complete;
			if (read_rebuild_mux(cell, mux, module_inputs, kept_bits)) {
				codes_complete = mux.choices == (1 << GetSize(mux.control));
			} else {
				if (!proof_by_cell.count(cell))
					continue;
				const PartitionProof &proof = *proof_by_cell.at(cell);
				if (!read_partition_proof(cell, proof, mux, module_inputs, kept_bits))
					continue;
				codes_complete = proof.exhaustive_selects;
			}
			SigSpec a = sigmap(cell->getPort(ID::A)), b = sigmap(cell->getPort(ID::B));
			std::map<Key, std::vector<Arm>> groups;
			for (int branch = codes_complete ? 0 : -1; branch < mux.choices; ++branch) {
				SigSpec data = branch < 0 ? a : b.extract(branch * mux.width, mux.width);
				if (!data[0].wire || ambiguous.count(data[0]) || !drivers.count(data[0]))
					continue;
				Cell *op = drivers.at(data[0]);
				bool common_type = allow_common && op->type.in(ID($xor), ID($xnor), ID($and), ID($or), ID($add));
				bool arith_type = allow_arith && op->type.in(ID($add), ID($sub));
				if ((!common_type && !arith_type) ||
						!selected.count(op) || op->get_bool_attribute(ID::keep) ||
						!PmuxSmallRom::exact_ports(op, {ID::A, ID::B, ID::Y}) ||
						!PmuxSmallRom::exact_integer(op, ID::A_WIDTH, mux.width) ||
						!PmuxSmallRom::exact_integer(op, ID::B_WIDTH, mux.width) ||
						!PmuxSmallRom::exact_integer(op, ID::Y_WIDTH, mux.width) ||
						(!PmuxSmallRom::exact_integer(op, ID::A_SIGNED, 0) && !PmuxSmallRom::exact_integer(op, ID::A_SIGNED, 1)) ||
						(!PmuxSmallRom::exact_integer(op, ID::B_SIGNED, 0) && !PmuxSmallRom::exact_integer(op, ID::B_SIGNED, 1)) ||
						op->getParam(ID::A_SIGNED).as_bool() != op->getParam(ID::B_SIGNED).as_bool() ||
						sigmap(op->getPort(ID::Y)) != data || GetSize(op->getPort(ID::A)) != mux.width ||
						GetSize(op->getPort(ID::B)) != mux.width)
					continue;
				bool valid = true;
				for (SigBit bit : data)
					if (!bit.wire || ambiguous.count(bit) || !drivers.count(bit) || drivers.at(bit) != op ||
							module_inputs.count(bit) || observed.count(bit))
						valid = false;
				SigSpec oa = sigmap(op->getPort(ID::A)), ob = sigmap(op->getPort(ID::B));
				for (const SigSpec &input : {oa, ob})
					for (SigBit bit : input) {
						if (!bit.wire && bit.data != State::S0 && bit.data != State::S1)
							valid = false;
						for (SigBit output : mux.output)
							if (bit == output)
								valid = false;
					}
				if (!valid)
					continue;
				int sign = op->getParam(ID::A_SIGNED).as_bool();
				if (common_type) {
					groups[Key(op->type, sign, oa)].push_back({branch, op, ob, SigSpec()});
					if (oa != ob)
						groups[Key(op->type, sign, ob)].push_back({branch, op, oa, SigSpec()});
				}
				if (arith_type)
					// Empty common operand distinguishes dual-operand resource
					// sharing. SUB keeps A and B in their original positions.
					groups[Key(op->type, sign, SigSpec())].push_back({branch, op, oa, ob});
			}
			Plan best;
			best.gain = 0;
			for (const auto &group : groups) {
				const auto &arms = group.second;
				bool dual = GetSize(std::get<2>(group.first)) == 0;
				if (GetSize(arms) < 2)
					continue;
				pool<int> branches;
				pool<Cell *> operators;
				std::set<std::pair<SigSpec, SigSpec>> expressions;
				for (const Arm &arm : arms) {
					branches.insert(arm.branch);
					operators.insert(arm.op);
					SigSpec left = sigmap(arm.op->getPort(ID::A)), right = sigmap(arm.op->getPort(ID::B));
					if (arm.op->type != ID($sub) && right < left)
						std::swap(left, right);
					// Type, width and signedness are already identical in this
					// group. Do not claim opt_merge's existing duplicate saving.
					expressions.emplace(left, right);
				}
				bool complete = GetSize(branches) == mux.choices + (codes_complete ? 0 : 1);
				bool includes_default = branches.count(-1);
				if (GetSize(expressions) < (dual ? 3 : 2))
					continue;
				bool removable = true;
				for (Cell *op : operators)
					for (SigBit bit : sigmap(op->getPort(ID::Y))) {
						for (const auto &user : users[bit])
							if (user.first != cell || (user.second != ID::B && !((complete || includes_default) && user.second == ID::A)))
								removable = false;
						for (int branch = 0; branch < mux.choices; ++branch)
							if (!branches.count(branch))
								for (SigBit other : b.extract(branch * mux.width, mux.width))
									if (other == bit)
										removable = false;
					}
				if (!removable)
					continue;
				int unit = std::get<0>(group.first).in(ID($add), ID($sub)) ? std::max(1, 2 * mux.width - 1) : mux.width;
				int gain = (GetSize(expressions) - 1) * unit;
				if (complete && dual)
					gain -= mux.width * std::min(mux.choices, sop_cost(mux.choices) + 1);
				if (!complete) {
					if (includes_default) {
						int count = GetSize(branches) - 1;
						int old_select = std::min(mux.choices, sop_cost(mux.choices) + 1);
						int outer_select = std::min(mux.choices - count, sop_cost(mux.choices - count) + 1);
						int inner_select = std::min(count, sop_cost(count) + 1);
						gain += mux.width * (old_select - outer_select - (dual ? 2 : 1) * inner_select);
					} else {
						int count = GetSize(branches);
						gain += mux.width * (sop_cost(mux.choices) - sop_cost(mux.choices - count + 1) - (dual ? 2 : 1) * sop_cost(count)) - 1;
					}
				}
				if (gain <= best.gain)
					continue;
				best.mux = cell;
				best.arms = arms;
				best.common = std::get<2>(group.first);
				best.complete = complete;
				best.dual = dual;
				best.exhaustive = codes_complete;
				best.gain = gain;
				best.width = mux.width;
				best.expressions = GetSize(expressions);
				best.operators.clear();
				for (Cell *op : operators)
					best.operators.push_back(op);
				std::sort(best.operators.begin(), best.operators.end(), [](Cell *x, Cell *y) { return x->name < y->name; });
			}
			if (best.gain > 0)
				plans.push_back(std::move(best));
		}
		for (const Plan &plan : plans) {
			Cell *cell = plan.mux;
			SigSpec original_y = cell->getPort(ID::Y), original_s = cell->getPort(ID::S);
			SigSpec original_b = cell->getPort(ID::B), default_data(State::S0, plan.width);
			SigSpec operand_data, group_selects, second_data, second_default(State::S0, plan.width);
			dict<int, SigSpec> operands, second_operands;
			for (const Arm &arm : plan.arms) {
				operands[arm.branch] = arm.other;
				second_operands[arm.branch] = arm.second;
			}
			SigSpec selected_operand, selected_second;
			if (plan.complete) {
				if (operands.count(-1)) {
					default_data = operands.at(-1);
					if (plan.dual)
						second_default = second_operands.at(-1);
				}
				for (int branch = 0; branch < GetSize(original_s); ++branch) {
					operand_data.append(operands.at(branch));
					if (plan.dual)
						second_data.append(second_operands.at(branch));
				}
				selected_operand = module->addWire(NEW_ID, plan.width);
				cell->setPort(ID::A, default_data);
				cell->setPort(ID::B, operand_data);
				cell->setPort(ID::Y, selected_operand);
				if (plan.dual)
					selected_second = module->Pmux(NEW_ID, second_default, second_data, original_s);
			} else {
				for (const Arm &arm : plan.arms) {
					if (arm.branch < 0) {
						default_data = arm.other;
						if (plan.dual)
							second_default = arm.second;
					} else {
						operand_data.append(arm.other);
						if (plan.dual)
							second_data.append(arm.second);
						group_selects.append(original_s[arm.branch]);
					}
				}
				selected_operand = module->Pmux(NEW_ID, default_data, operand_data, group_selects);
				if (plan.dual)
					selected_second = module->Pmux(NEW_ID, second_default, second_data, group_selects);
			}
			Cell *source = plan.operators.front();
			Cell *factored = module->addCell(NEW_ID, source->type);
			factored->parameters = source->parameters;
			if (source->attributes.count(ID::src))
				factored->attributes[ID::src] = source->attributes.at(ID::src);
			factored->setPort(ID::A, selected_operand);
			factored->setPort(ID::B, plan.dual ? selected_second : plan.common);
			SigSpec combined = plan.complete ? original_y : SigSpec(module->addWire(NEW_ID, plan.width));
			factored->setPort(ID::Y, combined);
			if (!plan.complete) {
				SigSpec data, select;
				for (int branch = 0; branch < GetSize(original_s); ++branch)
					if (!operands.count(branch)) {
						data.append(original_b.extract(branch * plan.width, plan.width));
						select.append(original_s[branch]);
					}
				if (operands.count(-1)) {
					// With no remaining outer selector active, the inner PMUX
					// selects either a grouped B operand or the exact old A operand.
					cell->setPort(ID::A, combined);
				} else {
					data.append(combined);
					select.append(module->ReduceOr(NEW_ID, group_selects));
				}
				cell->setPort(ID::B, data);
				cell->setPort(ID::S, select);
				cell->setParam(ID::S_WIDTH, GetSize(select));
				// This local algebraic rewrite preserves disjoint selectors.
				// Absorbing A removes some explicit coverage, so even a formerly
				// exhaustive S must now allow the reachable outer default.
				continuation_proofs.emplace_back(cell, operands.count(-1) ? false : plan.exhaustive);
			}
			for (Cell *op : plan.operators)
				module->remove(op);
			if (plan.dual) {
				shared_arith_muxes++;
				shared_arith_operators += GetSize(plan.operators) - 1;
				shared_arith_removed_cells += GetSize(plan.operators);
				shared_arith_distinct_expressions += plan.expressions;
				arith_estimated_gain += plan.gain;
			} else {
				factored_muxes++;
				factored_operators += GetSize(plan.operators) - 1;
				factor_estimated_gain += plan.gain;
			}
			if (verbose)
				log("  %s %d operators in %s.%s (%s table, estimated gain %d).\n",
					plan.dual ? "Shared arithmetic resource for" : "Factored common operand of", GetSize(plan.operators),
					log_id(module), log_id(cell), plan.complete ? "complete" : "partial", plan.gain);
		}
		// Retain only unchanged incoming proof snapshots. Complete rewrites
		// change the old cell into an operand mux and are deliberately not
		// continued; no arbitrary new cell gains a partition assumption.
		pool<Cell *> queued;
		for (const PartitionProof &proof : continuation_proofs)
			queued.insert(proof.cell);
		for (const PartitionProof &proof : proofs)
			if (module->cell(proof.name) == proof.cell && selected.count(proof.cell) && !queued.count(proof.cell) && proof.cell->type == ID($pmux) &&
					PmuxSmallRom::exact_ports(proof.cell, {ID::A, ID::B, ID::S, ID::Y}) && proof.cell->getPort(ID::A) == proof.a &&
					proof.cell->getPort(ID::B) == proof.b && proof.cell->getPort(ID::S) == proof.s &&
					proof.cell->getPort(ID::Y) == proof.y)
				continuation_proofs.push_back(proof);
	}

	void share_decoded()
	{
		std::vector<Cell *> possible;
		for (auto cell : module->selected_cells())
			if (cell->type == ID($pmux) && !cell->get_bool_attribute(ID::keep) &&
					GetSize(cell->getPort(ID::A)) >= 2 && GetSize(cell->getPort(ID::S)) >= 4)
				possible.push_back(cell);
		// Sharing requires two outputs. Avoid decoding a lone ROM/PMUX.
		if (GetSize(possible) < 2)
			return;
		std::vector<DecodedMux> candidates;
		for (auto cell : possible) {
			DecodedMux mux;
			if (read_decoded_mux(cell, mux)) {
				candidates.push_back(mux);
				if (GetSize(candidates) > 128)
					return;
			}
		}
		std::sort(candidates.begin(), candidates.end(), [](const DecodedMux &a, const DecodedMux &b) {
			return a.cell->name < b.cell->name;
		});
		pool<Cell *> consumed;
		// Greedy deterministic pairing. Limits above keep each table <= 256
		// entries; at most 128 wide muxes participate in the quadratic search.
		for (int i = 0; i < GetSize(candidates); i++) {
			auto &a = candidates[i];
			if (consumed.count(a.cell))
				continue;
			int best_j = -1, best_bit = -1, best_swap = 0, best_gain = 0;
			for (int j = i + 1; j < GetSize(candidates); j++) {
				auto &b = candidates[j];
				if (consumed.count(b.cell) || a.width != b.width || a.control != b.control)
					continue;
				// Whole-table duplicates are already handled by opt_merge. Do
				// not count their existing merge opportunity as a new saving.
				if (a.table == b.table)
					continue;
				int half = GetSize(a.table) / 2;
				// Include a conservative budget for the newly exposed decoders.
				int decoder_budget = half * ((GetSize(a.control) + 1) / 3);
				int original = (sop_cost(a.original_choices) + sop_cost(b.original_choices)) * a.width;
				for (int bit = 0; bit < GetSize(a.control); bit++) {
					// Count matches in-place. Materialize cofactors only after
					// choosing a profitable pair, not for every search candidate.
					int matches[2][2] = {{0, 0}, {0, 0}};
					int mask = 1 << bit;
					for (int code = 0; code < GetSize(a.table); code++) {
						if (code & mask)
							continue;
						int high = code | mask;
						matches[0][0] += a.table[code] == b.table[code];
						matches[0][1] += a.table[high] == b.table[high];
						matches[1][0] += a.table[code] == b.table[high];
						matches[1][1] += a.table[high] == b.table[code];
					}
					for (int swap = 0; swap < 2; swap++) {
						int m0 = matches[swap][0], m1 = matches[swap][1];
						if (m0 < 2 && m1 < 2)
							continue;
						int cost = (pair_cost(half, m0) + pair_cost(half, m1) + 2) * a.width + decoder_budget;
						int gain = original - cost;
						if (gain > best_gain) {
							best_j = j;
							best_bit = bit;
							best_swap = swap;
							best_gain = gain;
						}
					}
				}
			}
			if (best_j < 0)
				continue;
			auto &b = candidates[best_j];
			SigSpec residual_control = a.control;
			residual_control.remove(best_bit, 1);
			std::vector<SigBit> selectors;
			for (int code = 0; code < (1 << GetSize(residual_control)); code++)
				selectors.push_back(module->Eq(NEW_ID, residual_control, Const(code, GetSize(residual_control)))[0]);
			auto low = emit_pair(cofactor(a, best_bit, 0), cofactor(b, best_bit, best_swap), selectors, a.width);
			auto high = emit_pair(cofactor(a, best_bit, 1), cofactor(b, best_bit, 1 - best_swap), selectors, a.width);
			SigSpec out_a = module->Mux(NEW_ID, low.first, high.first, a.control[best_bit]);
			SigSpec out_b = module->Mux(NEW_ID, best_swap ? high.second : low.second,
				best_swap ? low.second : high.second, a.control[best_bit]);
			module->connect(a.cell->getPort(ID::Y), out_a);
			module->connect(b.cell->getPort(ID::Y), out_b);
			log("  Shared decoded muxes %s.%s and %s using control %s (estimated LUT gain %d).\n",
				log_id(module), log_id(a.cell), log_id(b.cell), log_signal(a.control[best_bit]), best_gain);
			consumed.insert(a.cell);
			consumed.insert(b.cell);
			module->remove(a.cell);
			module->remove(b.cell);
			shared_muxes += 2;
		}
	}

	void pack_priority()
	{
		pool<SigBit> used_as_default;
		pool<Cell *> selected;
		for (auto cell : module->selected_cells())
			if (cell->type == ID($mux) && !cell->get_bool_attribute(ID::keep) &&
					GetSize(cell->getPort(ID::Y)) >= 4)
				selected.insert(cell);
		if (GetSize(selected) < 3)
			return;
		for (auto cell : module->cells())
			if (cell->type == ID($mux))
				for (auto bit : sigmap(cell->getPort(ID::A)))
					used_as_default.insert(bit);
		for (auto root : module->selected_cells()) {
			if (root->type != ID($mux) || root->get_bool_attribute(ID::keep))
				continue;
			SigSpec output = sigmap(root->getPort(ID::Y));
			int width = GetSize(output);
			if (width < 4)
				continue;
			bool is_root = true;
			for (auto bit : output)
				if (used_as_default.count(bit))
					is_root = false;
			if (!is_root)
				continue;
			std::vector<Cell *> chain;
			pool<Cell *> seen;
			SigSpec default_data;
			Cell *current = root;
			while (GetSize(chain) < 8 && current->type == ID($mux) && selected.count(current) &&
					!current->get_bool_attribute(ID::keep) && !seen.count(current)) {
				if (GetSize(current->getPort(ID::A)) != width || GetSize(current->getPort(ID::B)) != width ||
						GetSize(current->getPort(ID::S)) != 1)
					break;
				seen.insert(current);
				chain.push_back(current);
				default_data = sigmap(current->getPort(ID::A));
				if (!default_data[0].wire || ambiguous.count(default_data[0]) || !drivers.count(default_data[0]))
					break;
				auto child = drivers.at(default_data[0]);
				if (!child->hasPort(ID::Y) || sigmap(child->getPort(ID::Y)) != default_data)
					break;
				current = child;
			}
			int n = GetSize(chain);
			if (n < 3)
				continue;
			int selector_cost = 0;
			for (int i = 2; i <= n; i++)
				selector_cost += (i + 1) / 3;
			int cost = width * (sop_cost(n) + (default_data.is_fully_zero() ? 0 : 1)) + selector_cost;
			if (cost >= width * n)
				continue;
			SigSpec selectors, data, prefix;
			for (int i = 0; i < n; i++) {
				prefix.append(sigmap(chain[i]->getPort(ID::S)));
				if (i == 0) {
					selectors.append(prefix);
				} else {
					Const pattern(State::S0, i + 1);
					pattern.set(i, State::S1);
					IdString decoder_name = NEW_ID;
					SigBit select = module->Eq(decoder_name, prefix, pattern)[0];
					selectors.append(select);
					drivers[select] = module->cell(decoder_name);
				}
				data.append(chain[i]->getPort(ID::B));
			}
			// Every selector retains the original priority: S_i AND the
			// negation of ALL earlier S_j. Multi-hot external inputs are valid.
			auto replacement = module->addPmux(NEW_ID, default_data, data, selectors, root->getPort(ID::Y));
			// No module connection/alias is added by priority packing. Fresh
			// decoder wires map to themselves, and the original root output
			// keeps its canonical bits. Update only these changed drivers so
			// the following phases can reuse the existing SigMap and graph.
			for (auto bit : output)
				if (bit.wire && drivers.count(bit) && drivers.at(bit) == root)
					drivers[bit] = replacement;
			log("  Packed %d priority branches in %s.%s into shared mutually exclusive selectors.\n",
				n, log_id(module), log_id(root));
			module->remove(root);
			// Inner chain cells may have other users. Standard opt_clean removes
			// only those made dead by replacement of the root.
			priority_chains++;
			priority_branches += n;
		}
	}
};

struct PmuxOptPass : public Pass
{
	PmuxOptPass() : Pass("pmux_opt", "cofactor mux data under decoded select conditions") { }

	void help() override
	{
		log("\n");
		log("    pmux_opt [-rom-functional] [-verbose] [-share] [-priority] [-prune-default] [-small-rom] [-project] [-rebuild-decode] [-factor-ops] [-share-arith] [-pre-factor] [-flatten-decode] [-flatten-binary] [-no-cofactor] [selection]\n");
		log("\n");
		log("Infer constants implied by $eq/$ne and Boolean select logic, then\n");
		log("substitute those constants in the corresponding $mux/$pmux data ports.\n");
		log("No case names, source paths, input patterns, or target devices are used.\n");
		log("Run after proc, before decoder techmapping and preferably before opt.\n");
		log("Run normal opt afterwards to merge constant/equal slices and registers.\n");
		log("With -share, consider one-level Shannon decomposition of decoded mux\n");
		log("pairs and share equal residual choices across their cofactors. Sparse\n");
		log("decoders preserve A in every missing table entry. A LUT4 cost estimate\n");
		log("rejects unprofitable pairs. -no-cofactor permits isolated ablation tests.\n");
		log("With -priority, pack profitable wide mux chains using shared mutually\n");
		log("exclusive selectors that preserve priority for multi-hot inputs.\n");
		log("With -prune-default, remove a fully defined constant PMUX branch only\n");
		log("when its select condition implies that A equals that branch. Other\n");
		log("symbolic branches stay unchanged. This accepts 1..64 data bits and\n");
		log("2..8 choices, and runs after cofactoring but before projection.\n");
		log("Contradictory facts, X/Z constants and incomplete matches are rejected.\n");
		log("With -small-rom, decompose complete constant ROMs with 1..5 address\n");
		log("bits and 1..64 data bits into LUT4 cofactors and shared muxes. This\n");
		log("requires whole-module selection, one asynchronous always-enabled read,\n");
		log("one full initialization, zero offset, no reset/write and no keep.\n");
		log("A table-derived complement-sharing estimate selects the split bit.\n");
		log("ROM decomposition runs before the mux analysis; its estimate does not\n");
		log("guarantee a mapped area improvement. Disabled unless requested.\n");
		log("With -project, replace small constant-data PMUXs by a mux whose data\n");
		log("bits match the decoder facts in every branch, preserving the original A.\n");
		log("Projection accepts 2..8 choices and 2..64 data bits, and runs after\n");
		log("cofactoring but before sharing. X/Z data and incomplete matches are skipped.\n");
		log("With -rebuild-decode, rebuild PMUXs decoded from one 1..4-bit control\n");
		log("as a reduced Shannon mux DAG with symbolic data and shared subtrees.\n");
		log("Missing codes preserve A; duplicate/overlapping or undefined decoders\n");
		log("and reachable X/Z data are rejected. Accepts 1..64 data bits, 2..16\n");
		log("choices, at most 128 candidates and 16384 expanded table bits per module.\n");
		log("Enumerates at most 24 control orders; a LUT4 estimate admits only\n");
		log("strict predicted savings, counting only truly removable decoders.\n");
		log("Runs last, before technology mapping; actual mapped area must be measured.\n");
		log("The conservative admission guard skips pure-control output columns\n");
		log("and data tables directly fed by another mux, preserving those cones.\n");
		log("With -pre-factor, run one enabled common-operand/arithmetic sharing\n");
		log("round on the original decoded PMUXs before nested flattening.\n");
		log("With -flatten-decode, collapse nested decoded PMUX graphs into a\n");
		log("disjoint partition of shared symbolic vector leaves. Missing codes\n");
		log("preserve their original defaults. Raw, overlapping and X/Z selectors\n");
		log("are rejected. Bounds: 11 control bits, 32 nodes, depth 8, 32 leaves,\n");
		log("64 data bits and 65536 expanded table bits per module. Shared exact\n");
		log("Internal equality-guarded muxes with one all-X arm are traversed only\n");
		log("under PMUX roots. Any reachable X/Z leaf rejects the entire root.\n");
		log("selector truth tables use LUT4 leaves and Shannon muxes; four orders\n");
		log("and every default leaf are costed. Only disposable inner nodes count\n");
		log("towards savings; outside fanout and all arithmetic remain intact.\n");
		log("With -flatten-binary, also compare a shared Shannon DAG of the full\n");
		log("root output table, using the same bounded orders and exact defaults.\n");
		log("Three distinct same-kind arithmetic leaves preserve the PMUX form.\n");
		log("This option implies -flatten-decode. Costs remain estimates only.\n");
		log("With -factor-ops, factor a common operand of equal-width XOR/XNOR/AND/OR\n");
		log("or ADD cells across proven mutually exclusive decoded branches. Complete\n");
		log("tables factor the default as well; partial groups retain the exact A\n");
		log("and all other branches behind an outer PMUX. Requires exclusive data\n");
		log("fanout, selected non-kept operators and a positive structural estimate.\n");
		log("With -share-arith, share a single ADD or SUB across at least three\n");
		log("distinct exclusive operators by selecting both operand vectors. No\n");
		log("common operand is required; SUB order, widths and signedness remain\n");
		log("exact. Partial groups retain default and other-branch behavior. Only\n");
		log("one profitable group per PMUX per round is selected, for at most three\n");
		log("rounds with fresh graph indices. Only validated local partition proofs\n");
		log("can continue across rounds; trusted flattened selectors allow up to 32\n");
		log("choices, while ordinary decoded PMUX admission remains capped at 16.\n");
		log("\n");
		log("The pass preserves Boolean behavior and Yosys undefined-value refinement;\n");
		log("it does not promise identical four-state simulation X propagation.\n");
		log("It never treats an X/Z decoder constant as a defined bit.\n");
		log("\n");
	}

	void execute(std::vector<std::string> args, Design *design) override
	{
		log_header(design, "Executing PMUX_OPT pass (decode-conditioned data cofactoring).\n");
		bool verbose = false;
		bool rom_functional = false;
		bool share = false, priority = false, project = false, cofactor = true;
		bool prune_default = false, small_rom = false, rebuild_decode = false, factor_ops = false, flatten_decode = false;
		bool share_arith = false;
		bool pre_factor = false, flatten_binary = false;
		size_t argidx = 1;
		for (; argidx < args.size(); argidx++) {
			if (args[argidx] == "-verbose")
				verbose = true;
			else if (args[argidx] == "-share")
				share = true;
			else if (args[argidx] == "-priority")
				priority = true;
			else if (args[argidx] == "-prune-default")
				prune_default = true;
			else if (args[argidx] == "-rom-functional")
				rom_functional = true;
			else if (args[argidx] == "-small-rom")
				small_rom = true;
			else if (args[argidx] == "-project")
				project = true;
			else if (args[argidx] == "-rebuild-decode")
				rebuild_decode = true;
			else if (args[argidx] == "-flatten-decode")
				flatten_decode = true;
			else if (args[argidx] == "-flatten-binary")
				flatten_binary = flatten_decode = true;
			else if (args[argidx] == "-pre-factor")
				pre_factor = true;
			else if (args[argidx] == "-factor-ops")
				factor_ops = true;
			else if (args[argidx] == "-share-arith")
				share_arith = true;
			else if (args[argidx] == "-no-cofactor")
				cofactor = false;
			else
				break;
		}
		extra_args(args, argidx, design);
		int cells = 0, bits = 0, decodes = 0, shared_muxes = 0, shared_choices = 0;
		int priority_chains = 0, priority_branches = 0;
		int projected_muxes = 0, projected_bits = 0;
		int default_pruned_muxes = 0, default_pruned_branches = 0;
		int rebuilt_muxes = 0, rebuilt_bits = 0, rebuilt_nodes = 0, rebuilt_decoders = 0, rebuild_estimated_gain = 0;
		int rebuild_control_skipped = 0, rebuild_nonleaf_skipped = 0;
		int factored_muxes = 0, factored_operators = 0, factor_estimated_gain = 0;
		int flattened_roots = 0, flattened_muxes = 0, flattened_leaves = 0, flattened_conditions = 0, flatten_estimated_gain = 0;
		int flattened_binary_roots = 0, flattened_binary_nodes = 0;
		int shared_arith_muxes = 0, shared_arith_operators = 0, arith_estimated_gain = 0;
		int shared_arith_removed_cells = 0, shared_arith_distinct_expressions = 0;
		int noop_modules = 0;
		PmuxSmallRom::Statistics small_rom_stats;
		for (auto module : design->selected_modules()) {
			if (module->has_processes_warn())
				continue;
			// ROM rewriting adds aliases. Run it before creating any reusable
			// worker SigMap/driver index, including for ROM-only modules.
			if (small_rom) {
				auto result = PmuxSmallRom::decompose(module, verbose, rom_functional);
				small_rom_stats.roms += result.roms;
				small_rom_stats.data_bits += result.data_bits;
				small_rom_stats.luts += result.luts;
				small_rom_stats.muxes += result.muxes;
				small_rom_stats.functional_outputs += result.functional_outputs;
			}
			bool has_mux = false;
			int wide_muxes = 0, share_candidates = 0;
			bool project_candidate = false;
			bool prune_candidate = false;
			bool rebuild_candidate = false;
			bool factor_candidate = false;
			for (auto cell : module->selected_cells()) {
				if (!cell->type.in(ID($mux), ID($pmux)) || cell->get_bool_attribute(ID::keep))
					continue;
				has_mux = true;
				if (priority && cell->type == ID($mux) && GetSize(cell->getPort(ID::Y)) >= 4)
					wide_muxes++;
				if (wide_muxes >= 3)
					break;
				if (cell->type != ID($pmux))
					continue;
				int width = GetSize(cell->getPort(ID::A));
				int choices = GetSize(cell->getPort(ID::S));
				if (share && width >= 2 && choices >= 4)
					share_candidates++;
				// Data need not yet be constant: cofactoring runs before
				// pruning/projection and can expose eligible constant slices.
				if (prune_default && width >= 1 && width <= 64 && choices >= 2 && choices <= 8 &&
						GetSize(cell->getPort(ID::B)) == width * choices && GetSize(cell->getPort(ID::Y)) == width)
					prune_candidate = true;
				if (project && width >= 2 && width <= 64 && choices >= 2 && choices <= 8 &&
						GetSize(cell->getPort(ID::B)) == width * choices && GetSize(cell->getPort(ID::Y)) == width)
					project_candidate = true;
				if ((rebuild_decode || flatten_decode) && width >= 1 && width <= 64 && choices >= 2 && choices <= 16 &&
						GetSize(cell->getPort(ID::B)) == width * choices && GetSize(cell->getPort(ID::Y)) == width)
					rebuild_candidate = true;
				if ((factor_ops || share_arith) && width >= 1 && width <= 64 && choices >= 1 && choices <= 16 &&
						GetSize(cell->getPort(ID::B)) == width * choices && GetSize(cell->getPort(ID::Y)) == width)
					factor_candidate = true;
				if (share_candidates >= 2 || project_candidate || prune_candidate || rebuild_candidate || factor_candidate)
					break;
			}
			// Technology libraries contain many empty/blackbox modules. Do
			// not build analysis tables for modules with no eligible mux.
			if (!has_mux || (!cofactor && !share && !priority && !project && !prune_default && !rebuild_decode && !factor_ops && !flatten_decode && !share_arith))
				continue;
			PmuxOptWorker worker(module, verbose);
			// Priority packing can create PMUXs for later phases, so any
			// possible structural rewrite keeps the complete analysis path.
			bool structural_candidate = wide_muxes >= 3 || share_candidates >= 2 || project_candidate || prune_candidate || rebuild_candidate || factor_candidate;
			if (!structural_candidate && (!cofactor || !worker.may_cofactor())) {
				noop_modules++;
				if (verbose)
					log("  Skipped %s: no structural candidate or cofactor dependency.\n", log_id(module));
				continue;
			}
			worker.index_drivers();
			if (priority)
				worker.pack_priority();
			if (cofactor)
				worker.run();
			if (prune_default)
				worker.prune_default_equivalent_branches();
			if (project)
				worker.project_constants();
			if (share)
				worker.share_decoded();
			cells += worker.changed_cells;
			bits += worker.replaced_bits;
			decodes += worker.decoded_conditions;
			shared_muxes += worker.shared_muxes;
			shared_choices += worker.shared_data_choices;
			priority_chains += worker.priority_chains;
			priority_branches += worker.priority_branches;
			projected_muxes += worker.projected_muxes;
			projected_bits += worker.projected_data_bits;
			default_pruned_muxes += worker.default_pruned_muxes;
			default_pruned_branches += worker.default_pruned_branches;
			int module_prefactored = 0;
			std::vector<PmuxOptWorker::PartitionProof> prefactor_proofs;
			if (pre_factor && (factor_ops || share_arith)) {
				PmuxOptWorker pre(module, verbose);
				pre.index_drivers();
				pre.factor_decoded_ops({}, factor_ops, share_arith);
				module_prefactored = pre.factored_muxes + pre.shared_arith_muxes;
				factored_muxes += pre.factored_muxes;
				factored_operators += pre.factored_operators;
				factor_estimated_gain += pre.factor_estimated_gain;
				shared_arith_muxes += pre.shared_arith_muxes;
				shared_arith_operators += pre.shared_arith_operators;
				shared_arith_removed_cells += pre.shared_arith_removed_cells;
				shared_arith_distinct_expressions += pre.shared_arith_distinct_expressions;
				arith_estimated_gain += pre.arith_estimated_gain;
				prefactor_proofs = std::move(pre.continuation_proofs);
			}
			std::vector<PmuxOptWorker::PartitionProof> partition_proofs;
			if (flatten_decode) {
				PmuxOptWorker fresh(module, verbose);
				fresh.index_drivers();
				fresh.flatten_decoded(partition_proofs, flatten_binary);
				flattened_binary_roots += fresh.flattened_binary_roots;
				flattened_binary_nodes += fresh.flattened_binary_nodes;
				flattened_roots += fresh.flattened_roots;
				flattened_muxes += fresh.flattened_muxes;
				flattened_leaves += fresh.flattened_leaves;
				flattened_conditions += fresh.flattened_conditions;
				flatten_estimated_gain += fresh.flatten_estimated_gain;
			}
			// Flatten may delete a previously factored root. Verify name and
			// identity before handing its snapshot to the next fresh worker.
			// read_partition_proof checks its ports and parameters again.
			for (auto &proof : prefactor_proofs)
				if (module->cell(proof.name) == proof.cell)
					partition_proofs.push_back(std::move(proof));
			int module_factored = module_prefactored;
			if (factor_ops || share_arith) {
				auto collect_factor = [&](const PmuxOptWorker &result) {
					module_factored += result.factored_muxes + result.shared_arith_muxes;
					factored_muxes += result.factored_muxes;
					factored_operators += result.factored_operators;
					factor_estimated_gain += result.factor_estimated_gain;
					shared_arith_muxes += result.shared_arith_muxes;
					shared_arith_operators += result.shared_arith_operators;
					shared_arith_removed_cells += result.shared_arith_removed_cells;
					shared_arith_distinct_expressions += result.shared_arith_distinct_expressions;
					arith_estimated_gain += result.arith_estimated_gain;
				};
				std::vector<PmuxOptWorker::PartitionProof> proofs = std::move(partition_proofs);
				for (int round = 0; round < 3; ++round) {
					PmuxOptWorker fresh(module, verbose);
					fresh.index_drivers();
					fresh.factor_decoded_ops(proofs, factor_ops, share_arith);
					collect_factor(fresh);
					if (!fresh.factored_muxes && !fresh.shared_arith_muxes)
						break;
					proofs = std::move(fresh.continuation_proofs);
				}
			}
			if (rebuild_decode) {
				auto collect_rebuild = [&](const PmuxOptWorker &result) {
					rebuilt_muxes += result.rebuilt_muxes;
					rebuilt_bits += result.rebuilt_data_bits;
					rebuilt_nodes += result.rebuilt_nodes;
					rebuilt_decoders += result.rebuilt_decoders;
					rebuild_estimated_gain += result.rebuild_estimated_gain;
					rebuild_control_skipped += result.rebuild_control_skipped;
					rebuild_nonleaf_skipped += result.rebuild_nonleaf_skipped;
				};
				if (worker.shared_muxes || module_factored || flatten_decode) {
					// Sharing creates aliases and deletes drivers. Only actual
					// sharing requires a fresh index; other phases maintain it.
					PmuxOptWorker fresh(module, verbose);
					fresh.index_drivers();
					fresh.rebuild_decoded();
					collect_rebuild(fresh);
				} else {
					worker.rebuild_decoded();
					collect_rebuild(worker);
				}
			}
		}
		design->scratchpad_set_int("pmux_opt.changed_cells", cells);
		design->scratchpad_set_int("pmux_opt.replaced_bits", bits);
		design->scratchpad_set_int("pmux_opt.decoded_conditions", decodes);
		design->scratchpad_set_int("pmux_opt.rom_functional_outputs", small_rom_stats.functional_outputs);
		design->scratchpad_set_int("pmux_opt.shared_muxes", shared_muxes);
		design->scratchpad_set_int("pmux_opt.shared_data_choices", shared_choices);
		design->scratchpad_set_int("pmux_opt.priority_chains", priority_chains);
		design->scratchpad_set_int("pmux_opt.priority_branches", priority_branches);
		design->scratchpad_set_int("pmux_opt.projected_muxes", projected_muxes);
		design->scratchpad_set_int("pmux_opt.projected_data_bits", projected_bits);
		design->scratchpad_set_int("pmux_opt.default_pruned_muxes", default_pruned_muxes);
		design->scratchpad_set_int("pmux_opt.default_pruned_branches", default_pruned_branches);
		design->scratchpad_set_int("pmux_opt.rebuilt_muxes", rebuilt_muxes);
		design->scratchpad_set_int("pmux_opt.rebuilt_data_bits", rebuilt_bits);
		design->scratchpad_set_int("pmux_opt.rebuilt_nodes", rebuilt_nodes);
		design->scratchpad_set_int("pmux_opt.rebuilt_decoders", rebuilt_decoders);
		design->scratchpad_set_int("pmux_opt.rebuild_estimated_gain", rebuild_estimated_gain);
		design->scratchpad_set_int("pmux_opt.rebuild_control_skipped", rebuild_control_skipped);
		design->scratchpad_set_int("pmux_opt.rebuild_nonleaf_skipped", rebuild_nonleaf_skipped);
		design->scratchpad_set_int("pmux_opt.factored_muxes", factored_muxes);
		design->scratchpad_set_int("pmux_opt.factored_operators", factored_operators);
		design->scratchpad_set_int("pmux_opt.factor_estimated_gain", factor_estimated_gain);
		design->scratchpad_set_int("pmux_opt.shared_arith_muxes", shared_arith_muxes);
		design->scratchpad_set_int("pmux_opt.shared_arith_operators", shared_arith_operators);
		design->scratchpad_set_int("pmux_opt.arith_estimated_gain", arith_estimated_gain);
		design->scratchpad_set_int("pmux_opt.shared_arith_removed_cells", shared_arith_removed_cells);
		design->scratchpad_set_int("pmux_opt.shared_arith_distinct_expressions", shared_arith_distinct_expressions);
		design->scratchpad_set_int("pmux_opt.flattened_binary_roots", flattened_binary_roots);
		design->scratchpad_set_int("pmux_opt.flattened_binary_nodes", flattened_binary_nodes);
		design->scratchpad_set_int("pmux_opt.flattened_roots", flattened_roots);
		design->scratchpad_set_int("pmux_opt.flattened_muxes", flattened_muxes);
		design->scratchpad_set_int("pmux_opt.flattened_leaves", flattened_leaves);
		design->scratchpad_set_int("pmux_opt.flattened_conditions", flattened_conditions);
		design->scratchpad_set_int("pmux_opt.flatten_estimated_gain", flatten_estimated_gain);
		design->scratchpad_set_int("pmux_opt.noop_modules", noop_modules);
		design->scratchpad_set_int("pmux_opt.small_roms", small_rom_stats.roms);
		design->scratchpad_set_int("pmux_opt.small_rom_data_bits", small_rom_stats.data_bits);
		design->scratchpad_set_int("pmux_opt.small_rom_luts", small_rom_stats.luts);
		design->scratchpad_set_int("pmux_opt.small_rom_muxes", small_rom_stats.muxes);
		if (bits || shared_muxes || priority_chains || projected_muxes || default_pruned_branches || small_rom_stats.roms || rebuilt_muxes || factored_muxes || flattened_roots || shared_arith_muxes)
			design->scratchpad_set_bool("opt.did_something", true);
		log("PMUX_OPT: replaced %d data bits in %d cells; inspected %d decoded conditions.\n", bits, cells, decodes);
		if (rom_functional) log("PMUX_OPT: exact two-LUT decomposition for %d ROM outputs.\n", small_rom_stats.functional_outputs);
		if (share)
			log("PMUX_OPT: shared %d decoded muxes and %d residual data choices.\n", shared_muxes, shared_choices);
		if (priority)
			log("PMUX_OPT: packed %d priority chains with %d branches.\n", priority_chains, priority_branches);
		if (project)
			log("PMUX_OPT: projected %d constant-data muxes and %d data bits.\n", projected_muxes, projected_bits);
		if (small_rom)
			log("PMUX_OPT: decomposed %d small ROMs (%d data bits), emitted %d LUTs and %d muxes.\n",
				small_rom_stats.roms, small_rom_stats.data_bits, small_rom_stats.luts, small_rom_stats.muxes);
		if (prune_default)
			log("PMUX_OPT: pruned %d default-equivalent branches in %d muxes.\n", default_pruned_branches, default_pruned_muxes);
		if (rebuild_decode)
			log("PMUX_OPT: rebuilt %d decoded muxes (%d data bits), emitted %d mux nodes, removed %d decoders; estimated gain %d.\n",
				rebuilt_muxes, rebuilt_bits, rebuilt_nodes, rebuilt_decoders, rebuild_estimated_gain);
		if (rebuild_decode)
			log("PMUX_OPT: conservative rebuild guard skipped %d control-only and %d non-leaf muxes.\n",
				rebuild_control_skipped, rebuild_nonleaf_skipped);
		if (factor_ops)
			log("PMUX_OPT: factored %d muxes, removed %d duplicate operators; estimated gain %d.\n",
				factored_muxes, factored_operators, factor_estimated_gain);
		if (share_arith)
			log("PMUX_OPT: shared arithmetic in %d muxes; removed %d old cells, net %d fewer operators; %d distinct expressions, estimated gain %d.\n",
				shared_arith_muxes, shared_arith_removed_cells, shared_arith_operators, shared_arith_distinct_expressions, arith_estimated_gain);
		if (flatten_binary)
			log("PMUX_OPT: binary alternative selected %d roots and emitted %d mux nodes.\n",
				flattened_binary_roots, flattened_binary_nodes);
		if (flatten_decode)
			log("PMUX_OPT: flattened %d roots across %d disposable decoded muxes into %d leaves and %d selector nodes; estimated gain %d.\n",
				flattened_roots, flattened_muxes, flattened_leaves, flattened_conditions, flatten_estimated_gain);
	}
} PmuxOptPass;

PRIVATE_NAMESPACE_END
