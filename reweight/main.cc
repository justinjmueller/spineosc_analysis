/*
 * @file main.cc
 * @brief Main entry point for the event reweighting code.
 * @details This code is intended to produce reweighted copies of the SPINE
 * ntuples used by the PROfit fits. At this stage it does the copying half of
 * that job only: it clones a sample's `selected` tree together with the
 * multisim and multisigma weight trees into a new file, preserving the
 * events/{key}/ directory layout the XML configs reference, and adds one
 * per-event column to the event tree: `selected`, the post-selection cuts
 * evaluated as a weight. The weights are carried through untouched, so the
 * output is a drop-in replacement for the input as far as PROfit is
 * concerned, and the reweighting itself can be layered on top of this without
 * disturbing the file plumbing.
 * @author Justin Mueller (FNAL)
 */
#include <TBranch.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TTree.h>
#include <TTreeFormula.h>

#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "configuration.h"

namespace reweight
{
    /**
     * @brief The trees copied out of a sample, relative to events/{key}/.
     * @details `selected` is the event tree the fits read; the other two are
     * the per-event weight trees the XML configs attach to it as friends.
     * All three are row-parallel -- one entry per selected event, in the same
     * order -- which is what makes the friend relationship valid, and what
     * lets a truncated copy stay self-consistent as long as every tree is cut
     * at the same entry count.
     */
    constexpr const char *k_selected_tree = "selected";
    constexpr const char *k_multisigma_tree = "selected_multisigmaTree";

    const std::vector<std::string> tree_names = {
        k_selected_tree,
        k_multisigma_tree,
        "selected_multisimTree",
    };

    /**
     * @brief Opens a ROOT file, exiting the process if it cannot be read.
     * @param path The path of the ROOT file to open.
     * @param mode The mode to open the file in, e.g. "READ" or "RECREATE".
     * @return A pointer to the newly-opened TFile (caller takes ownership).
     */
    TFile *open_file(const std::string &path, const std::string &mode)
    {
        TFile *file = TFile::Open(path.c_str(), mode.c_str());
        if (!file || file->IsZombie())
        {
            std::cerr << "ERROR: could not open " << path << std::endl;
            std::exit(1);
        }
        return file;
    }

    /**
     * @brief Fetches a tree from an open file, or nullptr if it is absent.
     * @param file The ROOT file to read from.
     * @param tree_path The in-file path of the tree, e.g. "events/sbnd/selected".
     * @return A pointer to the tree (owned by the file), or nullptr.
     */
    TTree *find_tree(TFile *file, const std::string &tree_path)
    {
        return static_cast<TTree *>(file->Get(tree_path.c_str()));
    }

    /**
     * @brief Fetches a tree from an open file, exiting if it is missing.
     * @param file The ROOT file to read from.
     * @param tree_path The in-file path of the tree, e.g. "events/sbnd/selected".
     * @return A pointer to the tree (owned by the file).
     */
    TTree *get_tree(TFile *file, const std::string &tree_path)
    {
        TTree *tree = find_tree(file, tree_path);
        if (!tree)
        {
            std::cerr << "ERROR: missing tree " << tree_path << std::endl;
            std::exit(1);
        }
        return tree;
    }

    /**
     * @brief Looks up a sample's configuration by its tree key.
     * @details An unrecognised key is fatal rather than defaulted: the
     * calibration and thresholds differ per sample, so guessing one would
     * quietly produce a selection that is wrong for the data in hand.
     * @param key The sample's tree key, e.g. "icarus_offbeam_run2".
     * @return The matching sample configuration.
     */
    SampleConfig find_sample(const std::string &key)
    {
        const std::vector<SampleConfig> samples = build_samples();
        for (const auto &sample : samples)
        {
            if (sample.key == key) return sample;
        }

        std::cerr << "ERROR: no configuration for sample \"" << key << "\". Known samples:" << std::endl;
        for (const auto &sample : samples) std::cerr << "  " << sample.key << std::endl;
        std::exit(1);
    }

    /**
     * @brief Clones one tree into a directory of the output file.
     * @details A full copy is fast-cloned, which moves compressed baskets
     * across without unzipping and rezipping them. A truncated copy cannot
     * take that path -- fast cloning works in whole baskets and would
     * overshoot the requested entry count -- so it falls back to an
     * entry-by-entry copy, which is slower but gives exactly max_entries.
     * @param tree The input tree to clone.
     * @param directory The output directory to clone into.
     * @param max_entries The number of entries to copy, or -1 for all of them.
     * @return The cloned tree, not yet written (owned by the directory).
     */
    TTree *clone_tree(TTree *tree, TDirectory *directory, Long64_t max_entries)
    {
        directory->cd();
        const bool full_copy = max_entries < 0;
        TTree *clone = tree->CloneTree(max_entries, full_copy ? "fast" : "");
        if (!clone)
        {
            std::cerr << "ERROR: could not clone " << tree->GetName() << std::endl;
            std::exit(1);
        }
        return clone;
    }

    /**
     * @brief A configured cut paired with the formula compiled from it.
     */
    struct CompiledCut
    {
        std::string name;       ///< The cut's name, for the summary line.
        TTreeFormula *formula;  ///< The compiled expression (owned by this struct).
        Long64_t n_failing;     ///< Events this cut rejects.
    };

    /**
     * @brief Compiles every enabled cut against the input tree.
     * @details An expression naming a branch that does not exist leaves
     * TTreeFormula with zero dimensions, which is checked here rather than
     * left to fill a column of silent zeros -- indistinguishable, downstream,
     * from a selection nothing passes. Compiling before the output file is
     * created means a bad expression fails without leaving a copy behind.
     *
     * Disabled cuts are skipped rather than compiled, so a cut can be turned
     * off in the configuration even when the sample in hand has no branch to
     * evaluate it against.
     * @param source The input tree the formulas are evaluated against.
     * @param selection The cuts to compile.
     * @return The compiled cuts (caller owns the formulas; see free_cuts()).
     */
    std::vector<CompiledCut> compile_cuts(TTree *source, const SelectionConfig &selection)
    {
        std::vector<CompiledCut> compiled;
        for (const auto &cut : selection.cuts)
        {
            if (!cut.enabled)
            {
                std::cout << "  " << cut.name << " DISABLED" << std::endl;
                continue;
            }

            TTreeFormula *formula = new TTreeFormula(cut.name.c_str(), cut.expression.c_str(), source);
            if (formula->GetNdim() == 0)
            {
                std::cerr << "ERROR: could not compile cut " << cut.name << " \"" << cut.expression
                          << "\" against " << source->GetName() << std::endl;
                std::exit(1);
            }
            std::cout << "  " << cut.name << " " << cut.expression << std::endl;
            compiled.push_back({cut.name, formula, 0});
        }

        if (compiled.empty())
        {
            std::cerr << "ERROR: every cut is disabled, so the selection would be vacuous" << std::endl;
            std::exit(1);
        }
        return compiled;
    }

    /**
     * @brief Frees the formulas held by a compiled cut list.
     * @param cuts The compiled cuts to free.
     */
    void free_cuts(std::vector<CompiledCut> &cuts)
    {
        for (auto &cut : cuts) delete cut.formula;
        cuts.clear();
    }

    /**
     * @brief Evaluates the configured selection into a new branch of a clone.
     * @details The branch is added after the tree has been cloned, so the
     * copy itself still gets the cheap basket-level path and only this one
     * branch is built by looping. The formula is evaluated against the INPUT
     * tree -- the clone's branches hold no addresses to read through -- and
     * the results are filled in entry order, which keeps the new branch
     * aligned with the rows already copied.
     *
     * Every cut is evaluated for every event rather than short-circuiting on
     * the first failure, so each cut's rejection count is its own marginal
     * effect. Those counts overlap and will not sum to the number of events
     * rejected overall, which is the point: it shows what each cut costs
     * independently of the order they happen to be listed in.
     * @param source The input tree the formulas are evaluated against.
     * @param clone The cloned tree to add the branch to.
     * @param selection The branch name to write.
     * @param cuts The compiled cuts, whose rejection counts this updates.
     * @param n_entries The number of entries to evaluate, matching the clone.
     * @return The number of entries passing every cut.
     */
    Long64_t add_selection_branch(TTree *source, TTree *clone, const SelectionConfig &selection,
                                  std::vector<CompiledCut> &cuts, Long64_t n_entries)
    {
        double value = 0.0;
        TBranch *branch = clone->Branch(selection.name.c_str(), &value, (selection.name + "/D").c_str());
        Long64_t n_passing = 0;
        for (Long64_t entry = 0; entry < n_entries; ++entry)
        {
            source->GetEntry(entry);

            bool passes = true;
            for (auto &cut : cuts)
            {
                // Required before EvalInstance: it loads the formula's leaves
                // for this entry and reports how many values it has.
                cut.formula->GetNdata();
                if (cut.formula->EvalInstance() == 0.0)
                {
                    ++cut.n_failing;
                    passes = false;
                }
            }

            value = passes ? 1.0 : 0.0;
            if (passes) ++n_passing;
            branch->Fill();
        }
        return n_passing;
    }

    /**
     * @brief The index of a knob within a systematic's per-event weight vector.
     * @details The knob values live in a parallel `<branch>_sigma` vector that
     * is the same for every event, so it is read once from entry 0. Returns -1
     * if the branch has no knob at the requested sigma, which is how a
     * systematic that carries no central value is distinguished from one that
     * does.
     * @param tree The multisigma tree.
     * @param branch The weight branch's name.
     * @param knob The sigma value to locate.
     * @return The index into the weight vector, or -1 if absent.
     */
    int find_knob_index(TTree *tree, const std::string &branch, double knob)
    {
        const std::string sigma_branch = branch + "_sigma";
        if (!tree->GetBranch(sigma_branch.c_str())) return -1;

        std::vector<double> *sigmas = nullptr;
        tree->SetBranchAddress(sigma_branch.c_str(), &sigmas);
        tree->GetEntry(0);

        int index = -1;
        if (sigmas)
        {
            for (std::size_t i = 0; i < sigmas->size(); ++i)
            {
                if ((*sigmas)[i] == knob)
                {
                    index = static_cast<int>(i);
                    break;
                }
            }
        }
        tree->ResetBranchAddresses();
        return index;
    }

    /**
     * @brief A CV reweight resolved against the multisigma tree in hand.
     */
    struct ResolvedReweight
    {
        std::string branch;             ///< The weight branch's name.
        int index;                      ///< Knob index, or -1 if the knob is absent.
        std::vector<double> *weights;   ///< Read buffer (owned by the tree).
        Long64_t n_guarded;             ///< Events whose weight was unusable and forced to 1.
    };

    /**
     * @brief Resolves the configured CV reweights against a multisigma tree.
     * @details A required systematic whose branch or knob is missing is fatal:
     * silently dropping a central-value correction would leave a prediction
     * that looks fine and is wrong by however much that correction was worth.
     * An optional one is reported and contributes a factor of 1.
     * @param tree The multisigma tree to read from.
     * @param config The reweighting configuration.
     * @return The resolved reweights, with read buffers attached to the tree.
     */
    std::vector<ResolvedReweight> resolve_reweights(TTree *tree, const ReweightConfig &config)
    {
        std::vector<ResolvedReweight> resolved;
        for (const auto &reweight : config.reweights)
        {
            const bool has_branch = tree->GetBranch(reweight.branch.c_str()) != nullptr;
            const int index = has_branch ? find_knob_index(tree, reweight.branch, reweight.knob) : -1;

            if (index < 0)
            {
                const std::string why = has_branch ? "no knob at sigma " + std::to_string(reweight.knob)
                                                   : "branch absent";
                if (reweight.required)
                {
                    std::cerr << "ERROR: required CV reweight " << reweight.branch << ": " << why << std::endl;
                    std::exit(1);
                }
                std::cout << "  " << reweight.branch << " SKIPPED (" << why << "), contributes 1" << std::endl;
                continue;
            }

            // Addresses are attached later, by attach_reweights(): resolving
            // the NEXT systematic calls find_knob_index(), which resets every
            // address on the tree, and the read buffer must outlive this loop
            // in storage that will not move again.
            std::cout << "  " << reweight.branch << " CV knob at index " << index << std::endl;
            resolved.push_back({reweight.branch, index, nullptr, 0});
        }
        return resolved;
    }

    /**
     * @brief Points each resolved reweight at its branch in the tree.
     * @details Split from resolve_reweights() because SetBranchAddress stores
     * the address of the member it is handed. Taking that address while the
     * vector is still being appended to would leave it dangling the moment a
     * push_back reallocated, and resolving a later systematic would reset it
     * regardless. Call this only once the list is final.
     * @param tree The multisigma tree to read through.
     * @param resolved The resolved reweights to attach, modified in place.
     */
    void attach_reweights(TTree *tree, std::vector<ResolvedReweight> &resolved)
    {
        for (auto &reweight : resolved)
            tree->SetBranchAddress(reweight.branch.c_str(), &reweight.weights);
    }

    /**
     * @brief Computes the per-event CV weight and adds it to the event tree.
     * @details The weight is the product of one knob per configured
     * systematic, read from the multisigma tree, which is row-parallel with
     * `selected` -- the same friend relationship the XML configs rely on, used
     * here directly.
     *
     * A weight that is not finite or not positive is forced to 1 and counted
     * rather than propagated: a zero would delete the event from the
     * prediction and a NaN would poison every bin it touched, and both are
     * input defects worth seeing a number for.
     * @param weights_tree The multisigma tree holding the weight vectors.
     * @param clone The cloned `selected` tree to add the column to.
     * @param config The reweighting configuration.
     * @param resolved The resolved reweights, whose guard counts this updates.
     * @param n_entries The number of entries to evaluate, matching the clone.
     * @return The mean weight over all entries, for the summary line.
     */
    double add_cv_weight_branch(TTree *weights_tree, TTree *clone, const ReweightConfig &config,
                                std::vector<ResolvedReweight> &resolved, Long64_t n_entries)
    {
        // Attached here, where the list is final and will not be reallocated.
        attach_reweights(weights_tree, resolved);

        double value = 1.0;
        TBranch *branch = clone->Branch(config.weight_branch.c_str(), &value,
                                        (config.weight_branch + "/D").c_str());
        double total = 0.0;
        for (Long64_t entry = 0; entry < n_entries; ++entry)
        {
            weights_tree->GetEntry(entry);

            value = 1.0;
            for (auto &reweight : resolved)
            {
                double w = 1.0;
                if (reweight.weights && reweight.index < static_cast<int>(reweight.weights->size()))
                    w = (*reweight.weights)[reweight.index];

                if (!std::isfinite(w) || w <= 0.0)
                {
                    ++reweight.n_guarded;
                    w = 1.0;
                }
                value *= w;
            }

            total += value;
            branch->Fill();
        }
        return n_entries > 0 ? total / static_cast<double>(n_entries) : 1.0;
    }

    /**
     * @brief Clones the multisigma tree, renormalising the configured branches.
     * @details Every branch named in `normalise` is divided, knob by knob, by
     * its own CV entry, so the resulting spline passes through exactly 1 at
     * zero sigma. The +-sigma variations survive as relative variations, which
     * is what keeps the component fittable after its central value has been
     * moved into the event weight.
     *
     * The rewritten branches are excluded from the clone and rebuilt
     * afterwards under the same names, because a cloned branch holds copied
     * baskets that cannot be edited in place. Everything else still travels as
     * a straight copy.
     * @param tree The input multisigma tree.
     * @param directory The output directory to clone into.
     * @param max_entries Entries to copy, or -1 for all of them.
     * @param config The reweighting configuration.
     * @return The cloned tree, not yet written (owned by the directory).
     */
    TTree *clone_weights_tree(TTree *tree, TDirectory *directory, Long64_t max_entries,
                              const ReweightConfig &config)
    {
        // Resolve first: a branch configured for renormalisation that this
        // sample does not have is skipped rather than fatal, so one
        // configuration covers samples with differing weight sets.
        std::vector<std::pair<std::string, int>> targets;
        for (const auto &entry : config.normalise)
        {
            if (!tree->GetBranch(entry.branch.c_str())) continue;
            const int index = find_knob_index(tree, entry.branch, entry.knob);
            if (index < 0)
            {
                std::cout << "  " << entry.branch << " has no knob at sigma " << entry.knob
                          << ", left as-is" << std::endl;
                continue;
            }
            targets.emplace_back(entry.branch, index);
        }

        for (const auto &[branch, index] : targets)
        {
            (void)index;
            tree->SetBranchStatus(branch.c_str(), 0);
        }

        // Fast cloning moves whole baskets and cannot honour a branch
        // selection, so a renormalising copy takes the slower entry-by-entry
        // path. A copy with nothing to rewrite keeps the fast path.
        directory->cd();
        const bool rewriting = !targets.empty();
        const bool full_copy = max_entries < 0 && !rewriting;
        TTree *clone = tree->CloneTree(max_entries, full_copy ? "fast" : "");
        if (!clone)
        {
            std::cerr << "ERROR: could not clone " << tree->GetName() << std::endl;
            std::exit(1);
        }
        if (!rewriting) return clone;

        // Re-enable for reading: branch status gates GetEntry as well as
        // CloneTree, so leaving these off would read nothing back.
        for (const auto &[branch, index] : targets)
        {
            (void)index;
            tree->SetBranchStatus(branch.c_str(), 1);
        }

        const Long64_t n_entries = clone->GetEntries();
        std::vector<std::vector<double> *> sources(targets.size(), nullptr);
        std::vector<std::vector<double>> buffers(targets.size());
        std::vector<Long64_t> guarded(targets.size(), 0);
        for (std::size_t i = 0; i < targets.size(); ++i)
        {
            tree->SetBranchAddress(targets[i].first.c_str(), &sources[i]);
            clone->Branch(targets[i].first.c_str(), &buffers[i]);
        }

        for (Long64_t entry = 0; entry < n_entries; ++entry)
        {
            tree->GetEntry(entry);
            for (std::size_t i = 0; i < targets.size(); ++i)
            {
                buffers[i].clear();
                if (!sources[i]) continue;

                const int index = targets[i].second;
                const bool usable = index < static_cast<int>(sources[i]->size()) &&
                                    std::isfinite((*sources[i])[index]) && (*sources[i])[index] > 0.0;
                // An unusable divisor leaves the universes untouched rather
                // than producing infinities; it is counted so a systematically
                // bad input cannot pass for a clean one.
                const double divisor = usable ? (*sources[i])[index] : 1.0;
                if (!usable) ++guarded[i];

                buffers[i].reserve(sources[i]->size());
                for (double w : *sources[i]) buffers[i].push_back(w / divisor);
            }
            for (std::size_t i = 0; i < targets.size(); ++i)
                clone->GetBranch(targets[i].first.c_str())->Fill();
        }

        for (std::size_t i = 0; i < targets.size(); ++i)
        {
            std::cout << "  " << targets[i].first << " renormalised at knob index " << targets[i].second;
            if (guarded[i] > 0) std::cout << " (" << guarded[i] << " entries had an unusable divisor)";
            std::cout << std::endl;
        }
        tree->ResetBranchAddresses();
        return clone;
    }

    /**
     * @brief Copies a sample's selected and weight trees into a new file.
     * @details The output keeps the input's events/{key}/ layout, so the
     * `treename` and `friend` paths in the XML configs resolve against the
     * copy unchanged. The output file inherits the input's compression
     * settings, both to keep the copy the same size as the original and
     * because fast cloning is only possible when the two agree.
     * @param input_path Path to the ROOT file to copy from.
     * @param output_path Path of the ROOT file to create.
     * @param key The sample's tree key, e.g. "sbnd" or "icarus_run2".
     * @param max_entries The number of entries to copy per tree, or -1 for all.
     */
    void copy_sample(const std::string &input_path, const std::string &output_path, const std::string &key,
                     Long64_t max_entries)
    {
        const SampleConfig sample = find_sample(key);
        TFile *input = open_file(input_path, "READ");

        // Every tree is resolved before the output file is created, so a
        // sample missing one of the weight trees fails without leaving a
        // half-written copy holding `selected` alone behind.
        //
        // The weight trees are required only of the samples configured to
        // have them. The off-beam samples are data-driven cosmics with no
        // generator weights to vary, so their absence there is the expected
        // state and not a failure -- but it is still reported, and a weight
        // tree that does turn up in such a sample is copied rather than
        // quietly dropped.
        std::vector<std::pair<std::string, TTree *>> trees;
        for (const auto &tree_name : tree_names)
        {
            const std::string tree_path = "events/" + key + "/" + tree_name;
            const bool required = tree_name == k_selected_tree || sample.expects_weight_trees;

            TTree *tree = required ? get_tree(input, tree_path) : find_tree(input, tree_path);
            if (tree)
                trees.emplace_back(tree_path, tree);
            else
                std::cout << tree_path << " absent (expected for this sample)" << std::endl;
        }

        // Compiled here, against the input tree and before the output file
        // exists, so a broken expression fails on its own rather than after a
        // copy has been started.
        const SelectionConfig selection = build_selection(sample);
        std::cout << "selection \"" << selection.name << "\":" << std::endl;
        std::vector<CompiledCut> cuts = compile_cuts(trees.front().second, selection);

        TFile *output = open_file(output_path, "RECREATE");
        output->SetCompressionSettings(input->GetCompressionSettings());

        // Built one level at a time rather than as "events/<key>": a slashed
        // path does create the nesting, but hands back the TOP directory, so
        // the trees would be written to events/ and the XML configs' friend
        // paths would no longer resolve.
        TDirectory *events_directory = output->mkdir("events");
        TDirectory *directory = events_directory ? events_directory->mkdir(key.c_str()) : nullptr;
        if (!directory)
        {
            std::cerr << "ERROR: could not create events/" << key << " in " << output_path << std::endl;
            std::exit(1);
        }

        // The multisigma tree supplies the CV weights that the event tree's
        // new column is built from, so it is located before the copy loop
        // rather than waiting for its turn in it.
        TTree *weights_tree = nullptr;
        for (const auto &[tree_path, tree] : trees)
        {
            (void)tree_path;
            if (std::string(tree->GetName()) == k_multisigma_tree) weights_tree = tree;
        }

        const ReweightConfig reweight = build_reweight();
        std::vector<ResolvedReweight> resolved;
        if (weights_tree)
        {
            std::cout << "central-value reweights folded into \"" << reweight.weight_branch << "\":" << std::endl;
            resolved = resolve_reweights(weights_tree, reweight);
        }
        else
        {
            std::cout << "no multisigma tree, so no central-value reweighting" << std::endl;
        }

        for (const auto &[tree_path, tree] : trees)
        {
            const bool is_selected_tree = std::string(tree->GetName()) == k_selected_tree;
            const bool is_weights_tree = std::string(tree->GetName()) == k_multisigma_tree;

            TTree *clone = is_weights_tree ? clone_weights_tree(tree, directory, max_entries, reweight)
                                           : clone_tree(tree, directory, max_entries);
            const Long64_t copied = clone->GetEntries();

            // Only the event tree carries the selection. The weight trees are
            // per-event universes of the same rows, so a flag on them would be
            // the same column stored three times.
            std::string selection_summary;
            if (is_selected_tree)
            {
                const Long64_t n_passing = add_selection_branch(tree, clone, selection, cuts, copied);
                selection_summary = ", " + std::to_string(n_passing) + " passing " + selection.name;

                // Built from the multisigma tree, which is row-parallel with
                // this one, so entry N of each describes the same event.
                if (weights_tree && !resolved.empty())
                {
                    const double mean = add_cv_weight_branch(weights_tree, clone, reweight, resolved, copied);
                    selection_summary += ", mean " + reweight.weight_branch + " " + std::to_string(mean);
                    weights_tree->ResetBranchAddresses();
                }
            }

            clone->Write("", TObject::kOverwrite);
            std::cout << tree_path << " " << copied << " of " << tree->GetEntries() << " entries"
                      << selection_summary << std::endl;

            // Per-cut rejection counts, printed once the selection has run.
            if (is_selected_tree)
            {
                for (const auto &cut : cuts)
                    std::cout << "  " << cut.name << " rejects " << cut.n_failing << std::endl;
            }
        }

        // Freed before the file they read through is closed.
        free_cuts(cuts);

        output->Close();
        input->Close();
        delete output;
        delete input;
    }
}  // namespace reweight

/**
 * @brief Program entry point.
 * @details Copies one sample's selected tree and its multisim/multisigma
 * weight trees from the input file into a newly-created output file.
 * @param argc Argument count.
 * @param argv Argument values: <input file> <output file> [tree key] [max entries].
 * Tree key defaults to "sbnd"; max entries defaults to -1, meaning every entry.
 * @return 0 on success, 1 if the arguments are not usable.
 */
int main(int argc, char **argv)
{
    if (argc < 3)
    {
        std::cerr << "Usage: " << argv[0] << " <input file> <output file> [tree key] [max entries]" << std::endl;
        std::cerr << "  tree key     sample key under events/, default \"sbnd\"" << std::endl;
        std::cerr << "  max entries  entries to copy per tree, default -1 (all)" << std::endl;
        return 1;
    }

    const std::string input_path = argv[1];
    const std::string output_path = argv[2];
    const std::string key = argc > 3 ? argv[3] : "sbnd";
    const Long64_t max_entries = argc > 4 ? std::stoll(argv[4]) : -1;

    if (input_path == output_path)
    {
        std::cerr << "ERROR: refusing to write the output over the input (" << input_path << ")" << std::endl;
        return 1;
    }

    reweight::copy_sample(input_path, output_path, key, max_entries);

    return 0;
}
