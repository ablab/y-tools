//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "stage.hpp"
//#include "graphio.hpp"

#include "logger/log_writers.hpp"

#include <algorithm>
#include <cstring>

namespace spades {

void AssemblyStage::load(const std::string &load_from,
                         const char* prefix) {
    std::string p = path::append_path(load_from, prefix == NULL ? id_ : prefix);
    INFO("Loading current state from " << p);

//    debruijn_graph::graphio::ScanAll(p, gp, false);
//    debruijn_graph::load_lib_data(p);
}

void AssemblyStage::save(const std::string &save_to,
                         const char* prefix) const {
    std::string p = path::append_path(save_to, prefix == NULL ? id_ : prefix);
    INFO("Saving current state to " << p);
}

class StageIdComparator {
  public:
    StageIdComparator(const char* id)
            : id_(id) {
        const char* pos = strstr(id, ":");
        len_ = (pos != NULL ? pos - id : strlen(id));
    }

    bool operator()(const std::unique_ptr<AssemblyStage> &stage) const {
        const char* sid = stage->id();
        return (0 == strncmp(id_, sid, len_) && sid[len_] == 0);
    }

  private:
    const char* id_;
    size_t len_;
};

class PhaseIdComparator {
  public:
    PhaseIdComparator(const char* id) {
        const char* pos = strstr(id, ":");
        VERIFY(pos != NULL);
        id_ = pos + 1;
    }

    bool operator()(const std::unique_ptr<CompositeStageBase::PhaseBase> &phase) const {
        return 0 == strcmp(id_, phase->id());
    }

  private:
    const char* id_;
};

void CompositeStageBase::run(const char* started_from) {
    VERIFY(parent_);
    auto start_phase = phases_.begin();
    if (started_from &&
        strstr(started_from, ":") &&
        started_from == strstr(started_from, id())) {
        start_phase = std::find_if(phases_.begin(), phases_.end(), PhaseIdComparator(started_from));
        if (start_phase == phases_.end()) {
            ERROR("Invalid start stage / phase combination specified: " << started_from);
            exit(-1);
        }
        if (start_phase != phases_.begin()) {
            PhaseBase * prev_phase = std::prev(start_phase)->get();
            std::string composite_id(id());
            composite_id += ":";
            composite_id += prev_phase->id();
            prev_phase->load(parent_->saves_policy().load_from_, composite_id.c_str());
        }
    }

    for (auto et = phases_.end(); start_phase != et; ++start_phase) {
        PhaseBase *phase = start_phase->get();

        INFO("PROCEDURE == " << phase->name());
        phase->run(started_from);

        if (parent_->saves_policy().make_saves_) {
            std::string composite_id(id());
            composite_id += ":";
            composite_id += phase->id();

            phase->save(parent_->saves_policy().save_to_, composite_id.c_str());
        }

    }
}

void StageManager::run(const char* start_from) {
    auto start_stage = stages_.begin();
    if (start_from) {
        start_stage = std::find_if(stages_.begin(), stages_.end(), StageIdComparator(start_from));
        if (start_stage == stages_.end()) {
            ERROR("Invalid start stage specified: " << start_from);
            exit(-1);
        }
        if (start_stage != stages_.begin())
            (*std::prev(start_stage))->load(saves_policy_.load_from_);
    }

    for (; start_stage != stages_.end(); ++start_stage) {
        AssemblyStage *stage = start_stage->get();

        INFO("STAGE == " << stage->name());
        stage->run(start_from);
        if (saves_policy_.make_saves_)
            stage->save(saves_policy_.save_to_);
    }
}

}
