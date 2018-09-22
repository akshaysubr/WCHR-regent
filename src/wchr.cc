/* Copyright 2018 Stanford University
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "wchr.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <map>
#include <vector>

#include "mappers/default_mapper.h"

using namespace Legion;
using namespace Legion::Mapping;

///
/// Mapper
///

static LegionRuntime::Logger::Category log_wchr("wchr");

class WCHRMapper : public DefaultMapper
{
public:
  WCHRMapper(MapperRuntime *rt, Machine machine, Processor local,
             const char *mapper_name,
             std::vector<Processor>* procs_list,
             std::vector<Memory>* sysmems_list,
             std::map<Memory, std::vector<Processor> >* sysmem_local_procs,
             std::map<Processor, Memory>* proc_sysmems,
             std::map<Processor, Memory>* proc_regmems);
  virtual int default_policy_select_garbage_collection_priority(
                                MapperContext ctx, 
                                MappingKind kind, Memory memory, 
                                const PhysicalInstance &instance,
                                bool meets_fill_constraints,bool reduction);
  virtual void map_inline(const MapperContext        ctx,
                          const InlineMapping&       inline_op,
                          const MapInlineInput&      input,
                                MapInlineOutput&     output);
private:
  // std::vector<Processor>& procs_list;
  // std::vector<Memory>& sysmems_list;
  //std::map<Memory, std::vector<Processor> >& sysmem_local_procs;
  //std::map<Processor, Memory>& proc_sysmems;
  // std::map<Processor, Memory>& proc_regmems;
};

WCHRMapper::WCHRMapper(MapperRuntime *rt, Machine machine, Processor local,
                             const char *mapper_name,
                             std::vector<Processor>* _procs_list,
                             std::vector<Memory>* _sysmems_list,
                             std::map<Memory, std::vector<Processor> >* _sysmem_local_procs,
                             std::map<Processor, Memory>* _proc_sysmems,
                             std::map<Processor, Memory>* _proc_regmems)
  : DefaultMapper(rt, machine, local, mapper_name)//,
    // procs_list(*_procs_list),
    // sysmems_list(*_sysmems_list),
    //sysmem_local_procs(*_sysmem_local_procs),
    //proc_sysmems(*_proc_sysmems)// ,
    // proc_regmems(*_proc_regmems)
{
}

int WCHRMapper::default_policy_select_garbage_collection_priority(
                            MapperContext ctx, MappingKind kind,
                            Memory memory, const PhysicalInstance &inst,
                            bool meets_fill_constraints, bool reduction)
{
  return GC_FIRST_PRIORITY;
}

void WCHRMapper::map_inline(const MapperContext        ctx,
                            const InlineMapping&       inline_op,
                            const MapInlineInput&      input,
                                  MapInlineOutput&     output)
{
  if (input.valid_instances.size() != 0)
  {
    DefaultMapper::map_inline(ctx, inline_op, input, output);
    return;
  }
  const RegionRequirement &req = inline_op.requirement;
  Memory target_memory =
    default_policy_select_target_memory(
        ctx, inline_op.parent_task->current_proc, req);
  LayoutConstraintSet creation_constraints;
  std::vector<FieldID> fields;
  default_policy_select_constraint_fields(ctx, req, fields);
  std::vector<DimensionKind> dimension_ordering(4);
  dimension_ordering[0] = DIM_X;
  dimension_ordering[1] = DIM_Y;
  dimension_ordering[2] = DIM_Z;
  dimension_ordering[3] = DIM_F;
  creation_constraints.add_constraint(SpecializedConstraint())
    .add_constraint(MemoryConstraint(target_memory.kind()))
    .add_constraint(FieldConstraint(fields, false/*contiguous*/,
          false/*inorder*/))
    .add_constraint(OrderingConstraint(dimension_ordering, 
          false/*contigous*/));
  creation_constraints.add_constraint(
      FieldConstraint(fields, false/*contig*/, false/*inorder*/));
  output.chosen_instances.resize(output.chosen_instances.size()+1);
  if (!default_make_instance(ctx, target_memory, creation_constraints,
        output.chosen_instances.back(), INLINE_MAPPING,
        true, true/*meets*/, req))
  {
    // If we failed to make it that is bad
    log_wchr.error("WCHR mapper failed allocation for region "
                   "requirement of inline mapping in task %s (UID %lld) "
                   "in memory " IDFMT "for processor " IDFMT ". This "
                   "means the working set of your application is too big "
                   "for the allotted capacity of the given memory under "
                   "the default mapper's mapping scheme. You have three "
                   "choices: ask Realm to allocate more memory, write a "
                   "custom mapper to better manage working sets, or find "
                   "a bigger machine. Good luck!", 
                   inline_op.parent_task->get_task_name(),
                   inline_op.parent_task->get_unique_id(),
                   target_memory.id,
                   inline_op.parent_task->current_proc.id);
    assert(false);
  }
}

static void create_mappers(Machine machine, HighLevelRuntime *runtime, const std::set<Processor> &local_procs)
{
  std::vector<Processor>* procs_list = new std::vector<Processor>();
  std::vector<Memory>* sysmems_list = new std::vector<Memory>();
  std::map<Memory, std::vector<Processor> >* sysmem_local_procs =
    new std::map<Memory, std::vector<Processor> >();
  std::map<Processor, Memory>* proc_sysmems = new std::map<Processor, Memory>();
  std::map<Processor, Memory>* proc_regmems = new std::map<Processor, Memory>();


  std::vector<Machine::ProcessorMemoryAffinity> proc_mem_affinities;
  machine.get_proc_mem_affinity(proc_mem_affinities);

  for (unsigned idx = 0; idx < proc_mem_affinities.size(); ++idx) {
    Machine::ProcessorMemoryAffinity& affinity = proc_mem_affinities[idx];
    if (affinity.p.kind() == Processor::LOC_PROC) {
      if (affinity.m.kind() == Memory::SYSTEM_MEM) {
        (*proc_sysmems)[affinity.p] = affinity.m;
        if (proc_regmems->find(affinity.p) == proc_regmems->end())
          (*proc_regmems)[affinity.p] = affinity.m;
      }
      else if (affinity.m.kind() == Memory::REGDMA_MEM)
        (*proc_regmems)[affinity.p] = affinity.m;
    }
  }

  for (std::map<Processor, Memory>::iterator it = proc_sysmems->begin();
       it != proc_sysmems->end(); ++it) {
    procs_list->push_back(it->first);
    (*sysmem_local_procs)[it->second].push_back(it->first);
  }

  for (std::map<Memory, std::vector<Processor> >::iterator it =
        sysmem_local_procs->begin(); it != sysmem_local_procs->end(); ++it)
    sysmems_list->push_back(it->first);

  for (std::set<Processor>::const_iterator it = local_procs.begin();
        it != local_procs.end(); it++)
  {
    WCHRMapper* mapper = new WCHRMapper(runtime->get_mapper_runtime(),
                                        machine, *it, "wchr_mapper",
                                        procs_list,
                                        sysmems_list,
                                        sysmem_local_procs,
                                        proc_sysmems,
                                        proc_regmems);
    runtime->replace_default_mapper(mapper, *it);
  }
}

void register_mappers()
{
  HighLevelRuntime::add_registration_callback(create_mappers);
}
