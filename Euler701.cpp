#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <cstdint>
#include <pthread.h>
#include <unistd.h>
#include <vector>

#define DEFAULT_WIDTH 7
#define DEFAULT_HEIGHT 7
#define MAX_CONNECTED    5

#define USE_128BIT_WEIGHT
#if defined(USE_128BIT_WEIGHT)
    using ProblemWeight = __int128 ;
    #define MAX_COLUMN_LIMIT    9
#else
    using ProblemWeight = uint64_t ;
    #define MAX_COLUMN_LIMIT    7
#endif


typedef uint32_t IndexValue  ;
typedef uint8_t SurfaceValue ;
typedef struct CombinationIndex {
    IndexValue N ;
    IndexValue k ;
    IndexValue *index ;
} CombinationIndex ;

CombinationIndex * combinationIndexAlloc(IndexValue N, IndexValue k) ;
CombinationIndex * combinationIndexDestroy(CombinationIndex * NS) ;
IndexValue combinationIndexSize(const CombinationIndex *NS,int ks) ;
IndexValue combinationIndexRank(const CombinationIndex *NS,int ks,SurfaceValue *sum) ;


typedef struct PartitionLayout {
    int id ;
    int source_mask ;
    int set_count ;
    uint16_t part_masks[MAX_CONNECTED]  ;
} PartitionLayout ;


typedef struct PartitionTransition {
    int     next_partition ;
    uint8_t source_index[MAX_CONNECTED] ;
    uint8_t added_surface[MAX_CONNECTED] ;
} PartitionTransition ;


typedef struct DpState {
    uint16_t partition_id ;
    SurfaceValue  surface[MAX_CONNECTED+1] ;
    ProblemWeight  prob ;
} DpState ;

typedef struct GeneratedState {
    IndexValue hash_key ;
    uint16_t partition_id ;
    SurfaceValue surface[MAX_CONNECTED+1] ;
    ProblemWeight prob ;
} GeneratedState ;

typedef struct BuildRangeTask {
    const DpState *states ;
    int state_start ;
    int state_end ;
    const PartitionLayout *partitions ;
    const PartitionTransition *transitions ;
    const CombinationIndex *hash_basis ;
    const IndexValue *partition_base_index ;
    int mask_count ;
    int added_mask ;
    int duplicate_factor ;
    int hash_table_size ;
    std::vector<GeneratedState> locals ;
} BuildRangeTask ;

typedef struct LayerSumTask {
    const DpState *states ;
    const PartitionLayout *partitions ;
    int state_start ;
    int state_end ;
    ProblemWeight partial ;
} LayerSumTask ;

typedef struct SymmetrySumTask {
    const DpState *states ;
    const PartitionLayout *partitions ;
    const PartitionTransition *transitions ;
    int mask_count ;
    int state_start ;
    int state_end ;
    ProblemWeight partial ;
} SymmetrySumTask ;

typedef struct PrecomputeData {
    int mask_count ;
    int partition_count ;
    int *symmetry ;
    PartitionLayout *partitions ;
    int *first_partition_by_mask ;
    PartitionTransition *transitions ;
} PrecomputeData ;

void buildPrecompute(PrecomputeData *solver_data, int width) ;

static int get_thread_count() {
    long n = sysconf(_SC_NPROCESSORS_ONLN);
    if(n < 1) return 1;
    return (int)n;
}

static void *buildTransitionWorker(void *arg) {
    BuildRangeTask *task = (BuildRangeTask *)arg;
    if(task->state_start >= task->state_end) {
        return nullptr;
    }
    task->locals.clear();
    task->locals.reserve((size_t)(task->state_end - task->state_start));
    std::vector<IndexValue> localHash(task->hash_table_size, 0);
    ProblemWeight addMul = task->duplicate_factor ? 2 : 1;

    for(int state_index = task->state_start; state_index < task->state_end; ++state_index) {
        const DpState *source_state = task->states + state_index;
        int source_partition = source_state->partition_id ;
        const PartitionTransition *trf = task->transitions + source_partition * task->mask_count + task->added_mask ;
        int target_partition = trf->next_partition ;
        SurfaceValue surf[MAX_CONNECTED+1] ;
        surf[0] = source_state->surface[0] ;
        for(int segment = 0; segment < task->partitions[target_partition].set_count; segment++) {
            surf[segment+1] = trf->added_surface[segment] ;
        }
        for(int segment = 0; segment < task->partitions[source_partition].set_count; segment++) {
            if(trf->source_index[segment]) {
                surf[trf->source_index[segment]] += source_state->surface[segment+1];
            } else {
                if(source_state->surface[segment+1] > surf[0]) surf[0] = source_state->surface[segment+1] ;
            }
        }
        IndexValue numHash = task->partition_base_index[target_partition] + combinationIndexRank(task->hash_basis, task->partitions[target_partition].set_count+1, surf) ;
        IndexValue loc = localHash[numHash] ;
        ProblemWeight add = addMul * source_state->prob ;
        if(loc == 0) {
            GeneratedState state ;
            state.hash_key = numHash ;
            state.partition_id = (uint16_t)target_partition ;
            memcpy(state.surface, surf, sizeof(state.surface)) ;
            state.prob = add ;
            task->locals.push_back(state) ;
            localHash[numHash] = (IndexValue) task->locals.size() ;
        } else {
            task->locals[loc - 1].prob += add ;
        }
    }
    return nullptr;
}

static void *accumulateLayerWorker(void *arg) {
    LayerSumTask *task = (LayerSumTask *)arg;
    ProblemWeight total = 0 ;
    for(int state_index = task->state_start; state_index < task->state_end; state_index++) {
        const DpState *new_state = task->states + state_index ;
        int maxS = new_state->surface[0];
        const PartitionLayout *pp = task->partitions + new_state->partition_id ;
        for(int segment = 0; segment < pp->set_count; segment++) {
            if(new_state->surface[segment+1] > maxS) maxS = new_state->surface[segment+1];
        }
        total += (ProblemWeight)maxS * new_state->prob ;
    }
    task->partial = total;
    return nullptr;
}

static void *accumulateSymmetryWorker(void *arg) {
    SymmetrySumTask *task = (SymmetrySumTask *)arg;
    ProblemWeight total = 0 ;
    for(int state_index = task->state_start; state_index < task->state_end; state_index++) {
        const DpState *source_state = task->states + state_index ;
        int source_partition = source_state->partition_id ;
        ProblemWeight state_prob = source_state->prob ;
        for(int column = 0; column < task->mask_count; column++) {
            const PartitionTransition *trf = task->transitions + source_partition*task->mask_count + column ;
            int target_partition = trf->next_partition ;
            SurfaceValue surf[MAX_CONNECTED] ;
            int maxS = source_state->surface[0] ;
            for(int segment = 0; segment < task->partitions[target_partition].set_count; segment++) {
                surf[segment] = trf->added_surface[segment] ;
            }
            for(int segment = 0; segment < task->partitions[source_partition].set_count; segment++) {
                if(trf->source_index[segment]) {
                    surf[trf->source_index[segment]-1] += source_state->surface[segment+1];
                } else {
                    if(source_state->surface[segment+1] > maxS ) maxS = source_state->surface[segment+1] ;
                }
            }
            for(int segment = 0; segment < task->partitions[target_partition].set_count; segment++) {
                if(surf[segment] > maxS) maxS = surf[segment]  ;
            }
            total += (ProblemWeight)maxS * state_prob;
        }
    }
    task->partial = total;
    return nullptr;
}

int main(int argc, const char * argv[]) {
    PrecomputeData solver{};
    int width = DEFAULT_WIDTH ;
    int maxNbColumn =  DEFAULT_HEIGHT ;
    if(argc > 1){
        width = atoi(argv[1]) ;
        if(width > MAX_COLUMN_LIMIT) width = MAX_COLUMN_LIMIT ;
        if(argc > 2) {
            maxNbColumn = atoi(argv[2]);
            if(maxNbColumn > MAX_COLUMN_LIMIT) maxNbColumn = MAX_COLUMN_LIMIT ;
            if(maxNbColumn < width) {
                int tmp = width ;
                width = maxNbColumn ;
                maxNbColumn = tmp ;
            }
        }
    }
    buildPrecompute(&solver, width) ;
    int maskCount = solver.mask_count ;
    const PartitionLayout *partitions = solver.partitions ;
    const PartitionTransition *transitions = solver.transitions;

    DpState *states = (DpState*) malloc(300000000 * sizeof(states[0])) ;
    int stateCount = 1;
    IndexValue *partitionBaseIndex = (IndexValue*) calloc(solver.partition_count, sizeof(partitionBaseIndex[0])) ;
    int firstStateInPrevLayer = 0 ;
    memset(states, 0, sizeof(states[0]));
    states[0].prob = 1 ;
    states[0].partition_id = 0 ;
    int firstStateInCurrentLayer = 1;
    double result ;
    for(int columns = 1; columns < maxNbColumn; columns++) {
        int totalWeight = width * columns ;
        CombinationIndex *hashIndex = combinationIndexAlloc(totalWeight, MAX_CONNECTED+1) ;
        ProblemWeight sumCurrent = 0 , sumSymmetric = 0 ;
        for(int addedMask = 0; addedMask < maskCount; addedMask++) {
            int symmetryMultiplier = 0 ;
            if(columns == maxNbColumn - 1) {
                if(solver.symmetry[addedMask] < 0) continue ;
                symmetryMultiplier = solver.symmetry[addedMask] ? 1 : 0 ;
            }
            int hashSize = 0 ;
            int firstStateForMask = stateCount ;
            for(int partitionId = solver.first_partition_by_mask[addedMask]; partitionId < solver.first_partition_by_mask[addedMask + 1]; partitionId++) {
                partitionBaseIndex[partitionId] = hashSize ;
                hashSize += combinationIndexSize(hashIndex, partitions[partitionId].set_count + 1)  ;
            }

            std::vector<pthread_t> buildThreads;
            std::vector<BuildRangeTask> buildTasks;
            int activeStateCount = firstStateInCurrentLayer - firstStateInPrevLayer ;
            int buildThreadCount = get_thread_count();
            if(activeStateCount < 4096) buildThreadCount = 1 ;
            if(buildThreadCount > activeStateCount) buildThreadCount = activeStateCount;
            if(buildThreadCount <= 0) buildThreadCount = 1;
            buildThreads.resize(buildThreadCount);
            buildTasks.resize(buildThreadCount);
            int block = (activeStateCount + buildThreadCount -1) / buildThreadCount ;
            for(int thread = 0; thread < buildThreadCount; thread++) {
                int start = firstStateInPrevLayer + thread * block ;
                int end = start + block ;
                if(end > firstStateInCurrentLayer) end = firstStateInCurrentLayer;
                BuildRangeTask &task = buildTasks[thread];
                task.states = states;
                task.state_start = start;
                task.state_end = end;
                task.partitions = partitions;
                task.transitions = transitions;
                task.hash_basis = hashIndex;
                task.partition_base_index = partitionBaseIndex;
                task.mask_count = maskCount;
                task.added_mask = addedMask;
                task.duplicate_factor = symmetryMultiplier;
                task.hash_table_size = hashSize;
                pthread_create(&buildThreads[thread], NULL, buildTransitionWorker, &task);
            }
            for(int thread = 0; thread < buildThreadCount; thread++) pthread_join(buildThreads[thread], NULL);

            std::vector<IndexValue> stateSlot(hashSize, 0);
            for(int thread = 0; thread < buildThreadCount; thread++) {
                const std::vector<GeneratedState> &local = buildTasks[thread].locals;
                for(size_t i = 0; i < local.size(); ++i) {
                    const GeneratedState &generated = local[i];
                    IndexValue existing = stateSlot[generated.hash_key];
                    if(existing == 0) {
                        stateSlot[generated.hash_key] = (IndexValue)(stateCount + 1);
                        DpState *nextState = states + stateCount++;
                        nextState->partition_id = generated.partition_id ;
                        nextState->prob = generated.prob ;
                        memcpy(nextState->surface, generated.surface, sizeof(nextState->surface));
                    } else {
                        states[existing - 1].prob += generated.prob ;
                    }
                }
            }

            int statesInSlice = stateCount - firstStateForMask ;
            if(statesInSlice > 0) {
                int sumThreads = get_thread_count();
                if(statesInSlice < 2048) sumThreads = 1 ;
                if(sumThreads > statesInSlice) sumThreads = statesInSlice;
                if(sumThreads <= 0) sumThreads = 1;
                std::vector<pthread_t> sumWorkers(sumThreads);
                std::vector<LayerSumTask> sumTasks(sumThreads);
                int per = (statesInSlice + sumThreads -1) / sumThreads ;
                for(int thread = 0; thread < sumThreads; thread++) {
                    int start = firstStateForMask + thread * per ;
                    int end = start + per ;
                    if(end > stateCount) end = stateCount;
                    LayerSumTask &task = sumTasks[thread];
                    task.states = states;
                    task.partitions = partitions;
                    task.state_start = start;
                    task.state_end = end;
                    task.partial = 0;
                    pthread_create(&sumWorkers[thread], NULL, accumulateLayerWorker, &task);
                }
                for(int thread = 0; thread < sumThreads; thread++) {
                    pthread_join(sumWorkers[thread], NULL);
                    sumCurrent += sumTasks[thread].partial ;
                }
            }

            if(columns == maxNbColumn - 1) {
                if(statesInSlice > 0) {
                    int symThreads = get_thread_count();
                    if(symThreads < 1) symThreads = 1 ;
                    if(symThreads > statesInSlice) symThreads = statesInSlice;
                    if(statesInSlice < 256) symThreads = 1 ;
                    std::vector<pthread_t> symWorkers(symThreads);
                    std::vector<SymmetrySumTask> symTasks(symThreads);
                    int per = (statesInSlice + symThreads -1) / symThreads ;
                    for(int thread = 0; thread < symThreads; thread++) {
                        int start = firstStateForMask + thread * per ;
                        int end = start + per ;
                        if(end > stateCount) end = stateCount;
                        SymmetrySumTask &task = symTasks[thread];
                        task.states = states;
                        task.partitions = partitions;
                        task.transitions = transitions;
                        task.mask_count = maskCount;
                        task.state_start = start;
                        task.state_end = end;
                        task.partial = 0;
                        pthread_create(&symWorkers[thread], NULL, accumulateSymmetryWorker, &task);
                    }
                    for(int thread = 0; thread < symThreads; thread++) {
                        pthread_join(symWorkers[thread], NULL);
                        sumSymmetric += symTasks[thread].partial ;
                    }
                }
                stateCount = firstStateForMask ;
            }
        }

        firstStateInPrevLayer = firstStateInCurrentLayer ;
        firstStateInCurrentLayer = stateCount ;
        hashIndex = combinationIndexDestroy(hashIndex);
        int exp2 = columns * width ;
        if(columns == maxNbColumn-1) {
            exp2 = (columns + 1) * width ;
            result = (double)sumSymmetric / (double) (1LL << (exp2/2 )) / (double)(1LL << ((exp2+1)/2)) ;
            printf("%.8f\n", result);
        }
    }
    free(partitionBaseIndex);
    free(states);
    return 0 ;
}

typedef struct MaskTransition {
    uint8_t surface_increment ;
    uint16_t contact_mask ;
} MaskTransition ;

typedef struct MaskProfile {
    int mask_union ;
    int set_count ;
    uint16_t part_masks[MAX_CONNECTED]  ;
} MaskProfile ;

void buildPrecompute(PrecomputeData *solver_data, int width) {
    solver_data->mask_count = 1 << width ;
    int *bit_count_by_mask = (int*) malloc(solver_data->mask_count * sizeof(bit_count_by_mask[0]));
    for(int mask = 0; mask < solver_data->mask_count; ++mask) {
        int bit_count = 0 ;
        for(int bit = 0; bit < width; ++bit) {
            if(mask & (1 << bit)) ++bit_count ;
        }
        bit_count_by_mask[mask] = bit_count ;
    }

    solver_data->symmetry = (int*) malloc(solver_data->mask_count * sizeof(solver_data->symmetry[0]));
    for(int mask = 0; mask < solver_data->mask_count; ++mask) {
        int mirror = 0 ;
        for(int bit = 0; bit < width; ++bit) {
            mirror |= ((mask >> bit) & 1) << (width - bit - 1) ;
        }
        if(mirror == mask) solver_data->symmetry[mask] = 0 ;
        else if(mask > mirror) solver_data->symmetry[mask] = 1;
        else solver_data->symmetry[mask] = -1;
    }

    MaskTransition *intersections = (MaskTransition*) malloc(solver_data->mask_count * solver_data->mask_count * sizeof(intersections[0]));
    for(int prev_mask = 0; prev_mask < solver_data->mask_count; ++prev_mask) {
        for(int cur_mask = 0; cur_mask < solver_data->mask_count; ++cur_mask) {
            int km = prev_mask & cur_mask ;
            if(km==0) {
                intersections[prev_mask*solver_data->mask_count + cur_mask].contact_mask = 0 ;
                intersections[prev_mask*solver_data->mask_count + cur_mask].surface_increment = 0 ;
                continue ;
            }
            for(int m = 0; m < width; m++) {
                if((m>0) && (km & (1 << (m - 1))) && (cur_mask & (1 << m)) ) km |= 1 << m ;
                if((m<width-1) && (km & (1 << (m + 1))) && (cur_mask & (1 << m)) ) km |= 1 << m ;
            }
            for(int m = width - 1; m >= 0; m--) {
                if((m>0) && (km & (1 << (m - 1))) && (cur_mask & (1 << m)) ) km |= 1 << m ;
                if((m<width-1) && (km & (1 << (m + 1))) && (cur_mask & (1 << m)) ) km |= 1 << m ;
            }
            intersections[prev_mask*solver_data->mask_count + cur_mask].surface_increment = (uint8_t) bit_count_by_mask[km] ;
            intersections[prev_mask*solver_data->mask_count + cur_mask].contact_mask = (uint16_t) km ;
        }
    }

    MaskProfile *profiles = (MaskProfile*) malloc(solver_data->mask_count * sizeof(profiles[0])) ;
    int partition_slots = 0 ;
    for(int mask = 0; mask < solver_data->mask_count; ++mask) {
        uint16_t component_masks[MAX_CONNECTED] ;
        memset(component_masks,0,sizeof(component_masks));
        int component_count = 0 ;
        int run = mask & 1 ;
        for(int bit = 1; bit <= width; ++bit) {
            if(mask & (1 << bit)) {
                ++run ;
            } else if(run) {
                component_masks[component_count++] = (uint16_t)((1 << bit) - (1 << (bit - run))) ;
                run = 0 ;
            }
        }
        assert(component_count <= width);
        int partitions_by_components[MAX_CONNECTED+1] = { 1,1,2,5,14,42} ;
        partition_slots += partitions_by_components[component_count] ;
        profiles[mask].mask_union = mask ;
        profiles[mask].set_count = component_count ;
        memcpy(profiles[mask].part_masks, component_masks, sizeof(component_masks)) ;
    }

    solver_data->partitions = (PartitionLayout*) malloc(partition_slots * sizeof(solver_data->partitions[0]));
    solver_data->first_partition_by_mask = (int*) malloc((solver_data->mask_count + 1) * sizeof(solver_data->first_partition_by_mask[0])) ;
    solver_data->partition_count = 0 ;
    solver_data->first_partition_by_mask[0] = 0 ;
    for(int mask = 0; mask < solver_data->mask_count; ++mask) {
        uint16_t component_masks[MAX_CONNECTED] ;
        memcpy(component_masks, profiles[mask].part_masks, sizeof(component_masks));
        int component_count = profiles[mask].set_count ;
        int first = solver_data->partition_count ;
        PartitionLayout *slot = nullptr ;

#define INIT_PARTITION(count) do { \
            slot = solver_data->partitions + solver_data->partition_count; \
            slot->id = solver_data->partition_count; \
            slot->source_mask = mask; \
            ++solver_data->partition_count; \
            memset(slot->part_masks,0,sizeof(slot->part_masks)); \
            slot->set_count = (count); \
        } while(0)

#define DST(i) slot->part_masks[i]
#define SRC(i) component_masks[i]
#define PART_A(i0)         INIT_PARTITION(1); DST(0)=SRC(i0)
#define PART_AB(i0,i1)     INIT_PARTITION(1); DST(0)=SRC(i0)|SRC(i1)
#define PART_A_B(i0,i1)    INIT_PARTITION(2); DST(0)=SRC(i0); DST(1)=SRC(i1)
#define PART_ABC(i0,i1,i2)     INIT_PARTITION(1); DST(0)=SRC(i0)|SRC(i1)|SRC(i2)
#define PART_AB_C(i0,i1,i2)    INIT_PARTITION(2); DST(0)=SRC(i0)|SRC(i1); DST(1)=SRC(i2)
#define PART_A_B_C(i0,i1,i2)   INIT_PARTITION(3); DST(0)=SRC(i0); DST(1)=SRC(i1); DST(2)=SRC(i2)
#define PART_ABCD(i0,i1,i2,i3)     INIT_PARTITION(1); DST(0)=SRC(i0)|SRC(i1)|SRC(i2)|SRC(i3)
#define PART_AB_CD(i0,i1,i2,i3)    INIT_PARTITION(2); DST(0)=SRC(i0)|SRC(i1); DST(1)=SRC(i2)|SRC(i3)
#define PART_ABC_D(i0,i1,i2,i3)    INIT_PARTITION(2); DST(0)=SRC(i0)|SRC(i1)|SRC(i2); DST(1)=SRC(i3)
#define PART_AB_C_D(i0,i1,i2,i3)   INIT_PARTITION(3); DST(0)=SRC(i0)|SRC(i1); DST(1)=SRC(i2); DST(2)=SRC(i3)
#define PART_A_B_C_D(i0,i1,i2,i3)  INIT_PARTITION(4); DST(0)=SRC(i0); DST(1)=SRC(i1); DST(2)=SRC(i2); DST(3)=SRC(i3)
#define PART_ABCDE(i0,i1,i2,i3,i4)     INIT_PARTITION(1); DST(0)=SRC(i0)|SRC(i1)|SRC(i2)|SRC(i3)|SRC(i4)
#define PART_ABCD_E(i0,i1,i2,i3,i4)    INIT_PARTITION(2); DST(0)=SRC(i0)|SRC(i1)|SRC(i2)|SRC(i3); DST(1)=SRC(i4)
#define PART_ABC_DE(i0,i1,i2,i3,i4)    INIT_PARTITION(2); DST(0)=SRC(i0)|SRC(i1)|SRC(i2); DST(1)=SRC(i3)|SRC(i4)
#define PART_ABC_D_E(i0,i1,i2,i3,i4)   INIT_PARTITION(3); DST(0)=SRC(i0)|SRC(i1)|SRC(i2); DST(1)=SRC(i3); DST(2)=SRC(i4)
#define PART_AB_CD_E(i0,i1,i2,i3,i4)   INIT_PARTITION(3); DST(0)=SRC(i0)|SRC(i1); DST(1)=SRC(i2)|SRC(i3); DST(2)=SRC(i4)
#define PART_AB_C_D_E(i0,i1,i2,i3,i4)  INIT_PARTITION(4); DST(0)=SRC(i0)|SRC(i1); DST(1)=SRC(i2); DST(2)=SRC(i3); DST(3)=SRC(i4)
#define PART_A_B_C_D_E(i0,i1,i2,i3,i4) INIT_PARTITION(5); DST(0)=SRC(i0); DST(1)=SRC(i1); DST(2)=SRC(i2); DST(3)=SRC(i3); DST(4)=SRC(i4)
        switch(component_count) {
            case 0:
                INIT_PARTITION(0);
                break;
            case 1:
                PART_A(0);
                break;
            case 2:
                PART_AB(0,1);
                PART_A_B(0,1);
                break;
            case 3:
                PART_ABC(0,1,2);
                PART_AB_C(0,1,2); PART_AB_C(0,2,1); PART_AB_C(1,2,0);
                PART_A_B_C(0,1,2);
                break;
            case 4:
                PART_ABCD(0,1,2,3);
                PART_ABC_D(0,1,2,3); PART_ABC_D(0,1,3,2); PART_ABC_D(0,2,3,1); PART_ABC_D(1,2,3,0);
                PART_AB_CD(0,1,2,3); PART_AB_CD(0,3,1,2);
                PART_AB_C_D(0,1,2,3); PART_AB_C_D(0,2,1,3); PART_AB_C_D(0,3,1,2); PART_AB_C_D(1,2,0,3); PART_AB_C_D(1,3,0,2); PART_AB_C_D(2,3,0,1);
                PART_A_B_C_D(0,1,2,3);
                break;
            case 5:
                PART_ABCDE(0,1,2,3,4);
                PART_ABCD_E(0,1,2,3,4); PART_ABCD_E(0,1,2,4,3); PART_ABCD_E(0,1,3,4,2); PART_ABCD_E(0,2,3,4,1); PART_ABCD_E(1,2,3,4,0);
                PART_ABC_DE(0,1,2,3,4); PART_ABC_DE(0,1,4,2,3); PART_ABC_DE(0,3,4,1,2); PART_ABC_DE(2,3,4,0,1); PART_ABC_DE(1,2,3,0,4);
                PART_ABC_D_E(0,1,2,3,4); PART_ABC_D_E(0,1,3,2,4); PART_ABC_D_E(0,2,3,1,4); PART_ABC_D_E(1,2,3,0,4); PART_ABC_D_E(0,1,4,2,3);
                PART_ABC_D_E(0,2,4,1,3); PART_ABC_D_E(1,2,4,0,3); PART_ABC_D_E(0,3,4,1,2); PART_ABC_D_E(1,3,4,0,2); PART_ABC_D_E(2,3,4,0,1);
                PART_AB_CD_E(0,1,2,3,4); PART_AB_CD_E(0,3,1,2,4); PART_AB_CD_E(0,1,2,4,3); PART_AB_CD_E(0,4,1,2,3); PART_AB_CD_E(0,1,3,4,2);
                PART_AB_CD_E(0,4,1,3,2); PART_AB_CD_E(0,2,3,4,1); PART_AB_CD_E(0,4,2,3,1); PART_AB_CD_E(1,2,3,4,0); PART_AB_CD_E(1,4,2,3,0);
                PART_AB_C_D_E(0,1,2,3,4); PART_AB_C_D_E(0,2,1,3,4); PART_AB_C_D_E(0,3,1,2,4); PART_AB_C_D_E(0,4,1,2,3); PART_AB_C_D_E(1,2,0,3,4);
                PART_AB_C_D_E(1,3,0,2,4); PART_AB_C_D_E(1,4,0,2,3); PART_AB_C_D_E(2,3,0,1,4); PART_AB_C_D_E(2,4,0,1,3); PART_AB_C_D_E(3,4,0,1,2);
                PART_A_B_C_D_E(0,1,2,3,4);
                break;
        }
        solver_data->first_partition_by_mask[mask+1] = solver_data->partition_count ;
        for(int ip = first; ip < solver_data->partition_count; ++ip) {
            int count = solver_data->partitions[ip].set_count;
            for(int i = 1; i < count; ++i) {
                for(int j = i; j < count; ++j) {
                    if(solver_data->partitions[ip].part_masks[i-1] > solver_data->partitions[ip].part_masks[j]) {
                        uint16_t tmp = solver_data->partitions[ip].part_masks[i-1] ;
                        solver_data->partitions[ip].part_masks[i-1] = solver_data->partitions[ip].part_masks[j] ;
                        solver_data->partitions[ip].part_masks[j] = tmp ;
                    }
                }
            }
        }
    }

    solver_data->transitions = (PartitionTransition*) calloc(solver_data->partition_count * solver_data->mask_count, sizeof(solver_data->transitions[0]));
    for(int partition_id = 0; partition_id < solver_data->partition_count; ++partition_id) {
        const PartitionLayout *from_partition = solver_data->partitions + partition_id ;
        for(int mask = 0; mask < solver_data->mask_count; ++mask) {
            int nbSet = 0 ;
            uint16_t merged[MAX_CONNECTED] ;
            memset(merged,0,sizeof(merged));
            int global_contact = 0 ;
            for(int i = 0; i < from_partition->set_count; ++i) {
                uint16_t new_contact = intersections[from_partition->part_masks[i] * solver_data->mask_count + mask].contact_mask ;
                if(new_contact) {
                    for(int j = 0; j < nbSet; ++j) {
                        if((new_contact & merged[j])) {
                            new_contact = intersections[(new_contact | merged[j]) * solver_data->mask_count + mask].contact_mask ;
                            for(int k = j + 1; k < nbSet; ++k) merged[k - 1] = merged[k] ;
                            --nbSet ;
                            --j ;
                            merged[nbSet] = 0 ;
                        }
                    }
                    int pos = nbSet - 1 ;
                    while(pos >= 0 && merged[pos] > new_contact) {
                        merged[pos + 1] = merged[pos] ;
                        --pos ;
                    }
                    merged[pos + 1] = new_contact ;
                    ++nbSet ;
                }
                global_contact |= new_contact ;
            }
            if(global_contact != mask) {
                global_contact ^= mask ;
                for(int k = 0; k < profiles[mask].set_count; ++k) {
                    if((profiles[mask].part_masks[k] & global_contact) == profiles[mask].part_masks[k] ) {
                        int pos = nbSet - 1 ;
                        while(pos >= 0 && merged[pos] > profiles[mask].part_masks[k]) {
                            merged[pos + 1] = merged[pos] ;
                            --pos ;
                        }
                        merged[pos + 1] = profiles[mask].part_masks[k] ;
                        ++nbSet ;
                    }
                }
            }

            int target_partition = solver_data->first_partition_by_mask[mask] ;
            while(target_partition < solver_data->first_partition_by_mask[mask + 1] && memcmp(merged, solver_data->partitions[target_partition].part_masks, sizeof(merged)) != 0) {
                ++target_partition ;
            }
            assert(target_partition != solver_data->first_partition_by_mask[mask + 1]) ;
            PartitionTransition *transition = solver_data->transitions + partition_id * solver_data->mask_count + mask ;
            transition->next_partition = (int)target_partition ;
            for(int id = 0; id < solver_data->partitions[target_partition].set_count; ++id) {
                transition->added_surface[id] = (uint8_t) bit_count_by_mask[solver_data->partitions[target_partition].part_masks[id]] ;
            }
            for(int i = 0; i < from_partition->set_count; ++i) {
                if(from_partition->part_masks[i] & mask) {
                    for(int id = 0; id < solver_data->partitions[target_partition].set_count; ++id){
                        if(from_partition->part_masks[i] & solver_data->partitions[target_partition].part_masks[id]) {
                            transition->source_index[i] = (uint8_t)(id + 1) ;
                            break ;
                        }
                        assert(id < solver_data->partitions[target_partition].set_count) ;
                    }
                }
            }
        }
    }
    free(profiles);
    free(bit_count_by_mask);
    free(intersections);
}

CombinationIndex * combinationIndexAlloc(IndexValue N, IndexValue k) {
    if(k<2) return NULL ;
    CombinationIndex * NS = (CombinationIndex*) calloc(1, sizeof(NS[0])) ;
    NS->N = N+1 ;
    NS->k = k ;
    NS->index = (IndexValue*) malloc((k-1) * (N+1) * sizeof(NS->index[0]));
    IndexValue *ids = NS->index+ ((k-2)*(N+1) );
    ids[0] = 0 ;
    for(IndexValue in=1;in<=N;in++){
        ids[in] = ids[in-1] + N + 2 - in ;
    }
    for(int ik=k-3;ik>=0;ik--) {
        IndexValue *idsAnt = ids ;
        ids = NS->index+ (ik*(N+1));
        ids[0] = 0 ;
        IndexValue indAnt = idsAnt[N]+1 ;
    for(IndexValue in=1;in<=N;in++){
        ids[in] = ids[in-1] + indAnt - idsAnt[in-1] ;
    }
    }
    return NS ;
}

CombinationIndex * combinationIndexDestroy(CombinationIndex * NS) {
    if(NS != NULL) {
        free(NS->index);
        free(NS);
    }
    return NULL ;
}
IndexValue combinationIndexSize(const CombinationIndex *NS,int ks) {
    if(ks<=1) {
        return ks ? NS->N : 1 ;
    }
    return NS->index[(NS->k-ks)*(NS->N) + (NS->N - 1)] + 1 ;
}

IndexValue combinationIndexRank(const CombinationIndex *NS,int ks,SurfaceValue *sum) {
    int k = NS->k ;
    int N = NS->N ;
    IndexValue * index = NS->index + (k-ks)*N ;
    switch(ks) {
        case 0: return 0;
        case 1: return sum[0] ;
        case 2: return index[sum[0]]+sum[1] ;
        case 3: return index[sum[0]] + index[N+sum[0]+sum[1]]-index[N+sum[0]]+sum[2];
        case 4: return index[sum[0]] + index[N+sum[0]+sum[1]]-index[N+sum[0]]
            + index[2*N+sum[0]+sum[1]+sum[2]]-index[2*N+sum[0]+sum[1]]  +sum[3];
        default : {
            IndexValue ksum = *sum++ ;
            IndexValue ind = index[ksum] ;
            for(int ik=1;ik< ks-1;ik++) {
                IndexValue ksumNext = ksum + *sum++ ;
                ind += index[N*ik+ksumNext]-index[N*ik+ksum];
                ksum=ksumNext ;
            }
            return ind+ *sum ;
        }
    }
}
