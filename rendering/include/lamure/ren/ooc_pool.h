// Copyright (c) 2014 Bauhaus-Universitaet Weimar
// This Software is distributed under the Modified BSD License, see license.txt.
//
// Virtual Reality and Visualization Research Group 
// Faculty of Media, Bauhaus-Universitaet Weimar
// http://www.uni-weimar.de/medien/vr

#ifndef REN_OOC_POOL_H_
#define REN_OOC_POOL_H_

#include <vector>
#include <map>
#include <thread>
#include <mutex>
#include <set>
#include <lamure/semaphore.h>
#include <lamure/utils.h>
#include <lamure/types.h>
#include <lamure/ren/config.h>
#include <lamure/ren/model_database.h>
#include <lamure/ren/lod_stream.h>
#include <lamure/ren/cache_queue.h>
#include <lamure/ren/cache_index.h>



namespace lamure {
namespace ren {

class ooc_pool
{
public:
                        ooc_pool(const uint32_t num_loader_threads,
                                const size_t size_of_slot_in_bytes);
    /*virtual*/         ~ooc_pool();

    const uint32_t      num_threads() const { return num_threads_; };

    bool                acknowledge_request(cache_queue::job job);
    void                acknowledge_update(const model_t model_id,
                                          const node_t node_id,
                                          int32_t priority);

    cache_queue::query_result acknowledge_query(const model_t model_id, const node_t node_id);

    void                resolve_cache_history(cache_index* index);
    void                perform_queue_maintenance(cache_index* index);

    void                lock();
    void                unlock();

    void                begin_measure();
    void                end_measure();

protected:

    void                run();
    bool                is_shutdown();

private:

    bool                locked_;
    semaphore           semaphore_;
    size_t              size_of_slot_;
    std::mutex          mutex_;

    uint32_t            num_threads_;
    std::vector<std::thread> threads_;

    bool                shutdown_;

    size_t              bytes_loaded_;

    std::vector<cache_queue::job> history_;

    cache_queue          priority_queue_;
};

} } // namespace lamure

#endif // REN_OOC_POOL_H_

