
/**
 * @file
 * @author  Aapo Kyrola <akyrola@cs.cmu.edu>
 * @version 1.0
 *
 * @section LICENSE
 *
 * Copyright [2012] [Aapo Kyrola, Guy Blelloch, Carlos Guestrin / Carnegie Mellon University]
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.

 *
 * @section DESCRIPTION
 *
 * Template for GraphChi applications. To create a new application, duplicate
 * this template.
 */

#include <string>
#include <climits>
#include "graphchi_basic_includes.hpp"

using namespace graphchi;

const unsigned INF = UINT32_MAX;

/**
 * Type definitions. Remember to create suitable graph shards using the
 * Sharder-program.
 */
typedef uint VertexDataType;
typedef uint EdgeDataType;

/**
 * SSSP
 */
struct SspProgram : public GraphChiProgram<VertexDataType, EdgeDataType>
{

    /**
     *  Vertex update function.
     */
    void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &gcontext)
    {

        if (gcontext.iteration == 0)
        {
            /* On first iteration, initialize vertex (and its edges). This is usually required, because
               on each run, GraphChi will modify the data files. To start from scratch, it is easiest
               do initialize the program in code. Alternatively, you can keep a copy of initial data files. */
            if (vertex.id() == 0)
            {
                vertex.set_data(1);
                for (int i = 0; i < vertex.num_outedges(); i++)
                {
                    graphchi_edge<uint> *edge = vertex.outedge(i);
                    edge->set_data(1);
                }
            }
        }
        else
        {
            /* Do computation */
            /* Loop over in-edges (example) */
            for (int i = 0; i < vertex.num_inedges(); i++)
            {
                // Do something
                //    value += vertex.inedge(i).get_data();
                graphchi_edge<uint> *edge = vertex.inedge(i);
                if (edge->get_data() != 0)
                {
                    vertex.set_data(1);
                }
            }

            /* Loop over out-edges (example) */
            for (int i = 0; i < vertex.num_outedges(); i++)
            {
                // Do something
                // vertex.outedge(i).set_data(x)
                vertex.outedge(i)->set_data(1);
            }

            /* Loop over all edges (ignore direction) */
            // for(int i=0; i < vertex.num_edges(); i++) {
            //     // vertex.edge(i).get_data()
            // }

            // v.set_data(new_value);
        }
    }

    /**
     * Called before an iteration starts.
     */
    void before_iteration(int iteration, graphchi_context &gcontext)
    {
    }

    /**
     * Called after an iteration has finished.
     */
    void after_iteration(int iteration, graphchi_context &gcontext)
    {
    }

    /**
     * Called before an execution interval is started.
     */
    void before_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &gcontext)
    {
    }

    /**
     * Called after an execution interval has finished.
     */
    void after_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &gcontext)
    {
    }
};

int main(int argc, const char **argv)
{
    /* GraphChi initialization will read the command line
       arguments and the configuration file. */
    graphchi_init(argc, argv);
    printf("INF number: %llu size: %d", INF, sizeof(INF));
    /* Metrics object for keeping track of performance counters
       and other information. Currently required. */
    metrics m("ssp");
    global_logger().set_log_level(LOG_DEBUG);

    /* Basic arguments for application */
    std::string filename = get_option_string("file"); // Base filename
    int niters = get_option_int("niters", 20);         // Number of iterations
    bool scheduler = get_option_int("scheduler", 0);  // Whether to use selective scheduling

    /* Detect the number of shards or preprocess an input to create them */
    int nshards = convert_if_notexists<EdgeDataType>(filename,
                                                     get_option_string("nshards", "auto"));

    /* Run */
    SspProgram program;
    graphchi_engine<VertexDataType, EdgeDataType> engine(filename, nshards, scheduler, m);
    engine.run(program, niters);

    /* Report execution metrics */
    metrics_report(m);
    return 0;
}
